"""Module D: EFT interpretation with anomalous HZZ couplings a2 (CP-even) and a3 (CP-odd).

The HZZ vertex in the JHUGen/MCFM convention is
    A(HZZ) ~ a1 mZ^2 eps1.eps2 + a2 f1_{mu nu} f2^{mu nu} + a3 f1_{mu nu} f~2^{mu nu},
with a1 = 1 in the SM normalisation of this MCFM version. a2 and a3 correspond to the dimension-6
operators of the Higgs basis (c_zz and c~_zz up to normalisation). Their amplitudes grow with the
momentum flowing through the vertex, so the off-shell tail is far more sensitive than the on-shell
rate. The width is fixed to its SM value here; the gg->H coupling strength mu_F and the qq->ZZ
normalisation are profiled.

Per bin, the signal and interference are quadratic in each coupling and are obtained from
MCFM runs at a = 0, +0.5, -0.5 (a2) and a = 0, 0.5 (a3, whose interference with the SM-like
amplitudes vanishes after integrating over angles).
"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mplhep as hep  # noqa: E402
import numpy as np  # noqa: E402
from iminuit import Minuit  # noqa: E402

from hwidth.offshell import mcfm_cross_section, poisson_nll, read_mcfm_hist, rebin_mcfm  # noqa: E402
from hwidth.samples import LUMI_FB  # noqa: E402
from run_moduleB import OFF_EDGES, load_inputs, mcfm  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
EFT = ROOT / "generators" / "mcfm" / "eft" / "output"
RES = ROOT / "results"
FIG = RES / "figures" / "mpl"  # quick-look matplotlib figures; publication figures: plots_root.py
A_GRID = np.linspace(-3.0, 3.0, 241)
K_GG = 1.7


def eft_off(point, run):
    return rebin_mcfm(*read_mcfm_hist(next(EFT.glob(f"*_e_{point}_{run}.dat"))), OFF_EDGES)[0]


def eft_on(point):
    return mcfm_cross_section(next(EFT.glob(f"*_e_{point}_on.dat")))[0]


def quadratic(v0, vp, vm, step=0.5):
    """Coefficients (c0, c1, c2) of v(a) = c0 + c1 a + c2 a^2 from values at 0, +step, -step."""
    c1 = (vp - vm) / (2 * step)
    c2 = (vp + vm - 2 * v0) / (2 * step**2)
    return v0, c1, c2


def build_model():
    on, on_obs, off, off_obs = load_inputs()
    fold = 2.0 * LUMI_FB * off.eff  # fb -> events per bin (all 4l channels, reconstruction efficiency)

    I0, SI0 = mcfm("ggH_I")[0], mcfm("ggH_SI")[0]
    S0 = SI0 - I0
    coeffs = {}
    # a2: three points.
    Sp, Ip = eft_off("a2p", "SI") - eft_off("a2p", "I"), eft_off("a2p", "I")
    Sm, Im = eft_off("a2m", "SI") - eft_off("a2m", "I"), eft_off("a2m", "I")
    coeffs["a2"] = {"S": quadratic(S0, Sp, Sm), "I": quadratic(I0, Ip, Im),
                    "on": quadratic(1.0, eft_on("a2p") / eft_on("sm"), eft_on("a2m") / eft_on("sm"))}
    # a3: CP-odd, no linear term.
    S3, I3 = eft_off("a3p", "SI") - eft_off("a3p", "I"), eft_off("a3p", "I")
    coeffs["a3"] = {"S": (S0, 0 * S0, (S3 - S0) / 0.25), "I": (I0, 0 * I0, (I3 - I0) / 0.25),
                    "on": (1.0, 0.0, (eft_on("a3p") / eft_on("sm") - 1.0) / 0.25)}
    return on, on_obs, off, off_obs, fold, coeffs


def nll_factory(on, on_obs, off, off_obs, fold, c, use_off=True):
    def nll(a, mu_f, mu_qq):
        r_on = c["on"][0] + c["on"][1] * a + c["on"][2] * a**2
        val = poisson_nll(on_obs, mu_f * r_on * on.S + on.B)
        if use_off:
            S = (c["S"][0] + c["S"][1] * a + c["S"][2] * a**2) * fold
            I = (c["I"][0] + c["I"][1] * a + c["I"][2] * a**2) * fold
            n = K_GG * (mu_f * S + np.sqrt(mu_f) * I + off.B) + mu_qq * off.Q + off.O
            val += poisson_nll(off_obs, n)
        return val
    return nll


def profile(nll, grid, use_off=True):
    vals = []
    for a in grid:
        m = Minuit(lambda mu_f, mu_qq: nll(a, mu_f, mu_qq), mu_f=1.0, mu_qq=1.0)
        m.errordef = Minuit.LIKELIHOOD
        m.limits["mu_f"] = (0.01, 10)
        m.limits["mu_qq"] = (0.2, 3)
        if not use_off:
            m.fixed["mu_qq"] = True
        m.migrad()
        vals.append(2 * m.fval)
    vals = np.array(vals)
    return vals - vals.min()


def interval(grid, q, level=3.84):
    inside = grid[q <= level]
    return [float(inside.min()), float(inside.max())] if len(inside) else [None, None]


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    on, on_obs, off, off_obs, fold, coeffs = build_model()
    out = {"grid": A_GRID.tolist(), "couplings": {}}

    hep.style.use("ATLAS")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, name, label in ((axes[0], "a2", r"$a_2$ (CP-even)"), (axes[1], "a3", r"$a_3$ (CP-odd)")):
        c = coeffs[name]
        res = {}
        for tag, use_off, data_on, data_off in (("obs_on", False, on_obs, off_obs), ("obs_comb", True, on_obs, off_obs)):
            q = profile(nll_factory(on, data_on, off, data_off, fold, c, use_off), A_GRID, use_off)
            res[tag] = {"scan": q.round(4).tolist(), "int68": interval(A_GRID, q, 1.0), "int95": interval(A_GRID, q)}
        # Expected (SM Asimov).
        n_on = on.S + on.B
        n_off = K_GG * ((c["S"][0] + c["I"][0]) * fold + off.B) + off.Q + off.O
        for tag, use_off in (("exp_on", False), ("exp_comb", True)):
            q = profile(nll_factory(on, n_on, off, n_off, fold, c, use_off), A_GRID, use_off)
            res[tag] = {"scan": q.round(4).tolist(), "int68": interval(A_GRID, q, 1.0), "int95": interval(A_GRID, q)}
        res["on_shell_rate_coeffs"] = [float(x) for x in c["on"]]
        res["offshell_signal_ratio_at_a=0.5"] = float(
            ((c["S"][0] + c["S"][1] * 0.5 + c["S"][2] * 0.25).sum()) / c["S"][0].sum())
        out["couplings"][name] = res

        for tag, color, ls, lab in (("obs_on", "#ff6b6b", "-", "On-shell only (obs.)"),
                                    ("obs_comb", "k", "-", "On + off shell (obs.)"),
                                    ("exp_on", "#ff6b6b", "--", "On-shell only (exp.)"),
                                    ("exp_comb", "k", "--", "On + off shell (exp.)")):
            ax.plot(A_GRID, res[tag]["scan"], color=color, ls=ls, lw=2, label=lab)
        ax.axhline(3.84, color="grey", ls=":", lw=1)
        ax.set_ylim(0, 10)
        ax.set_xlabel(label)
        ax.set_ylabel(r"$-2\,\Delta\ln L$")
        ax.legend(fontsize=11)
        print(name, {k: v["int95"] for k, v in res.items() if isinstance(v, dict)},
              "off-shell S ratio at 0.5:", round(res["offshell_signal_ratio_at_a=0.5"], 2),
              "on-shell coeffs:", np.round(res["on_shell_rate_coeffs"], 3))
    fig.savefig(FIG / "moduleD_eft.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Off-shell signal shape for a few a2 values (for the web page).
    c = coeffs["a2"]
    out["shapes_a2"] = {str(a): ((c["S"][0] + c["S"][1] * a + c["S"][2] * a**2) * fold * K_GG).round(4).tolist()
                        for a in (-0.5, 0.0, 0.5)}
    out["edges"] = OFF_EDGES.tolist()
    (RES / "moduleD.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
