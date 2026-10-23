"""Module C: stability of the width result against the assumptions of the measurement.

Each scenario relaxes or changes one assumption of the baseline fit (module B). For every scenario
the observed and expected 95% CL upper limits on Gamma_H / Gamma_H^SM are computed with
36.6 fb^-1, and the expected limit is projected to 3000 fb^-1 (HL-LHC) where the systematic
effects dominate over the statistical precision.
"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mplhep as hep  # noqa: E402
import numpy as np  # noqa: E402

from hwidth.offshell import Assumptions, WidthFit, asimov  # noqa: E402
from run_moduleB import GAMMA_SM_MEV, R_GRID, interval, load_inputs  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures" / "mpl"  # quick-look matplotlib figures; publication figures: plots_root.py
HL_LUMI = 3000.0

SCENARIOS = [
    ("baseline", "Baseline: K = 1.7 for S, B, I; qq->ZZ norm. free", Assumptions(), "none"),
    ("K_low", "gg K-factor 1.4 instead of 1.7", Assumptions(K_gg=1.4), "k"),
    ("K_high", "gg K-factor 2.0 instead of 1.7", Assumptions(K_gg=2.0), "k"),
    ("K_bkg_30", "Continuum K-factor uncertain by 30%, K_I = sqrt(K_S K_B)",
     Assumptions(K_bkg_unc=0.30, K_int_mode="sqrt"), "k"),
    ("qq_fixed", "qq->ZZ normalisation constrained to MC within 10%", Assumptions(qq_constraint=0.10), "bkg"),
    ("ff_1000", "New physics resolves the ggH loop: form factor, Lambda = 1 TeV",
     Assumptions(form_factor_scale=1000.0), "bsm"),
    ("ff_500", "New physics resolves the ggH loop: form factor, Lambda = 500 GeV",
     Assumptions(form_factor_scale=500.0), "bsm"),
]
# Ratio c = (kappa_g kappa_Z)^2 off shell / (kappa_g kappa_Z)^2 on shell, e.g. from couplings that run
# with energy. Then mu_off = c * mu_on * R, so a fit assuming c = 1 returns c * R instead of R.
COUPLING_RATIOS = [0.8, 0.9, 1.0, 1.1, 1.2]


def limits(on, on_obs, off, off_obs, a, mu_qq):
    """Observed and expected 95% CL upper limits on R for one scenario at the nominal luminosity,
    plus the expected limit at HL-LHC."""
    q_obs = WidthFit(on, on_obs, off, off_obs, a).scan(R_GRID)
    a_on, a_off = asimov(on, off, 1.0, 1.0, mu_qq, a)
    q_exp = WidthFit(on, a_on, off, a_off, a).scan(R_GRID)
    r = HL_LUMI / 36.6
    on_l, off_l = on.scaled(r), off.scaled(r)
    # Asimov data are generated with the *baseline* model so that a scenario's model change
    # shows up as a shift (bias) as well as a change in sensitivity.
    h_on, h_off = asimov(on_l, off_l, 1.0, 1.0, mu_qq, Assumptions())
    q_hl = WidthFit(on_l, h_on, off_l, h_off, a).scan(R_GRID)
    return {
        "obs95": WidthFit.upper_limit(R_GRID, q_obs),
        "exp95": WidthFit.upper_limit(R_GRID, q_exp),
        "hl_best": float(R_GRID[np.argmin(q_hl)]),
        "hl68": interval(R_GRID, q_hl, 1.0),
        "hl95": interval(R_GRID, q_hl, 3.84),
        "scan_obs": q_obs.round(4).tolist(),
    }


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    on, on_obs, off, off_obs = load_inputs()
    best, _ = WidthFit(on, on_obs, off, off_obs, Assumptions()).best_fit()
    mu_qq = best["mu_qq"]

    out = {"scenarios": [], "decoupled": {}, "coupling_bias": [], "indirect": {}}
    for key, label, a, group in SCENARIOS:
        res = limits(on, on_obs, off, off_obs, a, mu_qq)
        out["scenarios"].append({"key": key, "label": label, "group": group, **res})
        print(f"{key:10s} obs95={res['obs95']:.2f} exp95={res['exp95']:.2f} HL best={res['hl_best']:.2f} "
              f"HL68={res['hl68']}")

    # Relax the central assumption: on- and off-shell signal strengths independent.
    a = Assumptions(couple_on_off=False)
    grid = np.linspace(0, 25, 101)
    q = WidthFit(on, on_obs, off, off_obs, a).scan(grid)
    out["decoupled"] = {"mu_off_best": float(grid[np.argmin(q)]), "mu_off_68": interval(grid, q, 1.0),
                        "mu_off_95": interval(grid, q, 3.84), "grid": grid.tolist(), "scan": q.round(4).tolist()}

    # Bias from an off-shell/on-shell coupling ratio c (e.g. running couplings): HL-LHC Asimov.
    r = HL_LUMI / 36.6
    on_l, off_l = on.scaled(r), off.scaled(r)
    for c in COUPLING_RATIOS:
        n_on, n_off = asimov(on_l, off_l, 1.0, c, mu_qq, Assumptions())  # truth: mu_off = c * mu_on, R_true = 1
        q = WidthFit(on_l, n_on, off_l, n_off, Assumptions()).scan(R_GRID)
        out["coupling_bias"].append({"c": c, "R_fitted": float(R_GRID[np.argmin(q)]), "hl68": interval(R_GRID, q, 1.0)})

    # Indirect (coupling-fit) width: with universal kappa and B_BSM, mu_on = kappa^2 (1 - B_BSM) is all that is
    # measured, so Gamma/Gamma_SM = kappa^2 / (1 - B_BSM) is unbounded unless B_BSM = 0 or kappa <= 1 is assumed.
    mu_on = best["mu_on"]
    bbsm = np.linspace(0, 0.9, 46)
    kappa2 = mu_on / (1 - bbsm)
    out["indirect"] = {"mu_on": mu_on, "B_BSM": bbsm.round(3).tolist(),
                       "R_indirect": (kappa2 / (1 - bbsm)).round(4).tolist()}

    # Figure: assumption scan (HL-LHC 68% intervals) next to the observed 95% limits.
    hep.style.use("ATLAS")
    fig, ax = plt.subplots(figsize=(10, 6))
    labels = [s["label"] for s in out["scenarios"]]
    y = np.arange(len(labels))[::-1]
    for yi, s in zip(y, out["scenarios"]):
        lo, hi = s["hl68"]
        ax.plot([lo, hi], [yi, yi], color="#00b4d8", lw=6, solid_capstyle="butt")
        ax.plot(s["hl_best"], yi, "o", color="k")
    ax.axvline(1.0, color="k", lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel(r"Fitted $\Gamma_H/\Gamma_H^{SM}$ (HL-LHC 3000 fb$^{-1}$ expected, 68% CL; truth = 1)")
    ax.set_xlim(0, 3)
    fig.savefig(FIG / "moduleC_assumptions.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    out["gamma_sm_mev"] = GAMMA_SM_MEV
    (RES / "moduleC.json").write_text(json.dumps(out, indent=1))
    print("decoupled", {k: v for k, v in out["decoupled"].items() if k not in ("grid", "scan")})
    print("bias", out["coupling_bias"])


if __name__ == "__main__":
    main()
