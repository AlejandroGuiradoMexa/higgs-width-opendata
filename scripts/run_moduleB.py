"""Module B: Higgs width from the on-shell / off-shell ratio in ZZ -> 4l.

Writes results/moduleB.json and figures. Requires module A (processed parquet files) and the
MCFM outputs in generators/mcfm/output/.
"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mplhep as hep  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from hwidth.fit_onshell import histogram, smoothed_histogram  # noqa: E402
from hwidth.offshell import (Assumptions, OffshellTemplates, OnshellTemplates, WidthFit, asimov,  # noqa: E402
                             offshell_expectation, read_mcfm_hist, rebin_mcfm)
from hwidth.process import PROCESSED  # noqa: E402
from hwidth.samples import LUMI_FB  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
MCFM = ROOT / "generators" / "mcfm" / "output"
RES = ROOT / "results"
FIG = RES / "figures" / "mpl"  # quick-look matplotlib figures; publication figures: plots_root.py

# Off-shell bins lie on MCFM's 20 GeV bin boundaries (10 + 20k).
OFF_EDGES = np.array([230, 250, 290, 330, 370, 430, 510, 610, 810, 1010], dtype=float)
ON_EDGES = np.arange(110, 140 + 2.5, 2.5)
GAMMA_SM_MEV = 4.1
R_GRID = np.round(np.concatenate([np.arange(0, 5, 0.05), np.arange(5, 25.01, 0.25)]), 3)


def mcfm(name):
    path = next(MCFM.glob(f"*_{name}.dat"))
    return rebin_mcfm(*read_mcfm_hist(path), OFF_EDGES)


def load_inputs():
    dfs = {g: pd.read_parquet(PROCESSED / f"{g}.parquet") for g in ("data", "higgs", "zz", "other")}

    # Detector efficiency per bin: reconstructed qq->ZZ (2e2mu) / MCFM NLO fiducial qq->ZZ (2e2mu).
    zz_2e2mu = dfs["zz"][dfs["zz"]["channel"] == "2e2mu"]
    reco_2e2mu, _ = histogram(zz_2e2mu, OFF_EDGES)
    sigma_qq, _ = mcfm("qqZZ_NLO")
    eff = reco_2e2mu / (sigma_qq * LUMI_FB)

    sig_int, _ = mcfm("ggH_SI")
    interf, _ = mcfm("ggH_I")
    cont, _ = mcfm("ggZZ_B")
    signal = sig_int - interf

    # 2e2mu -> all channels (x2), fb -> events, and the reconstruction efficiency.
    fold = 2.0 * LUMI_FB * eff
    off = OffshellTemplates(
        edges=OFF_EDGES, S=signal * fold, I=interf * fold, B=cont * fold,
        Q=histogram(dfs["zz"], OFF_EDGES)[0], O=smoothed_histogram(dfs["other"], OFF_EDGES)[0], eff=eff)
    off_obs = np.histogram(dfs["data"]["m4l"], bins=OFF_EDGES)[0]

    on_bkg = histogram(dfs["zz"], ON_EDGES)[0] + smoothed_histogram(dfs["other"], ON_EDGES)[0]
    on = OnshellTemplates(S=histogram(dfs["higgs"], ON_EDGES)[0], B=on_bkg)
    on_obs = np.histogram(dfs["data"]["m4l"], bins=ON_EDGES)[0]
    return on, on_obs, off, off_obs


def interval(grid, q, level):
    imin = int(np.argmin(q))
    lo = grid[0] if q[:imin + 1].max() < level else float(np.interp(level, q[:imin + 1][::-1], grid[:imin + 1][::-1]))
    hi = WidthFit.upper_limit(grid, q, level)
    return [float(lo), float(hi)]


def summarize(grid, q):
    return {"best": float(grid[np.argmin(q)]), "int68": interval(grid, q, 1.0), "int95": interval(grid, q, 3.84),
            "scan": q.round(4).tolist()}


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    on, on_obs, off, off_obs = load_inputs()
    a = Assumptions()

    fit = WidthFit(on, on_obs, off, off_obs, a)
    best, _ = fit.best_fit()
    q_obs = fit.scan(R_GRID)

    # Expected: Asimov data at R = 1 with the qq->ZZ normalisation fitted to data.
    a_on, a_off = asimov(on, off, mu_on=1.0, R=1.0, mu_qq=best["mu_qq"], a=a)
    q_exp = WidthFit(on, a_on, off, a_off, a).scan(R_GRID)

    # Luminosity projections (expected, statistics + qq normalisation only).
    proj = {}
    for lumi in (36.6, 140.0, 300.0, 3000.0):
        r = lumi / LUMI_FB
        on_l, off_l = on.scaled(r), off.scaled(r)
        a_on, a_off = asimov(on_l, off_l, mu_on=1.0, R=1.0, mu_qq=best["mu_qq"], a=a)
        q = WidthFit(on_l, a_on, off_l, a_off, a).scan(R_GRID)
        proj[lumi] = WidthFit.upper_limit(R_GRID, q)

    # Yields table.
    k = a.K_gg
    yields = {"S": float(k * off.S.sum()), "I": float(k * off.I.sum()), "B_gg": float(k * off.B.sum()),
              "qqZZ": float(off.Q.sum()), "other": float(off.O.sum()), "data": int(off_obs.sum())}

    # Figure: off-shell m4l with SM and R = 10 expectations.
    hep.style.use("ATLAS")
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    centres = 0.5 * (OFF_EDGES[1:] + OFF_EDGES[:-1])
    widths = np.diff(OFF_EDGES)
    sm = offshell_expectation(off, best["mu_on"], best["mu_qq"], 0.0, a)
    bkg_only = offshell_expectation(off, 0.0, best["mu_qq"], 0.0, a)
    r10 = offshell_expectation(off, best["mu_on"] * 10, best["mu_qq"], 0.0, a)
    for y, lab, col, ls in ((sm, r"SM ($\Gamma_H = \Gamma_H^{SM}$)", "#00b4d8", "-"),
                            (r10, r"$\Gamma_H = 10\,\Gamma_H^{SM}$", "#ff6b6b", "--"),
                            (bkg_only, "no Higgs off shell", "#888888", ":")):
        ax.stairs(y / widths, OFF_EDGES, color=col, ls=ls, lw=2, label=lab)
    ax.errorbar(centres, off_obs / widths, yerr=np.sqrt(off_obs) / widths, fmt="ko", label="Data")
    ax.set_yscale("log")
    ax.set_ylabel("Events / GeV")
    ax.legend(fontsize=12)
    ax.text(0.04, 0.06, rf"ATLAS Open Data, {LUMI_FB} fb$^{{-1}}$ + MCFM gg$\to$ZZ", transform=ax.transAxes, fontsize=11)
    axr.plot(centres, r10 / sm, "o-", color="#ff6b6b", label=r"$10\,\Gamma^{SM}$ / SM")
    axr.axhline(1, color="k", lw=1)
    axr.set_ylabel("Ratio")
    axr.set_xlabel(r"$m_{4\ell}$ [GeV]")
    fig.savefig(FIG / "moduleB_offshell.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(R_GRID, q_obs, color="k", lw=2, label="Observed")
    ax.plot(R_GRID, q_exp, color="k", lw=2, ls="--", label="Expected (SM)")
    for lvl in (1.0, 3.84):
        ax.axhline(lvl, color="grey", lw=1, ls=":")
    ax.set_xlabel(r"$\Gamma_H / \Gamma_H^{SM}$")
    ax.set_ylabel(r"$-2\,\Delta\ln L$")
    ax.set_ylim(0, 8)
    ax.legend()
    fig.savefig(FIG / "moduleB_scan.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    out = {
        "lumi_fb": LUMI_FB, "gamma_sm_mev": GAMMA_SM_MEV, "K_gg": a.K_gg,
        "edges": OFF_EDGES.tolist(), "eff": off.eff.round(4).tolist(),
        "templates": {k_: getattr(off, k_).round(4).tolist() for k_ in ("S", "I", "B", "Q", "O")},
        "data": off_obs.tolist(), "yields": yields, "best_fit": best,
        "R_grid": R_GRID.tolist(), "observed": summarize(R_GRID, q_obs), "expected": summarize(R_GRID, q_exp),
        "projection_95CL_R": {str(k_): v for k_, v in proj.items()},
    }
    out["observed"]["width_95_mev"] = out["observed"]["int95"][1] * GAMMA_SM_MEV
    out["expected"]["width_95_mev"] = out["expected"]["int95"][1] * GAMMA_SM_MEV
    (RES / "moduleB.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k_: out[k_] for k_ in ("yields", "best_fit", "projection_95CL_R")}, indent=1))
    print("observed", {k_: v for k_, v in out["observed"].items() if k_ != "scan"})
    print("expected", {k_: v for k_, v in out["expected"].items() if k_ != "scan"})
    print("eff", out["eff"])


if __name__ == "__main__":
    main()
