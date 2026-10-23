"""Module A: H->ZZ*->4l on-shell analysis. Writes results/moduleA.json and figures."""
import json
import sys
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mplhep as hep  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from hwidth.fit_onshell import build_workspace, fit, histogram, smoothed_histogram  # noqa: E402
from hwidth.process import PROCESSED, process_all  # noqa: E402
from hwidth.samples import LUMI_FB  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures" / "mpl"  # quick-look matplotlib figures; publication figures: plots_root.py
GROUPS = ["other", "zz", "higgs"]
LABELS = {"other": r"$Z$+jets, $t\bar{t}$, $VVV$", "zz": r"$ZZ^{*}$", "higgs": r"$H\to ZZ^{*}$ ($m_H$=125 GeV)"}
COLORS = {"other": "#6b59d3", "zz": "#ff6b6b", "higgs": "#00b4d8"}


def template(dfs, group, edges):
    """Reducible background ('other') is smoothed; the other samples are plain histograms."""
    return smoothed_histogram(dfs[group], edges) if group == "other" else histogram(dfs[group], edges)


def load():
    return {g: pd.read_parquet(PROCESSED / f"{g}.parquet") for g in ["data"] + GROUPS}


def plot_stack(dfs, edges, fname, logy=False, xlabel=r"$m_{4\ell}$ [GeV]"):
    hep.style.use("ATLAS")
    fig, ax = plt.subplots(figsize=(8, 6))
    hists = [template(dfs, g, edges)[0] for g in GROUPS]
    hep.histplot(hists, edges, stack=True, histtype="fill", color=[COLORS[g] for g in GROUPS],
                 label=[LABELS[g] for g in GROUPS], ax=ax)
    d, _ = np.histogram(dfs["data"]["m4l"], bins=edges)
    centres = 0.5 * (edges[1:] + edges[:-1])
    ax.errorbar(centres[d > 0], d[d > 0], yerr=np.sqrt(d[d > 0]), fmt="ko", label="Data")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"Events / {edges[1] - edges[0]:g} GeV")
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(0.05, None)
    else:
        ax.set_ylim(0, None)
    ax.text(0.04, 0.93, "ATLAS Open Data", transform=ax.transAxes, fontsize=15, weight="bold")
    ax.text(0.04, 0.86, rf"$\sqrt{{s}}$ = 13 TeV, {LUMI_FB} fb$^{{-1}}$, education release",
            transform=ax.transAxes, fontsize=12)
    ax.legend(loc="upper right", fontsize=12)
    fig.savefig(FIG / fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return {"edges": edges.tolist(), "data": d.tolist(), **{g: h.tolist() for g, h in zip(GROUPS, hists)}}


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    if "--reprocess" in sys.argv or not (PROCESSED / "cutflow.json").exists():
        process_all()
    cutflows = json.loads((PROCESSED / "cutflow.json").read_text())
    dfs = load()

    full = plot_stack(dfs, np.arange(80, 250 + 2.5, 2.5), "moduleA_m4l.png")
    offshell = plot_stack(dfs, np.arange(220, 1020, 40.0), "moduleA_m4l_offshell.png", logy=True)

    edges = np.arange(110, 140 + 2.5, 2.5)
    templates = {g: template(dfs, g, edges) for g in GROUPS}
    data_counts, _ = np.histogram(dfs["data"]["m4l"], bins=edges)
    ws = build_workspace(templates, data_counts)
    res = fit(ws)
    (RES / "moduleA_workspace.json").write_text(json.dumps(dict(ws), indent=1))

    window = {g: float(dfs[g].query("115 < m4l < 130")["weight"].sum()) for g in GROUPS}
    window["data"] = int(((dfs["data"]["m4l"] > 115) & (dfs["data"]["m4l"] < 130)).sum())

    out = {
        "lumi_fb": LUMI_FB,
        "fit_range_gev": [float(edges[0]), float(edges[-1])],
        "fit": res,
        "yields_115_130": window,
        "channels_data_115_130": dfs["data"].query("115 < m4l < 130")["channel"].value_counts().to_dict(),
        "cutflow": cutflows,
        "hist_m4l": full,
        "hist_offshell": offshell,
    }
    (RES / "moduleA.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("yields_115_130", "channels_data_115_130")}, indent=1))
    print({k: res[k] for k in ("mu_hat", "mu_lo", "mu_hi", "Z_obs", "Z_exp")})


if __name__ == "__main__":
    main()
