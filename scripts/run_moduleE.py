"""Module E: electron shower shapes in the g4calo ATLAS-like LAr/Pb calorimeter.

Reads g4calo/out/*.txt and writes results/moduleE.json and figures.
"""
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mplhep as hep  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

from hwidth.showers import read, roc, shower_shapes  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "g4calo" / "out"
RES = ROOT / "results"
FIG = RES / "figures" / "mpl"  # quick-look matplotlib figures; publication figures: plots_root.py

PARTICLES = ["e-", "gamma", "pi0", "pi-"]
LABEL = {"e-": r"$e^-$", "gamma": r"$\gamma$", "pi0": r"$\pi^0$", "pi-": r"$\pi^-$"}
COLOR = {"e-": "#00b4d8", "gamma": "#ffb703", "pi0": "#ff6b6b", "pi-": "#6b59d3"}
ENERGIES = [10, 25, 50, 100]
VARS = {
    "R_eta": (r"$R_\eta$", (0.80, 1.0)),
    "R_phi": (r"$R_\phi$", (0.85, 1.0)),
    "w_eta2": (r"$w_{\eta 2}$", (0.006, 0.016)),
    "ws_tot": (r"$w_{s,\mathrm{tot}}$ [strips]", (0, 8)),
    "E_ratio": (r"$E_\mathrm{ratio}$", (0, 1)),
    "f1": (r"$f_1$", (0, 0.6)),
    "f3": (r"$f_3$", (0, 0.05)),
    "R_had": (r"$R_\mathrm{had}$ (visible)", (0, 0.5)),
    "E_over_p": (r"$E/p$", (0, 1.3)),
}
# Direction in which each variable is more electron-like.
ELECTRON_HIGH = {"R_eta": True, "R_phi": True, "w_eta2": False, "ws_tot": False, "E_ratio": True,
                 "f1": True, "f3": False, "R_had": False, "E_over_p": True}
SAMPLING = {}  # filled from 50 GeV electrons: visible LAr energy / beam energy


def em_visible(ev):
    """Visible energy in the LAr layers (presampler + accordion), excluding the hadronic scintillator."""
    return ev["active"] - ev["had"].sum(axis=(1, 2))


def load(p, E, x0=2.0):
    path = OUT / f"{p}_{E}GeV_X0{x0:.1f}.txt"
    if not path.exists():
        return None
    ev = read(path)
    v = shower_shapes(ev)
    if SAMPLING:
        # E/p: calibrated EM energy over the true momentum (what the tracker would measure for e/pi+-).
        v["E_over_p"] = em_visible(ev) / SAMPLING["f"] / (1000.0 * E)
    return ev, v


# Single events for the event display on the project page: the first N_EVENTS of each sample, unselected.
EVENT_SAMPLES = [("e-", 50, 2.0), ("gamma", 50, 2.0), ("pi0", 50, 2.0), ("pi-", 50, 2.0), ("e-", 50, 0.0), ("e-", 50, 6.0)]
N_EVENTS = 5


def example_events():
    out = []
    for p, E, x0 in EVENT_SAMPLES:
        path = OUT / f"{p}_{E}GeV_X0{x0:.1f}.txt"
        if not path.exists():
            continue
        ev = read(path)
        layers = ("ps", "s1", "s2", "s3", "had")
        out.append({"particle": p, "E": E, "x0": x0, "events": [
            {**{l: np.round(ev[l][i], 2).tolist() for l in layers}, "depth": np.round(ev["depth"][i], 1).tolist()}
            for i in range(N_EVENTS)]})
    (RES / "moduleE_events.json").write_text(json.dumps(out, separators=(",", ":")))
    return out


def hist_json(x, lo, hi, n=40):
    h, e = np.histogram(np.clip(x, lo, hi), bins=n, range=(lo, hi), density=True)
    return {"edges": e.round(6).tolist(), "density": h.round(6).tolist()}


def gamma_profile(t, a, b, norm):
    """Longitudinal profile dE/dt = norm * b (bt)^(a-1) e^(-bt) / Gamma(a) (Longo-Sestili)."""
    from scipy.special import gamma as G
    return norm * b * (b * t) ** (a - 1) * np.exp(-b * t) / G(a)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    hep.style.use("ATLAS")
    out = {"distributions_50GeV": {}, "means": {}, "resolution": {}, "profile": {}, "material": {}, "roc": {}}

    ref = load("e-", 50)
    SAMPLING["f"] = float(np.mean(em_visible(ref[0])) / 50000.0)

    # 1. Shower-shape distributions at 50 GeV.
    shapes50 = {p: load(p, 50) for p in PARTICLES}
    shapes25 = {p: load(p, 25) for p in PARTICLES}
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, var in zip(axes.flat, ["R_eta", "w_eta2", "ws_tot", "E_ratio", "R_had", "E_over_p"]):
        lab, (lo, hi) = VARS[var]
        for p in PARTICLES:
            if shapes50[p] is None:
                continue
            x = shapes50[p][1][var]
            ax.hist(np.clip(x, lo, hi), bins=40, range=(lo, hi), density=True, histtype="step", lw=2,
                    color=COLOR[p], label=LABEL[p])
            out["distributions_50GeV"].setdefault(var, {})[p] = hist_json(x, lo, hi)
        ax.set_xlabel(lab)
        ax.set_yticks([])
    axes.flat[0].legend()
    fig.suptitle("Shower shapes at 50 GeV, 2 $X_0$ upstream (g4calo, Geant4 11.4, FTFP_BERT)", fontsize=14)
    fig.savefig(FIG / "moduleE_shapes_50GeV.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 2. Means vs energy, and energy resolution of electrons.
    sampling = SAMPLING["f"]
    for p in PARTICLES:
        for E in ENERGIES:
            r = load(p, E)
            if r is None:
                continue
            ev, v = r
            out["means"].setdefault(p, {})[E] = {k: [float(np.mean(v[k])), float(np.std(v[k]))] for k in VARS}
    for p in ("e-", "gamma"):
        es, res, err = [], [], []
        for E in ENERGIES:
            r = load(p, E)
            if r is None:
                continue
            ev, _ = r
            reco = em_visible(ev) / sampling / 1000.0
            core = reco[np.abs(reco - np.median(reco)) < 3 * np.std(reco)]
            es.append(E)
            res.append(float(np.std(core) / np.mean(core)))
            err.append(float(np.std(core) / np.mean(core) / np.sqrt(2 * len(core))))
        if len(es) >= 3:
            f = lambda E, a, b: np.sqrt((a / np.sqrt(E)) ** 2 + b**2)
            (a, b), cov = curve_fit(f, np.array(es), np.array(res), p0=[0.1, 0.01], sigma=np.array(err))
            out["resolution"][p] = {"E": es, "sigma_over_E": res, "err": err, "stochastic": float(abs(a)),
                                    "constant": float(abs(b))}

    # 3. Longitudinal profiles at 50 GeV with a Gamma-function fit.
    t = np.arange(30) + 0.5
    fig, ax = plt.subplots(figsize=(8, 6))
    for p in PARTICLES:
        if shapes50[p] is None:
            continue
        prof = shapes50[p][0]["depth"].mean(axis=0)
        prof = prof / prof.sum()
        ax.step(t, prof, where="mid", color=COLOR[p], label=LABEL[p], lw=2)
        entry = {"t": t.tolist(), "dEdt": prof.round(6).tolist()}
        if p in ("e-", "gamma"):
            (a, b, n), _ = curve_fit(gamma_profile, t, prof, p0=[4, 0.5, 1])
            entry["t_max"] = float((a - 1) / b)
            ax.plot(t, gamma_profile(t, a, b, n), color=COLOR[p], ls="--", lw=1)
        out["profile"][p] = entry
    ax.set_xlabel(r"Depth in calorimeter [$X_0$]")
    ax.set_ylabel("Fraction of active energy / $X_0$")
    ax.legend()
    fig.savefig(FIG / "moduleE_profile.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 4. Dependence on upstream material (the classic source of data/MC shower-shape shifts).
    for p in ("e-", "gamma"):
        rows = {}
        for x0 in (0.0, 1.0, 2.0, 4.0, 6.0):
            r = load(p, 50, x0)
            if r is None:
                continue
            ev, v = r
            reco = em_visible(ev) / sampling / 1000.0
            rows[x0] = {"R_eta": float(np.mean(v["R_eta"])), "w_eta2": float(np.mean(v["w_eta2"])),
                        "f1": float(np.mean(v["f1"])), "f_ps": float(np.mean(v["f_ps"])),
                        "E_mean_GeV": float(np.mean(reco)), "E_res": float(np.std(reco) / np.mean(reco))}
        out["material"][p] = rows

    # 5. Single-variable ROC curves at 25 GeV, where the two pi0 photons are ~11 mm (2.3 strips)
    # apart: photons vs pi0 with the strip variables (photon ID), electrons vs pi- with hadronic
    # leakage and E/p (electron ID).
    if shapes25["gamma"] is not None:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for ax, sig, bkg, variables in ((axes[0], "gamma", "pi0", ["E_ratio", "ws_tot", "w_eta2", "R_eta"]),
                                        (axes[1], "e-", "pi-", ["R_had", "E_over_p", "R_eta", "w_eta2"])):
            if shapes25[bkg] is None or shapes25[sig] is None:
                continue
            for var in variables:
                es, eb = roc(shapes25[sig][1][var], shapes25[bkg][1][var], ELECTRON_HIGH[var])
                rej = 1.0 / np.maximum(eb, 1.0 / len(eb) / 10)
                ax.plot(es, rej, lw=2, label=VARS[var][0])
                out["roc"].setdefault(f"{sig}_vs_{bkg}", {})[var] = {"eff_sig": es.round(4).tolist(),
                                                                     "eff_bkg": eb.round(5).tolist()}
            ax.set_yscale("log")
            ax.set_xlim(0.5, 1.0)
            ax.set_xlabel(f"{LABEL[sig]} efficiency")
            ax.set_ylabel(f"{LABEL[bkg]} rejection")
            ax.set_title("25 GeV", fontsize=13)
            ax.legend()
        fig.savefig(FIG / "moduleE_roc.png", dpi=130, bbox_inches="tight")
        plt.close(fig)

    out["sampling_fraction"] = sampling
    (RES / "moduleE.json").write_text(json.dumps(out, indent=1))
    example_events()
    print(json.dumps({"sampling": sampling, "resolution": {k: (v["stochastic"], v["constant"]) for k, v in out["resolution"].items()},
                      "material_e": out["material"].get("e-")}, indent=1))


if __name__ == "__main__":
    main()
