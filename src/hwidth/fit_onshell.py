"""On-shell H->ZZ*->4l signal-strength fit with pyhf."""
import numpy as np
import pyhf
from scipy.special import erf as _erf

pyhf.set_backend("numpy", "minuit")

# Normalisation uncertainties (fractional), typical of the published H->4l analyses.
NORM_UNC = {"zz": 0.10, "other": 0.50}
LUMI_UNC = 0.017


def histogram(df, edges):
    """Weighted histogram and sqrt(sum w^2)."""
    h, _ = np.histogram(df["m4l"], bins=edges, weights=df["weight"])
    w2, _ = np.histogram(df["m4l"], bins=edges, weights=df["weight"] ** 2)
    return h, np.sqrt(w2)


def smoothed_histogram(df, edges, bandwidth=10.0):
    """Histogram from a per-sample Gaussian kernel estimate.

    Each sample's m4l shape is taken from its unweighted events and scaled to the sample's
    total weight. This removes spikes from the few high-weight events of the Z+jets and ttbar
    samples. The statistical uncertainty is set to the relative MC statistics of the sample.
    """
    h = np.zeros(len(edges) - 1)
    var = np.zeros(len(edges) - 1)
    for _, s in df.groupby("dsid"):
        x = s["m4l"].to_numpy()
        cdf = lambda e: 0.5 * (1 + _erf((e[None, :] - x[:, None]) / (np.sqrt(2) * bandwidth)))
        frac = np.diff(cdf(np.asarray(edges, dtype=float)), axis=1).mean(axis=0)
        total = s["weight"].sum()
        h += total * frac
        var += (total * frac) ** 2 / max(len(s), 1)
    return h, np.sqrt(var)



def build_workspace(templates, data_counts):
    """templates: {name: (counts, stat_unc)} with keys 'higgs', 'zz', 'other'."""
    samples = []
    for name, (counts, unc) in templates.items():
        mods = [
            {"name": "lumi", "type": "normsys", "data": {"hi": 1 + LUMI_UNC, "lo": 1 - LUMI_UNC}},
            {"name": "staterror_sr", "type": "staterror", "data": list(np.maximum(unc, 1e-6))},
        ]
        if name == "higgs":
            mods.append({"name": "mu", "type": "normfactor", "data": None})
        else:
            u = NORM_UNC[name]
            mods.append({"name": f"norm_{name}", "type": "normsys", "data": {"hi": 1 + u, "lo": 1 - u}})
        samples.append({"name": name, "data": list(np.maximum(counts, 1e-6)), "modifiers": mods})
    spec = {
        "channels": [{"name": "sr", "samples": samples}],
        "observations": [{"name": "sr", "data": list(map(float, data_counts))}],
        "measurements": [{"name": "meas", "config": {"poi": "mu", "parameters": []}}],
        "version": "1.0.0",
    }
    return pyhf.Workspace(spec)


def crossings(grid, scan, level, best):
    """Linear interpolation of where the scan crosses `level` on each side of `best`."""
    left, right = grid <= best, grid >= best
    lo = np.interp(level, scan[left][::-1], grid[left][::-1]) if scan[left].max() > level else grid[0]
    hi = np.interp(level, scan[right], grid[right]) if scan[right].max() > level else grid[-1]
    return float(lo), float(hi)


def fit(ws):
    """Return best-fit mu with profile-likelihood interval, plus significances."""
    model = ws.model()
    data = ws.data(model)
    poi = model.config.poi_index

    bestfit, twice_nll = pyhf.infer.mle.fit(data, model, return_fitted_val=True)
    mu_hat = float(bestfit[poi])

    # Profile-likelihood scan of -2 Delta lnL in mu.
    grid = np.linspace(0.0, max(3.0, 3 * mu_hat), 301)
    scan = []
    for mu in grid:
        _, t = pyhf.infer.mle.fixed_poi_fit(mu, data, model, return_fitted_val=True)
        scan.append(float(t - twice_nll))
    scan = np.array(scan)
    lo, hi = crossings(grid, scan, 1.0, mu_hat)

    # Discovery significance (q0) observed and expected (Asimov with mu = 1).
    init, bounds = model.config.suggested_init(), model.config.suggested_bounds()
    bounds[poi] = (0.0, 10.0)
    q0_obs = pyhf.infer.test_statistics.q0(0.0, data, model, init, bounds, model.config.suggested_fixed())
    asimov = pyhf.infer.calculators.generate_asimov_data(
        1.0, data, model, init, bounds, model.config.suggested_fixed())
    q0_exp = pyhf.infer.test_statistics.q0(0.0, asimov, model, init, bounds, model.config.suggested_fixed())

    return {
        "mu_hat": mu_hat,
        "mu_lo": lo,
        "mu_hi": hi,
        "Z_obs": float(np.sqrt(max(q0_obs, 0))),
        "Z_exp": float(np.sqrt(max(q0_exp, 0))),
        "scan_mu": grid.tolist(),
        "scan_2dnll": scan.tolist(),
        "pulls": {name: bestfit[model.config.par_slice(name)].tolist() for name in model.config.par_order},
    }
