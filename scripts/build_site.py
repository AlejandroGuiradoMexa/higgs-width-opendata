"""Bundle results/*.json into site/data.json (only what the web page plots) and copy figures."""
import json
import pathlib
import shutil

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
SITE = ROOT / "site"


def load(name):
    path = RES / f"{name}.json"
    return json.loads(path.read_text()) if path.exists() else None


def main():
    a, b, c, d, e = (load(n) for n in ("moduleA", "moduleB", "moduleC", "moduleD", "moduleE"))
    data = {}
    if a:
        data["A"] = {"fit": {k: a["fit"][k] for k in ("mu_hat", "mu_lo", "mu_hi", "Z_obs", "Z_exp", "scan_mu", "scan_2dnll")},
                     "yields": a["yields_115_130"], "channels": a["channels_data_115_130"],
                     "hist": a["hist_m4l"], "lumi": a["lumi_fb"], "cutflow": a["cutflow"].get("data")}
    if b:
        data["B"] = {k: b[k] for k in ("edges", "templates", "data", "yields", "best_fit", "R_grid", "K_gg",
                                       "projection_95CL_R", "gamma_sm_mev", "eff")}
        data["B"]["observed"] = b["observed"]
        data["B"]["expected"] = b["expected"]
    if c:
        data["C"] = {"scenarios": [{k: s[k] for k in ("key", "label", "group", "obs95", "exp95", "hl_best", "hl68", "hl95")}
                                   for s in c["scenarios"]],
                     "decoupled": {k: c["decoupled"][k] for k in ("mu_off_best", "mu_off_68", "mu_off_95")},
                     "coupling_bias": c["coupling_bias"], "indirect": c["indirect"]}
    if d:
        data["D"] = d
    if e:
        data["E"] = e
    SITE.mkdir(exist_ok=True)
    (SITE / "data.json").write_text(json.dumps(data, separators=(",", ":")))
    figs = SITE / "figures"
    figs.mkdir(exist_ok=True)
    for f in (RES / "figures").glob("*.png"):
        shutil.copy(f, figs / f.name)
    print(f"site/data.json: {(SITE / 'data.json').stat().st_size / 1024:.0f} kB, modules: {sorted(data)}")


if __name__ == "__main__":
    main()
