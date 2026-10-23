"""Read the exactly4lep skims, apply the selection and store selected events as parquet."""
import json
import pathlib
import re

import awkward as ak
import numpy as np
import pandas as pd
import uproot

from .samples import LUMI_FB, SAMPLES
from .selection import select

ROOT = pathlib.Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"

LEP_BRANCHES = ["lep_pt", "lep_eta", "lep_phi", "lep_e", "lep_charge", "lep_type",
                "lep_isLooseID", "lep_isMediumID", "lep_isLooseIso", "lep_isTrigMatched"]
EVENT_BRANCHES = ["trigE", "trigM", "lep_n", "eventNumber", "runNumber"]
WEIGHT_BRANCHES = ["mcWeight", "xsec", "filteff", "kfac", "sum_of_weights", "ScaleFactor_PILEUP",
                   "ScaleFactor_ELE", "ScaleFactor_MUON", "ScaleFactor_LepTRIGGER"]


def files_for(group):
    """Local files belonging to a sample group."""
    if group == "data":
        return sorted(RAW.glob("*_data1[56]_*.root"))
    out = []
    for dsid in SAMPLES[group]:
        out += sorted(RAW.glob(f"*_mc_{dsid}.*.root"))
    return out


def mc_weight(arr):
    """Per-event weight normalised to LUMI_FB (signed mcWeight, so negative-weight generators stay unbiased)."""
    w = LUMI_FB * 1000.0 * arr["xsec"] * arr["filteff"] * arr["kfac"] * arr["mcWeight"] / arr["sum_of_weights"]
    for sf in WEIGHT_BRANCHES[5:]:
        w = w * arr[sf]
    return np.asarray(w, dtype=float)


def process_file(path, is_mc):
    branches = LEP_BRANCHES + EVENT_BRANCHES + (WEIGHT_BRANCHES if is_mc else [])
    tree = uproot.open(f"{path}:analysis")
    frames, cutflow = [], {}
    for arr in tree.iterate(branches, library="ak", step_size="200 MB"):
        arr = arr[arr["lep_n"] == 4]
        ev = {b: ak.to_numpy(arr[b]) for b in LEP_BRANCHES}
        ev["lep_isTrigMatched"] = ev["lep_isTrigMatched"].astype(bool)
        for b in ("lep_isLooseID", "lep_isMediumID", "lep_isLooseIso"):
            ev[b] = ev[b].astype(bool)
        ev["trigE"] = ak.to_numpy(arr["trigE"]).astype(bool)
        ev["trigM"] = ak.to_numpy(arr["trigM"]).astype(bool)
        res = select(ev)
        for k, v in res["cutflow"].items():
            cutflow[k] = cutflow.get(k, 0) + v
        m = res["mask"]
        df = pd.DataFrame({
            "m4l": res["m4l"][m], "m12": res["m12"][m], "m34": res["m34"][m],
            "channel": res["channel"][m],
            "weight": mc_weight(arr)[m] if is_mc else np.ones(m.sum()),
            "run": ak.to_numpy(arr["runNumber"])[m], "event": ak.to_numpy(arr["eventNumber"])[m],
        })
        for i in range(4):
            for b in ("lep_pt", "lep_eta", "lep_phi", "lep_e", "lep_charge", "lep_type"):
                df[f"{b}_{i}"] = ev[b][m, i]
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    dsid = re.search(r"_mc_(\d+)\.", path.name)
    df["dsid"] = int(dsid.group(1)) if dsid else 0
    return df, cutflow


def process_all():
    PROCESSED.mkdir(parents=True, exist_ok=True)
    cutflows = {}
    for group in SAMPLES:
        frames, total = [], {}
        for path in files_for(group):
            df, cf = process_file(path, is_mc=group != "data")
            frames.append(df)
            for k, v in cf.items():
                total[k] = total.get(k, 0) + v
        if not frames:
            print(f"{group}: no files found")
            continue
        df = pd.concat(frames, ignore_index=True)
        df.to_parquet(PROCESSED / f"{group}.parquet")
        cutflows[group] = total
        print(f"{group}: {len(df)} events, sum of weights {df['weight'].sum():.2f}")
    (PROCESSED / "cutflow.json").write_text(json.dumps(cutflows, indent=1))
    return cutflows
