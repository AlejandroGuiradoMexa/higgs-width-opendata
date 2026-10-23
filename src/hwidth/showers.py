"""ATLAS electron/photon shower-shape variables computed from the g4calo readout windows.

Cell windows (eta x phi): presampler 5x3, strips 41x3, middle 11x11, back 5x11, hadronic 3x3,
all centred on the impact point. Definitions follow the ATLAS e/gamma identification papers
(e.g. Eur. Phys. J. C 79 (2019) 639).
"""
import numpy as np

SHAPES = {"ps": (5, 3), "s1": (41, 3), "s2": (11, 11), "s3": (5, 11), "had": (3, 3)}
N_DEPTH = 30
DETA_S1 = 0.003125
DETA_S2 = 0.025


def read(path):
    """Read a g4calo output file into a dict of arrays (energies in MeV)."""
    ncol = sum(nx * ny for nx, ny in SHAPES.values()) + N_DEPTH + 2
    with open(path) as f:
        # Skip the header and any line cut short while the simulation is still writing.
        rows = [line.split() for line in f if not line.startswith("#")]
    raw = np.array([r for r in rows if len(r) == ncol], dtype=float).reshape(-1, ncol)
    out, col = {}, 0
    for name, (nx, ny) in SHAPES.items():
        out[name] = raw[:, col:col + nx * ny].reshape(-1, nx, ny)
        col += nx * ny
    out["depth"] = raw[:, col:col + N_DEPTH]
    col += N_DEPTH
    out["active"], out["passive"] = raw[:, col], raw[:, col + 1]
    return out


def window(cells, neta, nphi):
    """Sum of a neta x nphi window centred on the middle cell, for every event."""
    _, nx, ny = cells.shape
    cx, cy = nx // 2, ny // 2
    return cells[:, cx - neta // 2:cx + neta // 2 + 1, cy - nphi // 2:cy + nphi // 2 + 1].sum(axis=(1, 2))


def second_maximum(strips):
    """Energy of the second-highest local maximum in each strip profile (0 if none)."""
    left = np.pad(strips, ((0, 0), (1, 0)))[:, :-1]
    right = np.pad(strips, ((0, 0), (0, 1)))[:, 1:]
    is_max = (strips > left) & (strips >= right) & (strips > 0)
    vals = np.sort(np.where(is_max, strips, 0.0), axis=1)
    return vals[:, -2]


def shower_shapes(ev):
    """Return a dict of per-event shower-shape variables."""
    s2 = ev["s2"]
    e_s2_37 = window(s2, 3, 7)
    e_s2_77 = window(s2, 7, 7)
    e_s2_33 = window(s2, 3, 3)

    # Lateral width in eta in a 3x5 window of the middle layer.
    sub = s2[:, 4:7, 3:8].sum(axis=2)
    eta = (np.arange(3) - 1) * DETA_S2
    esum = sub.sum(axis=1)
    mean = (sub * eta).sum(axis=1) / esum
    w_eta2 = np.sqrt(np.maximum((sub * eta**2).sum(axis=1) / esum - mean**2, 0.0))

    # Strip profile (sum over phi), total width over 20 strips around the maximum.
    strips = ev["s1"].sum(axis=2)
    imax = strips.argmax(axis=1)
    idx = np.arange(strips.shape[1])
    dist = idx[None, :] - imax[:, None]
    in20 = np.abs(dist) <= 10
    e20 = np.where(in20, strips, 0.0)
    ws_tot = np.sqrt((e20 * dist**2).sum(axis=1) / np.maximum(e20.sum(axis=1), 1e-9))
    e1 = strips.max(axis=1)
    e2 = second_maximum(strips)
    e_ratio = (e1 - e2) / np.maximum(e1 + e2, 1e-9)

    e_s1 = strips.sum(axis=1)
    e_s3 = ev["s3"].sum(axis=(1, 2))
    e_ps = ev["ps"].sum(axis=(1, 2))
    e_acc = e_s1 + s2.sum(axis=(1, 2)) + e_s3
    return {
        "R_eta": e_s2_37 / np.maximum(e_s2_77, 1e-9),
        "R_phi": e_s2_33 / np.maximum(e_s2_37, 1e-9),
        "w_eta2": w_eta2,
        "ws_tot": ws_tot,
        "E_ratio": e_ratio,
        "f1": e_s1 / np.maximum(e_acc, 1e-9),
        "f3": e_s3 / np.maximum(e_acc, 1e-9),
        "f_ps": e_ps / np.maximum(e_acc + e_ps, 1e-9),
        # Hadronic leakage: visible scintillator energy over visible LAr energy (uncalibrated ratio).
        "R_had": ev["had"].sum(axis=(1, 2)) / np.maximum(e_acc + e_ps, 1e-9),
        "E_active": ev["active"],
    }


def roc(signal, background, larger_is_signal=True, n=200):
    """Signal efficiency vs background rejection for a single variable."""
    s, b = (signal, background) if larger_is_signal else (-signal, -background)
    cuts = np.quantile(np.concatenate([s, b]), np.linspace(0, 1, n))
    eff_s = np.array([(s >= c).mean() for c in cuts])
    eff_b = np.array([(b >= c).mean() for c in cuts])
    return eff_s, eff_b
