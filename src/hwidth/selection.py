"""Event selection for ZZ(*) -> 4l, following the structure of the ATLAS H->4l analyses.

All functions work on numpy arrays of shape (n_events, 4) for the four leptons of the
exactly4lep skim. Energies and momenta are in GeV.
"""
import numpy as np

MZ = 91.1876

# The three ways of splitting four leptons into two pairs.
PAIRINGS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


def four_momenta(pt, eta, phi, e):
    """Return (px, py, pz, E) arrays, each of shape (n, 4)."""
    return pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta), e


def inv_mass(p4, idx):
    """Invariant mass of the sum of leptons `idx` for every event."""
    px, py, pz, e = (comp[:, list(idx)].sum(axis=1) for comp in p4)
    return np.sqrt(np.maximum(e**2 - px**2 - py**2 - pz**2, 0.0))


def is_sfos(lep_type, charge, i, j):
    """Same flavour, opposite sign."""
    return (lep_type[:, i] == lep_type[:, j]) & (charge[:, i] + charge[:, j] == 0)


def m34_threshold(m4l):
    """Lower cut on the subleading pair mass: 12 GeV below 140 GeV, rising linearly to 50 GeV at 190 GeV."""
    return np.clip(12.0 + (m4l - 140.0) * (50.0 - 12.0) / (190.0 - 140.0), 12.0, 50.0)


def pair_leptons(lep_type, charge, p4):
    """Pick the pairing whose leading pair is closest to the Z mass.

    Returns (m12, m34, valid), where the leading pair Z1 is the SFOS pair closest to mZ
    and Z2 is the remaining pair, which must also be SFOS.
    """
    n = len(lep_type)
    best_dist = np.full(n, np.inf)
    m12 = np.zeros(n)
    m34 = np.zeros(n)
    for a, b in PAIRINGS:
        ok = is_sfos(lep_type, charge, *a) & is_sfos(lep_type, charge, *b)
        ma, mb = inv_mass(p4, a), inv_mass(p4, b)
        for z1, z2 in ((ma, mb), (mb, ma)):
            dist = np.where(ok, np.abs(z1 - MZ), np.inf)
            better = dist < best_dist
            best_dist = np.where(better, dist, best_dist)
            m12 = np.where(better, z1, m12)
            m34 = np.where(better, z2, m34)
    return m12, m34, np.isfinite(best_dist)


def delta_r_min(eta, phi):
    """Smallest Delta R between any two of the four leptons."""
    dr = np.full(len(eta), np.inf)
    for i in range(4):
        for j in range(i + 1, 4):
            dphi = np.mod(phi[:, i] - phi[:, j] + np.pi, 2 * np.pi) - np.pi
            dr = np.minimum(dr, np.hypot(eta[:, i] - eta[:, j], dphi))
    return dr


def jpsi_veto(lep_type, charge, p4, mmin=5.0):
    """True if no SFOS pair has a mass below `mmin` GeV."""
    keep = np.ones(len(lep_type), dtype=bool)
    for i in range(4):
        for j in range(i + 1, 4):
            low = is_sfos(lep_type, charge, i, j) & (inv_mass(p4, (i, j)) < mmin)
            keep &= ~low
    return keep


def select(ev):
    """Apply the 4l selection to a dict of numpy arrays and return a dict of selected quantities.

    `ev` holds per-lepton (n, 4) arrays: lep_pt, lep_eta, lep_phi, lep_e, lep_charge, lep_type,
    lep_isLooseID, lep_isMediumID, lep_isLooseIso, lep_isTrigMatched; and per-event arrays
    trigE, trigM.
    """
    t, q = ev["lep_type"], ev["lep_charge"]
    p4 = four_momenta(ev["lep_pt"], ev["lep_eta"], ev["lep_phi"], ev["lep_e"])

    pt_sorted = -np.sort(-ev["lep_pt"], axis=1)
    cut = {}
    cut["trigger"] = (ev["trigE"] | ev["trigM"]) & ev["lep_isTrigMatched"].any(axis=1)
    ident = np.where(t == 11, ev["lep_isLooseID"], ev["lep_isMediumID"]) & ev["lep_isLooseIso"]
    cut["id_iso"] = ident.all(axis=1)
    cut["pt"] = (pt_sorted[:, 0] > 20) & (pt_sorted[:, 1] > 15) & (pt_sorted[:, 2] > 10)
    cut["charge_flavour"] = (q.sum(axis=1) == 0) & np.isin(t.sum(axis=1), (44, 48, 52))

    m12, m34, paired = pair_leptons(t, q, p4)
    m4l = inv_mass(p4, (0, 1, 2, 3))
    cut["pairing"] = paired
    cut["m12"] = (m12 > 50) & (m12 < 106)
    cut["m34"] = (m34 > m34_threshold(m4l)) & (m34 < 115)
    cut["dr"] = delta_r_min(ev["lep_eta"], ev["lep_phi"]) > 0.1
    cut["jpsi"] = jpsi_veto(t, q, p4)

    passed = np.ones(len(t), dtype=bool)
    cutflow = {}
    for name, c in cut.items():
        passed &= c
        cutflow[name] = int(passed.sum())

    channel = np.select([t.sum(axis=1) == 44, t.sum(axis=1) == 52], ["4e", "4mu"], "2e2mu")
    return {
        "mask": passed,
        "m4l": m4l,
        "m12": m12,
        "m34": m34,
        "channel": channel,
        "cutflow": cutflow,
    }
