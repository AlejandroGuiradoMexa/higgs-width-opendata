import numpy as np

from hwidth.selection import MZ, four_momenta, inv_mass, m34_threshold, pair_leptons


def back_to_back(m, pt_scale=40.0):
    """Two massless leptons along x with invariant mass m (in their rest frame, boosted along nothing)."""
    e = m / 2.0
    return [e, e], [0.0, 0.0], [0.0, np.pi], [e, e]  # pt, eta, phi, E


def build(pairs):
    """Build one event from two (mass, eta-offset) pairs; returns p4 and arrays."""
    pt, eta, phi, e = [], [], [], []
    for m, eta0 in pairs:
        p, _, ph, en = back_to_back(m)
        pt += p
        eta += [eta0, eta0]
        phi += ph
        e += [x * np.cosh(eta0) for x in en]
    arr = lambda x: np.array([x], dtype=float)
    return four_momenta(arr(pt), arr(eta), arr(phi), arr(e))


def test_inv_mass_of_pair():
    p4 = build([(MZ, 0.0), (30.0, 1.0)])
    assert np.isclose(inv_mass(p4, (0, 1))[0], MZ)
    assert np.isclose(inv_mass(p4, (2, 3))[0], 30.0)


def test_pairing_2e2mu_picks_flavour_pairs():
    p4 = build([(MZ, 0.0), (30.0, 1.0)])
    t = np.array([[11, 11, 13, 13]])
    q = np.array([[1, -1, 1, -1]])
    m12, m34, ok = pair_leptons(t, q, p4)
    assert ok[0] and np.isclose(m12[0], MZ) and np.isclose(m34[0], 30.0)


def test_pairing_4mu_puts_z_like_pair_first():
    p4 = build([(30.0, 1.0), (MZ, 0.0)])
    t = np.array([[13, 13, 13, 13]])
    q = np.array([[1, -1, 1, -1]])
    m12, m34, ok = pair_leptons(t, q, p4)
    assert ok[0] and np.isclose(m12[0], MZ, atol=1e-6)


def test_same_sign_event_is_not_paired():
    p4 = build([(MZ, 0.0), (30.0, 1.0)])
    t = np.array([[11, 11, 13, 13]])
    q = np.array([[1, 1, -1, -1]])
    assert not pair_leptons(t, q, p4)[2][0]


def test_m34_threshold_endpoints():
    assert np.allclose(m34_threshold(np.array([100.0, 140.0, 165.0, 190.0, 500.0])), [12, 12, 31, 50, 50])
