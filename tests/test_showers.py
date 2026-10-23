import numpy as np

from hwidth.showers import SHAPES, second_maximum, shower_shapes, window


def empty_event():
    ev = {name: np.zeros((1, nx, ny)) for name, (nx, ny) in SHAPES.items()}
    ev["active"] = np.ones(1)
    return ev


def test_point_like_shower_is_fully_contained():
    ev = empty_event()
    ev["s2"][0, 5, 5] = 1000.0
    ev["s1"][0, 20, 1] = 300.0
    v = shower_shapes(ev)
    assert np.isclose(v["R_eta"][0], 1.0)
    assert np.isclose(v["R_phi"][0], 1.0)
    assert np.isclose(v["w_eta2"][0], 0.0)
    assert np.isclose(v["E_ratio"][0], 1.0)
    assert np.isclose(v["f1"][0], 300.0 / 1300.0)


def test_energy_outside_core_lowers_r_eta():
    ev = empty_event()
    ev["s2"][0, 5, 5] = 800.0
    ev["s2"][0, 2, 5] = 200.0  # inside 7x7, outside 3x7
    assert np.isclose(shower_shapes(ev)["R_eta"][0], 0.8)


def test_two_photon_strip_profile_has_low_e_ratio():
    strips = np.zeros((1, 41))
    strips[0, 18] = 100.0
    strips[0, 24] = 80.0
    assert np.isclose(second_maximum(strips)[0], 80.0)


def test_window_sums_centre():
    cells = np.arange(121, dtype=float).reshape(1, 11, 11)
    assert window(cells, 1, 1)[0] == cells[0, 5, 5]
