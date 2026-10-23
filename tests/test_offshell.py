import numpy as np

from hwidth.offshell import (Assumptions, OffshellTemplates, OnshellTemplates, WidthFit, asimov,
                             offshell_expectation, rebin_mcfm)

EDGES = np.array([230.0, 290.0, 410.0, 1010.0])


def toy(scale=1.0):
    off = OffshellTemplates(EDGES, S=np.array([4.0, 3.0, 2.0]) * scale, I=np.array([-6.0, -5.0, -3.0]) * scale,
                            B=np.array([40.0, 30.0, 10.0]) * scale, Q=np.array([300.0, 200.0, 60.0]) * scale,
                            O=np.array([2.0, 1.0, 0.5]) * scale)
    on = OnshellTemplates(S=np.array([5.0, 20.0, 5.0]) * scale, B=np.array([3.0, 4.0, 3.0]) * scale)
    return on, off


def test_expectation_scales_with_sqrt_mu_for_interference():
    _, off = toy()
    a = Assumptions(K_gg=1.0)
    n0 = offshell_expectation(off, 0.0, 1.0, 0.0, a)
    n1 = offshell_expectation(off, 1.0, 1.0, 0.0, a)
    n4 = offshell_expectation(off, 4.0, 1.0, 0.0, a)
    assert np.allclose(n1 - n0, off.S + off.I)
    assert np.allclose(n4 - n0, 4 * off.S + 2 * off.I)


def test_fit_recovers_injected_width_on_asimov():
    on, off = toy(scale=20.0)
    a = Assumptions()
    n_on, n_off = asimov(on, off, mu_on=1.0, R=3.0, mu_qq=1.0, a=a)
    best, _ = WidthFit(on, n_on, off, n_off, a).best_fit()
    assert abs(best["x"] - 3.0) < 0.05
    assert abs(best["mu_on"] - 1.0) < 0.01


def test_decoupled_fit_measures_mu_off_directly():
    on, off = toy(scale=20.0)
    a = Assumptions(couple_on_off=False)
    n_on, n_off = asimov(on, off, mu_on=2.0, R=1.5, mu_qq=1.0, a=Assumptions())
    best, _ = WidthFit(on, n_on, off, n_off, a).best_fit()
    assert abs(best["x"] - 3.0) < 0.05  # mu_off = mu_on * R


def test_rebin_integrates_mcfm_bins():
    centres = np.array([240.0, 260.0, 280.0])
    dsdm = np.array([0.1, 0.2, 0.3])
    out, _ = rebin_mcfm(centres, dsdm, np.zeros(3), np.array([230.0, 250.0, 290.0]))
    assert np.allclose(out, [0.1 * 20, (0.2 + 0.3) * 20])
