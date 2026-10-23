"""Off-shell H*->ZZ->4l model and the joint on-shell/off-shell width fit.

Parton-level gg->ZZ->2e2mu templates come from MCFM (signal |S|^2, signal+interference,
continuum |B|^2). They are folded with a detector efficiency per m4l bin obtained by dividing
the reconstructed ATLAS Open Data qq->ZZ (2e2mu) yield by the MCFM NLO qq->ZZ fiducial cross
section in the same bin. 4e + 4mu are taken equal to 2e2mu, which holds well above 2 mZ.

Expected off-shell yield in bin i:
    n_i = K_gg * (mu_off * S_i + sqrt(mu_off) * I_i + B_i) + mu_qq * Q_i + O_i
with mu_off = mu_on * R under the baseline assumption that couplings are equal on and off shell,
where R = Gamma_H / Gamma_H^SM.
"""
import pathlib
import re
from dataclasses import dataclass, field

import numpy as np
from iminuit import Minuit

MCFM_BIN_WIDTH = 20.0


def read_mcfm_hist(path, hist=1):
    """Return (centres, dsigma/dm [fb/GeV], error) for histogram `hist` of an MCFM .dat file."""
    lines = pathlib.Path(path).read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if re.match(rf"\s*HIST =\s+{hist}\s*$", l))
    rows = []
    for line in lines[start + 1:]:
        if "HIST =" in line or "IS EMPTY" in line:
            break
        parts = line.split()
        if len(parts) == 3:
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2]


def mcfm_cross_section(path):
    m = re.search(r"Cross-section is:\s+([0-9.E+-]+)\s+\+/-\s+([0-9.E+-]+)", pathlib.Path(path).read_text())
    return float(m.group(1)), float(m.group(2))


def rebin_mcfm(centres, dsdm, err, edges):
    """Integrate dsigma/dm into `edges` (which must lie on MCFM bin boundaries). Returns fb per bin."""
    lo = centres - MCFM_BIN_WIDTH / 2
    out, var = np.zeros(len(edges) - 1), np.zeros(len(edges) - 1)
    for k in range(len(edges) - 1):
        sel = (lo >= edges[k] - 1e-6) & (lo + MCFM_BIN_WIDTH <= edges[k + 1] + 1e-6)
        out[k] = (dsdm[sel] * MCFM_BIN_WIDTH).sum()
        var[k] = ((err[sel] * MCFM_BIN_WIDTH) ** 2).sum()
    return out, np.sqrt(var)


@dataclass
class OffshellTemplates:
    """Expected event counts per off-shell bin (all 4l channels, before K-factors)."""
    edges: np.ndarray
    S: np.ndarray          # gg -> H* -> ZZ, |S|^2
    I: np.ndarray          # interference (negative)
    B: np.ndarray          # gg -> ZZ continuum
    Q: np.ndarray          # qq -> ZZ (reconstructed MC)
    O: np.ndarray          # reducible + triboson
    eff: np.ndarray = field(default=None)

    def scaled(self, lumi_ratio):
        return OffshellTemplates(self.edges, *(x * lumi_ratio for x in (self.S, self.I, self.B, self.Q, self.O)),
                                 eff=self.eff)


@dataclass
class OnshellTemplates:
    S: np.ndarray
    B: np.ndarray

    def scaled(self, lumi_ratio):
        return OnshellTemplates(self.S * lumi_ratio, self.B * lumi_ratio)


@dataclass
class Assumptions:
    """Switches for the stability tests (module C). Defaults reproduce the baseline analysis."""
    couple_on_off: bool = True          # mu_off = mu_on * R  (equal couplings on and off shell)
    K_gg: float = 1.7                   # NNLO/LO K-factor applied to S, I and B
    K_bkg_unc: float = 0.0              # extra relative uncertainty on the gg continuum K-factor
    K_int_mode: str = "same"            # 'same': K_I = K_gg; 'sqrt': K_I = sqrt(K_S * K_B) with K_B varied
    qq_constraint: float | None = None  # None: mu_qq free; else Gaussian width on mu_qq around 1
    form_factor_scale: float | None = None  # Lambda [GeV]: S -> S/(1+(m/Lambda)^2)^2 off shell (new physics in the loop)


def bin_centres(edges):
    return 0.5 * (edges[1:] + edges[:-1])


def offshell_expectation(t: OffshellTemplates, mu_off, mu_qq, theta_k, a: Assumptions):
    S = t.S
    I = t.I
    if a.form_factor_scale:
        ff = 1.0 / (1.0 + (bin_centres(t.edges) / a.form_factor_scale) ** 2)
        S, I = S * ff**2, I * ff
    k_b = a.K_gg * (1.0 + a.K_bkg_unc * theta_k)
    k_s = a.K_gg
    k_i = np.sqrt(k_s * k_b) if a.K_int_mode == "sqrt" else k_s
    mu_off = max(mu_off, 0.0)
    return k_s * mu_off * S + k_i * np.sqrt(mu_off) * I + k_b * t.B + mu_qq * t.Q + t.O


def poisson_nll(n_obs, n_exp):
    n_exp = np.maximum(n_exp, 1e-9)
    return float(np.sum(n_exp - n_obs * np.log(n_exp)))


class WidthFit:
    """Joint on-shell + off-shell likelihood in (mu_on, R or mu_off, mu_qq, theta_k)."""

    def __init__(self, on: OnshellTemplates, on_obs, off: OffshellTemplates, off_obs, a: Assumptions = Assumptions()):
        self.on, self.on_obs = on, np.asarray(on_obs, float)
        self.off, self.off_obs = off, np.asarray(off_obs, float)
        self.a = a

    def nll(self, mu_on, x, mu_qq, theta_k):
        """x is R = Gamma/Gamma_SM if couple_on_off, else mu_off."""
        mu_off = mu_on * x if self.a.couple_on_off else x
        n_on = mu_on * self.on.S + self.on.B
        n_off = offshell_expectation(self.off, mu_off, mu_qq, theta_k, self.a)
        val = poisson_nll(self.on_obs, n_on) + poisson_nll(self.off_obs, n_off)
        if self.a.K_bkg_unc > 0:
            val += 0.5 * theta_k**2
        if self.a.qq_constraint:
            val += 0.5 * ((mu_qq - 1.0) / self.a.qq_constraint) ** 2
        return val

    def _minuit(self, x_fixed=None):
        m = Minuit(lambda mu_on, x, mu_qq, theta_k: self.nll(mu_on, x, mu_qq, theta_k),
                   mu_on=1.0, x=1.0 if x_fixed is None else x_fixed, mu_qq=1.0, theta_k=0.0)
        m.errordef = Minuit.LIKELIHOOD
        m.limits["mu_on"] = (0.01, 5)
        m.limits["x"] = (0.0, 100)
        m.limits["mu_qq"] = (0.2, 3)
        m.limits["theta_k"] = (-5, 5)
        if self.a.K_bkg_unc == 0:
            m.fixed["theta_k"] = True
        if x_fixed is not None:
            m.fixed["x"] = True
        m.migrad()
        return m

    def best_fit(self):
        m = self._minuit()
        return {k: float(m.values[k]) for k in ("mu_on", "x", "mu_qq", "theta_k")}, 2 * float(m.fval)

    def scan(self, grid):
        _, nll0 = self.best_fit()
        prof = np.array([2 * float(self._minuit(x).fval) for x in grid])
        best = min(nll0, prof.min())
        return prof - best

    @staticmethod
    def upper_limit(grid, q, level=3.84):
        """Value where the profile -2 Delta lnL first crosses `level` above its minimum."""
        imin = int(np.argmin(q))
        right = q[imin:]
        if right.max() < level:
            return float(grid[-1])
        j = int(np.argmax(right >= level)) + imin
        return float(np.interp(level, [q[j - 1], q[j]], [grid[j - 1], grid[j]]))


def asimov(on: OnshellTemplates, off: OffshellTemplates, mu_on=1.0, R=1.0, mu_qq=1.0, a: Assumptions = Assumptions()):
    """Expected data for the given parameters (no fluctuations)."""
    n_on = mu_on * on.S + on.B
    n_off = offshell_expectation(off, mu_on * R, mu_qq, 0.0, a)
    return n_on, n_off
