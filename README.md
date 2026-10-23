# Higgs boson width with ATLAS Open Data

A small-scale reproduction of the four tasks of a PhD project on the Higgs boson width in H→ZZ:

| Module | Task | What is done here |
|---|---|---|
| A | (foundation) | H→ZZ*→4ℓ selection on 36.6 fb⁻¹ of 13 TeV ATLAS Open Data, pyhf fit of the on-shell signal strength |
| B | Width from H→ZZ | Joint on-shell/off-shell fit with MCFM gg→ZZ signal, interference and continuum templates |
| C | Stability tests | The width fit repeated with one assumption changed at a time (K-factors, background normalisation, new physics in the ggH loop, on/off-shell coupling equality, B_BSM in the indirect method) |
| D | EFT interpretation | Anomalous HZZ couplings a₂ (CP-even) and a₃ (CP-odd) from on- and off-shell rates |
| E | Electron shower shapes | Geant4 (C++) ATLAS-like LAr calorimeter, ATLAS shower-shape variables, material effects, e/γ/π⁰/π± separation |

Results are shown on the project page, https://higgs.alejandroguirado.com (`site/`, deployed with Cloudflare Pages).

## Main results

- **A:** μ_on = 0.93 +0.21/−0.19, observed (expected) significance 6.1σ (7.2σ).
- **B:** Γ_H/Γ_H^SM < 12.2 observed, < 13.3 expected at 95% CL (Γ_H < 50 MeV), using m₄ℓ only. Expected reach: < 7.3 at 140 fb⁻¹, < 2.6 at 3000 fb⁻¹.
- **C:** at HL-LHC luminosity, a gg K-factor of 1.4 or 2.0 instead of 1.7 moves the fitted width from 1.0 to 0.4 or 1.8; a form factor in the ggH loop at Λ = 500 GeV moves it to 2.6; freeing μ_off from μ_on leaves the width undetermined.
- **D:** |a₃| < 1.1 and a₂ > −0.52 at 95% CL from on + off shell; the on-shell rate alone gives no constraint.
- **E:** stochastic term 8.4%/√E; R_η drops from 0.969 to 0.933 and f₁ rises from 0.10 to 0.52 as upstream material grows from 0 to 6 X₀.

These are educational results from open data and simplified tools, not ATLAS results.

## Layout

```
src/hwidth/        selection, on-shell fit, off-shell model, shower shapes (Python package)
scripts/           download, one script per module, site builder
generators/mcfm/   MCFM cards, run scripts, outputs, and the patch applied to MCFM
g4calo/            Geant4 application (C++17, CMake) and production script
tests/             pytest unit tests
results/           JSON results; figures/ (ROOT, ATLAS style) and figures/mpl/ (matplotlib quick-look)
site/              one-page website (static HTML + Plotly)
```

## Reproduce

Requirements: Python ≥ 3.11, ROOT ≥ 6.30 with PyROOT, Ghostscript and poppler (figures and their web copies), gfortran, CMake, Geant4 ≥ 11 (built with `GEANT4_INSTALL_DATA=ON`), and MCFM-7.0_JHUGen with `generators/mcfm/patches/mcfm-7.0-jhugen.patch` (see the patch README) in `~/opt/src/MCFM-7.0_JHUGen`.

```bash
python -m venv .venv && .venv/bin/pip install -e . pytest
.venv/bin/python scripts/download.py                 # exactly4lep skims, ~1 GB
.venv/bin/python scripts/run_moduleA.py --reprocess
python3 generators/mcfm/make_cards.py && generators/mcfm/run_mcfm.sh
python3 generators/mcfm/make_cards.py --eft && generators/mcfm/run_eft.sh
.venv/bin/python scripts/run_moduleB.py
.venv/bin/python scripts/run_moduleC.py
.venv/bin/python scripts/run_moduleD.py
(cd g4calo && mkdir -p build && cd build && cmake .. && make) && g4calo/run_production.sh 1000
.venv/bin/python scripts/run_moduleE.py
# ROOT figures (ATLAS style) in results/figures/; PyROOT needs the Python version ROOT was built with
$(root-config --python-version | cut -d. -f1,2 | sed 's/^/python/') -m venv .venv-root && .venv-root/bin/pip install numpy iminuit && .venv-root/bin/pip install --no-deps -e .
PYTHONPATH=$(root-config --libdir) .venv-root/bin/python scripts/plots_root.py
g4calo/vis/make_views.sh                              # 3D event views (Geant4 VRML) for the page
.venv/bin/python scripts/build_site.py && (cd site && python3 -m http.server)
.venv/bin/pytest
```

## Simplifications

- Education release of ATLAS Open Data (2015–2016, 36.6 fb⁻¹); 4e + 4μ off-shell templates taken equal to 2e2μ.
- Off-shell discrimination uses m₄ℓ only (no matrix-element discriminant); no VBF off-shell component.
- gg templates at LO with a flat K-factor of 1.7; detector efficiency per m₄ℓ bin from the open-data qq→ZZ sample divided by MCFM NLO.
- The 128 process of MCFM (|H|² alone) samples a narrow Breit–Wigner at m_H and returns zero above 220 GeV, so the signal is taken as (S+I) − I from processes 130 and 129.
- The calorimeter is planar, with the impact point at a cell centre and no magnetic field.

## Data and tools

ATLAS Open Data (CC0), MCFM (Campbell, Ellis, Williams; JHUGen fork by U. Sarica et al.), Geant4, pyhf, iminuit, uproot, scikit-hep.
