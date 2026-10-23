# MCFM patch

Applies to https://github.com/usarica/MCFM-7.0_JHUGen at commit 7c6db1fb79fc7e3cbf6f83fb74239cef9a33c77f.

```bash
git clone https://github.com/usarica/MCFM-7.0_JHUGen.git && cd MCFM-7.0_JHUGen
git checkout 7c6db1fb79fc7e3cbf6f83fb74239cef9a33c77f
git apply /path/to/mcfm-7.0-jhugen.patch
for d in QCDLoop/ff QCDLoop/ql TensorReduction/ov TensorReduction/pv TensorReduction/recur/small{Y,P,G,F} TensorReduction/recur; do (cd $d && mkdir -p obj && make); done
mkdir -p obj && make
```

What it changes:

- **makefiles**: `-std=legacy -fallow-argument-mismatch -w` so the Fortran 77 code builds with gfortran 16 (Apple silicon).
- **src/Parton/mrst2004.f, jeppe02.f**: declare `ff(nnx,*)` (assumed size) instead of `ff(nnx,mmy)`; gfortran 16 rejects passing smaller arrays to the fixed-size dummy.
- **src/User/usercode_f77.f, src/Need/mcfm_init.f**: new `read_anom_env` routine, called after the input card is read. It sets the anomalous HZZ couplings `ghz2` (a2) and `ghz4` (a3) from the `MCFM_GHZ2` / `MCFM_GHZ4` environment variables, switches on `AllowAnomalousCouplings`, and sets `kappa_top = kappa_bot = 1` (their defaults are 0, which would remove the SM gg->H loop once anomalous couplings are on). With `MCFM_GHZ2=0` the result reproduces the SM run exactly (-0.0723075 fb for process 130).
