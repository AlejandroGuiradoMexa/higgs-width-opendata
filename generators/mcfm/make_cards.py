"""Write MCFM-7.0 (JHUGen fork) input cards for the off-shell gg->ZZ->2e2mu templates.

Runs: 129 = interference only, 130 = |H|^2 + interference, 132 = |box|^2 (gg continuum), 81 = qq->ZZ at NLO
(used to calibrate the detector efficiency against the ATLAS Open Data qq->ZZ sample).
Scales are dynamic: mu_R = mu_F = m4l/2.
"""
import pathlib
import re
import sys

MCFM_BIN = pathlib.Path.home() / "opt/src/MCFM-7.0_JHUGen/Bin"
HERE = pathlib.Path(__file__).resolve().parent

COMMON = {
    "writeroot": ".false.", "writedat": ".true.", "writetop": ".false.",
    "sqrts in GeV": "13000d0", "hmass": "125d0",
    "scale:QCD scale choice": "0.5d0", "facscale:QCD fac_scale choice": "0.5d0", "dynamicscale": "'m(3456)'",
    "itmx1, number of iterations for pre-conditioning": "5", "ncall1": "100000",
    "itmx2, number of iterations for final run": "8", "ncall2": "400000",
    "m3456min": "220d0", "m3456max": "2000d0",
    "makecuts": ".true.", "ptlepton_min": "6d0", "|etalepton|_max": "2.5d0",
    "ptlepton(2nd+)_min": "6d0", "|etalepton(2nd+)|_max": "2.5d0", "R(lept,lept)_min": "0.1d0",
    "m34min": "12d0", "m56min": "12d0",
}
RUNS = {
    # 128 (|H|^2 alone) samples m4l with a narrow Breit-Wigner at mH and returns zero above
    # 220 GeV, so the signal is taken as S = (130) - (129).
    "ggH_I": {"nproc": "129", "part": "'lord'", "pdlabel": "'cteq6l1'"},
    "ggH_SI": {"nproc": "130", "part": "'lord'", "pdlabel": "'cteq6l1'"},
    "ggZZ_B": {"nproc": "132", "part": "'lord'", "pdlabel": "'cteq6l1'"},
    "qqZZ_NLO": {"nproc": "81", "part": "'tota'", "pdlabel": "'CT10.00'", "omitgg": ".true.",
                 "ncall2": "200000"},
}


def card(name, extra):
    s = (MCFM_BIN / "input.DAT").read_text()
    settings = {**COMMON, **extra, "runstring": f"'{name}'"}
    for key, val in settings.items():
        if key == "runstring":
            pat = r"^.*(\['runstring'\])"
        elif key == "part":
            pat = r"^.*(\[part 'lord','real' or 'virt','tota'\])"
        else:
            pat = r"^.*(\[%s\])" % re.escape(key)
        s, n = re.subn(pat, lambda m: f"{val}\t\t{m.group(1)}", s, count=1, flags=re.M)
        if n != 1:
            sys.exit(f"setting not found: {key}")
    return s


if __name__ == "__main__" and "--eft" not in sys.argv:
    for name, extra in RUNS.items():
        (HERE / f"{name}.DAT").write_text(card(name, extra))
        print(name)


# EFT templates (module D): anomalous HZZ couplings a2 (ghz2) and a3 (ghz4), set at run time through
# the MCFM_GHZ2 / MCFM_GHZ4 environment variables read by our read_anom_env patch.
# Off shell: processes 129 and 130 as above. On shell: 128 in a 105-140 GeV window gives the on-shell
# rate, which also changes with the couplings.
EFT_POINTS = {"a2p": {"MCFM_GHZ2": "0.5"}, "a2m": {"MCFM_GHZ2": "-0.5"}, "a3p": {"MCFM_GHZ4": "0.5"}}
EFT_RUNS = {
    "I": {"nproc": "129", "part": "'lord'", "pdlabel": "'cteq6l1'"},
    "SI": {"nproc": "130", "part": "'lord'", "pdlabel": "'cteq6l1'"},
    "on": {"nproc": "128", "part": "'lord'", "pdlabel": "'cteq6l1'", "m3456min": "105d0", "m3456max": "140d0",
           "ncall2": "200000"},
}


def write_eft_cards():
    eft = HERE / "eft"
    eft.mkdir(exist_ok=True)
    jobs = []
    for point, env in [("sm", {})] + list(EFT_POINTS.items()):
        for run, extra in EFT_RUNS.items():
            if point == "sm" and run != "on":
                continue  # SM off-shell templates already exist
            name = f"e_{point}_{run}"
            (eft / f"{name}.DAT").write_text(card(name, extra))
            jobs.append(" ".join([f"{k}={v}" for k, v in env.items()] + [name]))
    (eft / "jobs.txt").write_text("\n".join(jobs) + "\n")
    print("\n".join(jobs))

if __name__ == "__main__" and "--eft" in sys.argv:
    write_eft_cards()
