"""Publication figures for modules A-E drawn with ROOT in the ATLAS style.

Reads results/*.json (written by run_module{A..E}.py) and writes PNG and PDF files to
results/figures/. Run after the module scripts:
    .venv/bin/python scripts/plots_root.py
"""
import json
import pathlib

import numpy as np
import ROOT

from hwidth.offshell import Assumptions, OffshellTemplates, offshell_expectation

ROOT.gROOT.SetBatch(True)

BASE = pathlib.Path(__file__).resolve().parents[1]
RES = BASE / "results"
FIG = RES / "figures"

LUMI_TEXT = "#sqrt{s} = 13 TeV, 36.6 fb^{-1}, education release"
SIM_TEXT = "Geant4 11.4, FTFP_BERT, 2 X_{0} upstream"
COLORS = {"higgs": "#00b4d8", "zz": "#ff6b6b", "other": "#6b59d3", "e-": "#00b4d8", "gamma": "#ffb703",
          "pi0": "#ff6b6b", "pi-": "#6b59d3", "grey": "#888888"}
_keep = []  # ROOT objects must stay alive until the canvas is saved


def color(name):
    return ROOT.TColor.GetColor(COLORS.get(name, name))


def keep(*objs):
    _keep.extend(objs)
    return objs[0]


def atlas_style():
    s = ROOT.TStyle("ATLAS", "ATLAS style")
    for f in (s.SetFrameBorderMode, s.SetCanvasBorderMode, s.SetPadBorderMode):
        f(0)
    for f in (s.SetFrameFillColor, s.SetCanvasColor, s.SetPadColor, s.SetStatColor):
        f(0)
    s.SetPadTopMargin(0.05)
    s.SetPadRightMargin(0.05)
    s.SetPadBottomMargin(0.16)
    s.SetPadLeftMargin(0.16)
    s.SetTitleXOffset(1.4)
    s.SetTitleYOffset(1.4)
    for axis in "xyz":
        s.SetLabelFont(42, axis)
        s.SetTitleFont(42, axis)
        s.SetLabelSize(0.05, axis)
        s.SetTitleSize(0.05, axis)
    s.SetTextFont(42)
    s.SetTextSize(0.05)
    s.SetMarkerStyle(20)
    s.SetMarkerSize(1.2)
    s.SetHistLineWidth(2)
    s.SetLineStyleString(2, "[12 12]")
    s.SetEndErrorSize(0.0)
    s.SetOptTitle(0)
    s.SetOptStat(0)
    s.SetOptFit(0)
    s.SetPadTickX(1)
    s.SetPadTickY(1)
    s.SetLegendBorderSize(0)
    s.SetLegendFillColor(0)
    s.SetImageScaling(2.5)
    ROOT.gROOT.SetStyle("ATLAS")
    ROOT.gROOT.ForceStyle()


def label(x, y, title="ATLAS Open Data", sub=LUMI_TEXT, extra=None, size=0.04):
    t = keep(ROOT.TLatex())
    t.SetNDC()
    t.SetTextFont(62)
    t.SetTextSize(size * 1.2)
    t.DrawLatex(x, y, title)
    t2 = keep(ROOT.TLatex())
    t2.SetNDC()
    t2.SetTextFont(42)
    t2.SetTextSize(size * 0.85)
    lines = [sub] + ([extra] if isinstance(extra, str) else list(extra or []))
    for i, line in enumerate(lines):
        if line:
            t2.DrawLatex(x, y - 1.45 * size * (i + 1), line)


def legend(x1, y1, x2, y2, size=0.036):
    leg = keep(ROOT.TLegend(x1, y1, x2, y2))
    leg.SetTextFont(42)
    leg.SetTextSize(size)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    return leg


def canvas(name, w=800, h=600):
    return keep(ROOT.TCanvas(name, name, w, h))


def save(c, name):
    c.SaveAs(str(FIG / f"{name}.png"))
    c.SaveAs(str(FIG / f"{name}.pdf"))
    _keep.clear()


def th1(name, edges, values):
    edges = np.asarray(edges, dtype="d")
    h = keep(ROOT.TH1D(name, "", len(edges) - 1, edges))
    for i, v in enumerate(values):
        h.SetBinContent(i + 1, v)
    return h


def data_graph(edges, counts, per_gev=False):
    """Data points with sqrt(N) errors at bin centres, skipping empty bins."""
    edges = np.asarray(edges, dtype=float)
    centres, widths = 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)
    g = keep(ROOT.TGraphErrors())
    for x, w, n in zip(centres, widths, counts):
        if n > 0:
            scale = 1.0 / w if per_gev else 1.0
            i = g.GetN()
            g.SetPoint(i, x, n * scale)
            g.SetPointError(i, 0.0, np.sqrt(n) * scale)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.3)
    g.SetLineWidth(4)
    return g


def graph(x, y, col, width=4, style=1):
    g = keep(ROOT.TGraph(len(x), np.asarray(x, dtype="d"), np.asarray(y, dtype="d")))
    g.SetLineColor(col)
    g.SetLineWidth(width)
    g.SetLineStyle(style)
    return g


def frame(xmin, xmax, ymin, ymax, xtitle, ytitle, logy=False):
    pad = ROOT.gPad
    if logy:
        pad.SetLogy()
    h = keep(pad.DrawFrame(xmin, ymin, xmax, ymax))
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetTitleOffset(1.35)
    return h


def nice_max(v, fraction=1.3):
    return float(np.ceil(v * fraction / 5.0) * 5.0) if v > 5 else v * fraction


# --------------------------------------------------------------------------------------------- A
STACK = [("other", "Z+jets, t#bar{t}, VVV"), ("zz", "ZZ*"), ("higgs", "H #rightarrow ZZ* (m_{H}=125 GeV)")]


def stack_plot(hist, name, logy=False):
    edges = np.asarray(hist["edges"], dtype=float)
    width = edges[1] - edges[0]
    c = canvas(name)
    if logy:
        c.SetLogy()
    stack = keep(ROOT.THStack(f"st_{name}", ""))
    hs = {}
    for key, _ in STACK:
        h = th1(f"{name}_{key}", edges, np.maximum(hist[key], 0.0) if logy else hist[key])
        h.SetFillColor(color(key))
        h.SetLineColor(color(key))
        h.SetLineWidth(0)
        stack.Add(h)
        hs[key] = h
    data = np.asarray(hist["data"])
    ymax = max((data + np.sqrt(data)).max(), stack.GetMaximum())
    if logy:
        f = frame(edges[0], edges[-1], 0.05, ymax * 300, "m_{4#font[12]{l}} [GeV]", f"Events / {width:g} GeV")
        c.SetLogy()
    else:
        f = frame(edges[0], edges[-1], 0, nice_max(ymax), "m_{4#font[12]{l}} [GeV]", f"Events / {width:g} GeV")
    stack.Draw("hist same")
    g = data_graph(edges, data)
    g.Draw("P same")
    f.Draw("axis same")
    leg = legend(0.6, 0.7, 0.93, 0.92)
    for key, text in reversed(STACK):
        leg.AddEntry(hs[key], text, "f")
    leg.AddEntry(g, "Data", "pe")
    leg.Draw()
    label(0.2, 0.87)
    save(c, name)


def module_a():
    a = json.loads((RES / "moduleA.json").read_text())
    stack_plot(a["hist_m4l"], "moduleA_m4l")
    stack_plot(a["hist_offshell"], "moduleA_m4l_offshell", logy=True)

    fit = a["fit"]
    c = canvas("moduleA_scan", 800, 700)
    c.SetTopMargin(0.22)
    mu, q = np.array(fit["scan_mu"]), np.array(fit["scan_2dnll"])
    frame(0, 2.5, 0, 10, "#mu_{on}", "-2 #Delta ln L")
    for lvl, txt in ((1.0, "68% CL"), (3.84, "95% CL")):
        line = keep(ROOT.TLine(0, lvl, 2.5, lvl))
        line.SetLineColor(color("grey"))
        line.SetLineStyle(3)
        line.SetLineWidth(2)
        line.Draw()
        t = keep(ROOT.TLatex(2.12, lvl + 0.15, txt))
        t.SetTextSize(0.035)
        t.SetTextColor(color("grey"))
        t.Draw()
    graph(mu, q, ROOT.kBlack).Draw("L same")
    label(0.16, 0.94, extra=f"#hat{{#mu}}_{{on}} = {fit['mu_hat']:.2f}^{{+{fit['mu_hi'] - fit['mu_hat']:.2f}}}"
          f"_{{-{fit['mu_hat'] - fit['mu_lo']:.2f}}},  Z = {fit['Z_obs']:.1f}#sigma ({fit['Z_exp']:.1f}#sigma exp.)")
    save(c, "moduleA_scan")


# --------------------------------------------------------------------------------------------- B
def module_b():
    b = json.loads((RES / "moduleB.json").read_text())
    edges = np.array(b["edges"], dtype=float)
    widths = np.diff(edges)
    t = b["templates"]
    off = OffshellTemplates(edges, *(np.array(t[k]) for k in ("S", "I", "B", "Q", "O")))
    a, best = Assumptions(), b["best_fit"]
    curves = [
        (offshell_expectation(off, best["mu_on"], best["mu_qq"], 0.0, a), "SM (#Gamma_{H} = #Gamma_{H}^{SM})", "higgs", 1),
        (offshell_expectation(off, 10 * best["mu_on"], best["mu_qq"], 0.0, a), "#Gamma_{H} = 10 #Gamma_{H}^{SM}", "zz", 2),
        (offshell_expectation(off, 0.0, best["mu_qq"], 0.0, a), "no off-shell Higgs", "grey", 3),
    ]

    c = canvas("moduleB_offshell", 800, 800)
    top = keep(ROOT.TPad("top", "", 0, 0.3, 1, 1))
    bot = keep(ROOT.TPad("bot", "", 0, 0, 1, 0.3))
    top.SetBottomMargin(0.02)
    top.SetTopMargin(0.06)
    bot.SetTopMargin(0.03)
    bot.SetBottomMargin(0.38)
    top.Draw()
    bot.Draw()
    top.cd()
    top.SetLogy()
    f = frame(edges[0], edges[-1], 2e-3, 100, "", "Events / GeV")
    f.GetXaxis().SetLabelSize(0)
    f.GetYaxis().SetTitleSize(0.065)
    f.GetYaxis().SetLabelSize(0.06)
    f.GetYaxis().SetTitleOffset(1.1)
    leg = legend(0.2, 0.06, 0.55, 0.32, 0.045)
    hs = []
    for y, text, col, style in curves:
        h = th1(f"b_{col}", edges, y / widths)
        h.SetLineColor(color(col))
        h.SetLineStyle(style)
        h.SetLineWidth(4)
        h.Draw("hist same")
        hs.append(h)
        leg.AddEntry(h, text, "l")
    g = data_graph(edges, b["data"], per_gev=True)
    g.Draw("P same")
    leg.AddEntry(g, "Data", "pe")
    leg.Draw()
    label(0.42, 0.85, sub=LUMI_TEXT, extra="+ MCFM gg #rightarrow ZZ (K = 1.7)", size=0.048)
    bot.cd()
    ratio = curves[1][0] / curves[0][0]
    fr = frame(edges[0], edges[-1], 0.95, 1.3, "m_{4#font[12]{l}} [GeV]", "10 #Gamma^{SM} / SM")
    for ax in (fr.GetXaxis(), fr.GetYaxis()):
        ax.SetTitleSize(0.12)
        ax.SetLabelSize(0.12)
    fr.GetYaxis().SetTitleOffset(0.55)
    fr.GetYaxis().SetNdivisions(504)
    fr.GetXaxis().SetTitleOffset(1.2)
    one = keep(ROOT.TLine(edges[0], 1, edges[-1], 1))
    one.SetLineStyle(2)
    one.Draw()
    hr = th1("b_ratio", edges, ratio)
    hr.SetLineColor(color("zz"))
    hr.SetLineWidth(4)
    hr.Draw("hist same")
    save(c, "moduleB_offshell")

    c = canvas("moduleB_scan", 800, 700)
    c.SetTopMargin(0.22)
    grid = np.array(b["R_grid"])
    frame(0, 25, 0, 8, "#Gamma_{H} / #Gamma_{H}^{SM}", "-2 #Delta ln L")
    for lvl, txt in ((1.0, "68% CL"), (3.84, "95% CL")):
        line = keep(ROOT.TLine(0, lvl, 25, lvl))
        line.SetLineColor(color("grey"))
        line.SetLineStyle(3)
        line.SetLineWidth(2)
        line.Draw()
        tx = keep(ROOT.TLatex(21.3, lvl + 0.12, txt))
        tx.SetTextSize(0.035)
        tx.SetTextColor(color("grey"))
        tx.Draw()
    go = graph(grid, b["observed"]["scan"], ROOT.kBlack)
    ge = graph(grid, b["expected"]["scan"], ROOT.kBlack, style=2)
    go.Draw("L same")
    ge.Draw("L same")
    leg = legend(0.66, 0.8, 0.95, 0.94)
    leg.AddEntry(go, "Observed", "l")
    leg.AddEntry(ge, "Expected (SM)", "l")
    leg.Draw()
    label(0.16, 0.94, extra=f"#Gamma_{{H}} < {b['observed']['width_95_mev']:.0f} MeV "
          f"({b['expected']['width_95_mev']:.0f} MeV exp.) at 95% CL")
    save(c, "moduleB_scan")


# --------------------------------------------------------------------------------------------- C
SHORT = {"baseline": "Baseline", "K_low": "K_{gg} = 1.4", "K_high": "K_{gg} = 2.0",
         "K_bkg_30": "K_{B} #pm 30%, K_{I} = #sqrt{K_{S}K_{B}}", "qq_fixed": "#mu_{qq} constrained to 10%",
         "ff_1000": "Form factor, #Lambda = 1 TeV", "ff_500": "Form factor, #Lambda = 500 GeV"}


def module_c():
    cdat = json.loads((RES / "moduleC.json").read_text())
    sc = cdat["scenarios"]
    n = len(sc)
    c = canvas("moduleC_assumptions", 900, 600)
    c.SetLeftMargin(0.36)
    c.SetTopMargin(0.2)
    h = keep(ROOT.TH2D("c_frame", "", 10, 0, 3, n, 0, n))
    for i, s in enumerate(sc):
        h.GetYaxis().SetBinLabel(n - i, SHORT.get(s["key"], s["key"]))
    h.GetYaxis().SetLabelSize(0.045)
    h.GetXaxis().SetTitle("Fitted #Gamma_{H}/#Gamma_{H}^{SM}  (HL-LHC 3000 fb^{-1}, truth = 1)")
    h.GetXaxis().SetTitleOffset(1.3)
    h.GetXaxis().SetTitleSize(0.045)
    h.Draw("axis")
    for i, s in enumerate(sc):
        y = n - i - 0.5
        lo, hi = s["hl68"]
        box = keep(ROOT.TBox(lo, y - 0.18, min(hi, 3), y + 0.18))
        box.SetFillColor(color("higgs"))
        box.Draw()
        m = keep(ROOT.TMarker(s["hl_best"], y, 20))
        m.SetMarkerSize(1.4)
        m.Draw()
    line = keep(ROOT.TLine(1, 0, 1, n))
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()
    h.Draw("axis same")
    leg = legend(0.66, 0.81, 0.95, 0.95, 0.038)
    leg.AddEntry(box, "68% CL interval", "f")
    leg.AddEntry(m, "best fit", "p")
    leg.Draw()
    label(0.04, 0.93, sub="Asimov data with baseline model, #sqrt{s} = 13 TeV", size=0.042)
    save(c, "moduleC_assumptions")


# --------------------------------------------------------------------------------------------- D
def module_d():
    d = json.loads((RES / "moduleD.json").read_text())
    grid = np.array(d["grid"])
    styles = [("obs_comb", ROOT.kBlack, 1, "On + off shell (obs.)"), ("exp_comb", ROOT.kBlack, 2, "On + off shell (exp.)"),
              ("obs_on", color("zz"), 1, "On shell only (obs.)"), ("exp_on", color("zz"), 2, "On shell only (exp.)")]
    for name, title in (("a2", "a_{2} (CP-even)"), ("a3", "a_{3} (CP-odd)")):
        res = d["couplings"][name]
        c = canvas(f"moduleD_eft_{name}", 800, 700)
        c.SetTopMargin(0.25)
        frame(-3, 3, 0, 10, title, "-2 #Delta ln L")
        line = keep(ROOT.TLine(-3, 3.84, 3, 3.84))
        line.SetLineColor(color("grey"))
        line.SetLineStyle(3)
        line.SetLineWidth(2)
        line.Draw()
        leg = legend(0.63, 0.77, 0.95, 0.99, 0.032)
        for tag, col, style, text in styles:
            g = graph(grid, res[tag]["scan"], col, style=style)
            g.Draw("L same")
            leg.AddEntry(g, text, "l")
        leg.Draw()
        label(0.16, 0.94, extra=["+ MCFM-JHUGen, #Gamma_{H} = #Gamma_{H}^{SM}"])
        save(c, f"moduleD_eft_{name}")


# --------------------------------------------------------------------------------------------- E
PARTICLES = [("e-", "e^{-}"), ("gamma", "#gamma"), ("pi0", "#pi^{0}"), ("pi-", "#pi^{-}")]
VARS = {"R_eta": "R_{#eta}", "w_eta2": "w_{#eta 2}", "ws_tot": "w_{s,tot} [strips]", "E_ratio": "E_{ratio}",
        "R_had": "R_{had} (visible)", "E_over_p": "E/p"}


def module_e():
    e = json.loads((RES / "moduleE.json").read_text())
    for var, title in VARS.items():
        dist = e["distributions_50GeV"][var]
        c = canvas(f"moduleE_{var}")
        ymax = max(max(dist[p]["density"]) for p, _ in PARTICLES if p in dist)
        edges = dist["e-"]["edges"]
        logy = var in ("R_had",)
        f = frame(edges[0], edges[-1], 1e-3 * ymax if logy else 0, ymax * (20 if logy else 1.5), title,
                  "1/N dN/dx", logy=logy)
        leg = legend(0.72, 0.62, 0.93, 0.9)
        for p, ptitle in PARTICLES:
            if p not in dist:
                continue
            h = th1(f"e_{var}_{p}", dist[p]["edges"], dist[p]["density"])
            h.SetLineColor(color(p))
            h.SetLineWidth(4)
            h.Draw("hist same")
            leg.AddEntry(h, ptitle, "l")
        f.Draw("axis same")
        leg.Draw()
        label(0.2, 0.87, title="g4calo simulation", sub=SIM_TEXT, extra="50 GeV, under/overflow in the edge bins")
        save(c, f"moduleE_{var}")

    c = canvas("moduleE_profile")
    frame(0, 30, 0, 0.16, "Depth in calorimeter [X_{0}]", "Fraction of energy / X_{0}")
    leg = legend(0.66, 0.55, 0.93, 0.9)
    for p, ptitle in PARTICLES:
        prof = e["profile"].get(p)
        if not prof:
            continue
        h = th1(f"prof_{p}", np.arange(31, dtype=float), prof["dEdt"])
        h.SetLineColor(color(p))
        h.SetLineWidth(4)
        h.Draw("hist same")
        leg.AddEntry(h, ptitle, "l")
        if p in ("e-", "gamma"):
            fn = keep(ROOT.TF1(f"gam_{p}", "[2]*[1]*pow([1]*x,[0]-1)*exp(-[1]*x)/TMath::Gamma([0])", 0.5, 29.5))
            fn.SetParameters(4.0, 0.5, 1.0)
            h.Fit(fn, "QNR0")
            fn.SetLineColor(color(p))
            fn.SetLineStyle(2)
            fn.SetLineWidth(2)
            fn.Draw("same")
            tmax = (fn.GetParameter(0) - 1) / fn.GetParameter(1)
            leg.AddEntry(fn, f"#Gamma fit, t_{{max}} = {tmax:.1f} X_{{0}}", "l")
    leg.Draw()
    label(0.2, 0.87, title="g4calo simulation", sub=SIM_TEXT, extra="50 GeV")
    save(c, "moduleE_profile")

    rocs = {"gamma_vs_pi0": ("#gamma efficiency", "#pi^{0} rejection", ["E_ratio", "ws_tot", "w_eta2", "R_eta"]),
            "e-_vs_pi-": ("e^{-} efficiency", "#pi^{-} rejection", ["R_had", "E_over_p", "R_eta", "w_eta2"])}
    palette = ["#00b4d8", "#ff6b6b", "#6b59d3", "#ffb703"]
    for key, (xt, yt, variables) in rocs.items():
        r = e["roc"][key]
        c = canvas(f"moduleE_roc_{key.replace('-', 'm')}")
        pim = key == "e-_vs_pi-"
        frame(0.5, 1.0, 1, 5000 if pim else 100, xt, yt, logy=True)
        leg = legend(0.2, 0.36, 0.5, 0.6) if pim else legend(0.62, 0.68, 0.93, 0.92)
        for var, col in zip(variables, palette):
            es, eb = np.array(r[var]["eff_sig"]), np.array(r[var]["eff_bkg"])
            ok = (es >= 0.5) & (eb > 0)
            g = graph(es[ok], 1.0 / eb[ok], color(col))
            g.Draw("L same")
            leg.AddEntry(g, VARS[var].replace(" [strips]", ""), "l")
        leg.Draw()
        label(0.45 if pim else 0.2, 0.88, title="g4calo simulation", sub=SIM_TEXT, extra="25 GeV, 1000 events per particle")
        save(c, f"moduleE_roc_{key.replace('-', 'm')}")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    atlas_style()
    for fn in (module_a, module_b, module_c, module_d, module_e):
        fn()


if __name__ == "__main__":
    main()
