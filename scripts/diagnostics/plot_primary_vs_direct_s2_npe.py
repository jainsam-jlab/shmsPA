#!/usr/bin/env python3
"""Compare event-level Cherenkov NPE for primary- and daughter-S2 truth classes."""
import argparse
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="secondary-pion truth diagnostic ROOT file")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)

    source = ROOT.TFile.Open(args.input, "READ")
    tree = source.Get("EventSummary")
    primary_cut = "PrimaryReachedS2"
    daughter_cut = "!PrimaryReachedS2&&HasDirectDaughterReachedS2"
    classes = (
        ("Primary #pi^{+} reached S2", primary_cut, ROOT.kBlue + 1),
        ("Primary missed S2; direct daughter #pi^{+} reached S2",
         daughter_cut, ROOT.kRed + 1),
    )

    for variable, bins, low, high, threshold in (
            ("AGCNPE", 120, 0, 600, 20),
            ("HGCNPE", 120, 0, 3000, 50)):
        canvas = ROOT.TCanvas(f"c_{variable}", "", 1050, 760)
        legend = ROOT.TLegend(0.39, 0.70, 0.89, 0.88)
        histograms = []
        for index, (label, cut, color) in enumerate(classes):
            hist = ROOT.TH1D(f"h_{variable}_{index}",
                             f";Event-level {variable};Normalized events",
                             bins, low, high)
            tree.Draw(f"{variable}>>{hist.GetName()}", cut, "goff")
            entries = int(hist.GetEntries())
            if hist.Integral():
                hist.Scale(1.0 / hist.Integral())
            hist.SetLineColor(color)
            hist.SetLineWidth(3)
            hist.SetDirectory(0)
            histograms.append(hist)
            legend.AddEntry(hist, f"{label} (N={entries:,})", "l")
        maximum = 1.25 * max(hist.GetMaximum() for hist in histograms)
        histograms[0].SetMaximum(maximum)
        histograms[0].Draw("hist")
        histograms[1].Draw("hist same")
        line = ROOT.TLine(threshold, 0, threshold, maximum)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        legend.AddEntry(line, f"Old threshold: {variable} = {threshold}", "l")
        legend.Draw()
        canvas.SaveAs(str(output / f"{variable}_primary_vs_direct_daughter_S2.png"))

    for index, (label, cut, _) in enumerate(classes):
        canvas = ROOT.TCanvas(f"c_2d_{index}", "", 950, 780)
        title = ("Primary #pi^{+} reached S2" if index == 0
                 else "#pi^{+} direct daughter reached S2")
        hist = ROOT.TH2D(f"h_2d_{index}",
                         f"{title};Event-level HGCNPE;Event-level AGCNPE",
                         120, 0, 3000, 120, 0, 600)
        tree.Draw(f"AGCNPE:HGCNPE>>{hist.GetName()}", cut, "colz")
        hgc_line = ROOT.TLine(50, 0, 50, 600)
        agc_line = ROOT.TLine(0, 20, 3000, 20)
        for line in (agc_line, hgc_line):
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()
        suffix = "primary_reached_S2" if index == 0 else "direct_daughter_reached_S2"
        canvas.SaveAs(str(output / f"AGC_vs_HGC_{suffix}.png"))

    print(f"Primary reached S2: {tree.GetEntries(primary_cut)}")
    print(f"Primary missed S2, direct daughter reached S2: {tree.GetEntries(daughter_cut)}")
    print(f"Plots saved in {output}")


if __name__ == "__main__":
    main()
