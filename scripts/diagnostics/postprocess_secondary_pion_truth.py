#!/usr/bin/env python3
"""Create the compact CSV summary and truth-category plots."""
import argparse
import csv
import time
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--particle-label", default="pi+")
    args = parser.parse_args()

    start = time.time()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    root_file = ROOT.TFile.Open(args.input)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {args.input}")
    events = root_file.Get("EventSummary")
    daughters = root_file.Get("DirectDaughters")
    if not events or not daughters:
        raise RuntimeError("EventSummary or DirectDaughters tree is missing")

    fields = (
        "PrimaryReachedS2",
        "HasDirectDaughter",
        "DirectDaughterReachedHGC",
        "DirectDaughterReachedAGC",
        "DirectDaughterReachedS2",
        "DirectDaughterReachedCal",
        "DirectDaughterExitedCal",
    )
    rows = [{name: int(getattr(event, name)) for name in fields} for event in events]
    generated = len(rows)
    missed = [row for row in rows if not row["PrimaryReachedS2"]]

    def count(selected_rows, key):
        return sum(row[key] for row in selected_rows)

    primary_reached = count(rows, "PrimaryReachedS2")
    direct_recovered = count(missed, "DirectDaughterReachedS2")
    recovered = primary_reached + direct_recovered

    csv_rows = []

    def add(section, metric, value, denominator):
        percent = 100.0 * value / denominator if denominator else 0.0
        csv_rows.append(
            {
                "particle": args.particle_label,
                "section": section,
                "metric": metric,
                "count": value,
                "denominator": denominator,
                "percent": f"{percent:.6f}",
            }
        )

    add("all_generated_events", "generated", generated, generated)
    add("all_generated_events", "primary_reached_s2", primary_reached, generated)
    add("all_generated_events", "primary_missed_s2", len(missed), generated)
    add(
        "all_generated_events",
        "primary_or_direct_daughter_reached_s2",
        recovered,
        generated,
    )
    add("all_generated_events", "truth_not_recovered_at_s2", generated - recovered, generated)
    add("primary_missed_s2", "direct_daughter_produced", count(missed, "HasDirectDaughter"), len(missed))
    add("primary_missed_s2", "direct_daughter_reached_hgc", count(missed, "DirectDaughterReachedHGC"), len(missed))
    add("primary_missed_s2", "direct_daughter_reached_agc", count(missed, "DirectDaughterReachedAGC"), len(missed))
    add("primary_missed_s2", "direct_daughter_reached_s2", direct_recovered, len(missed))
    add("primary_missed_s2", "direct_daughter_reached_cal", count(missed, "DirectDaughterReachedCal"), len(missed))
    add("primary_missed_s2", "direct_daughter_exited_cal", count(missed, "DirectDaughterExitedCal"), len(missed))

    csv_path = output_dir / "absorption_summary.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("particle", "section", "metric", "count", "denominator", "percent"),
        )
        writer.writeheader()
        writer.writerows(csv_rows)

    categories = {
        "primary_reached_S2": "PrimaryReachedS2",
        "direct_daughter_reached_S2": (
            "!PrimaryReachedS2&&DirectDaughterReachedS2"
        ),
        "later_descendant_reached_S2": (
            "!PrimaryReachedS2&&!DirectDaughterReachedS2&&AnyDescendantReachedS2"
        ),
    }
    colors = (ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7)

    for variable, bins, low, high in (
        ("AGCNPE", 120, 0, 600),
        ("HGCNPE", 120, 0, 3000),
    ):
        canvas = ROOT.TCanvas(f"c_{variable}", f"c_{variable}", 1000, 750)
        legend = ROOT.TLegend(0.56, 0.66, 0.88, 0.88)
        histograms = []
        for index, (label, selection) in enumerate(categories.items()):
            histogram = ROOT.TH1D(
                f"h_{variable}_{index}",
                f"{args.particle_label} truth categories;Event-level {variable};Normalized events",
                bins,
                low,
                high,
            )
            events.Draw(f"{variable}>>{histogram.GetName()}", selection, "goff")
            histogram.SetDirectory(0)
            histogram.SetLineColor(colors[index])
            histogram.SetLineWidth(3)
            if histogram.Integral():
                histogram.Scale(1.0 / histogram.Integral())
            histograms.append(histogram)
            legend.AddEntry(histogram, label.replace("_", " "), "l")
        maximum = max(histogram.GetMaximum() for histogram in histograms) * 1.25
        histograms[0].SetMaximum(maximum)
        histograms[0].Draw("hist")
        for histogram in histograms[1:]:
            histogram.Draw("hist same")
        legend.Draw()
        canvas.SaveAs(str(output_dir / f"{variable}_truth_categories.png"))

    for label, selection in categories.items():
        canvas = ROOT.TCanvas(f"c2_{label}", f"c2_{label}", 900, 750)
        histogram = ROOT.TH2D(
            f"h2_{label}",
            f"{args.particle_label} {label.replace('_', ' ')};"
            "Event-level AGCNPE;Event-level HGCNPE",
            100,
            0,
            600,
            100,
            0,
            3000,
        )
        events.Draw(f"HGCNPE:AGCNPE>>{histogram.GetName()}", selection, "colz")
        canvas.SaveAs(str(output_dir / f"AGC_vs_HGC_{label}.png"))

    canvas = ROOT.TCanvas("c_stage", "c_stage", 1000, 750)
    stage_histogram = ROOT.TH1D(
        "stage",
        f"{args.particle_label} direct daughter furthest stage;Stage;Tracks",
        7,
        -0.5,
        6.5,
    )
    daughters.Draw("FurthestStage>>stage", "", "goff")
    for index, label in enumerate(
        ("BeforeS1", "S1", "HGC", "AGC", "S2", "Cal", "ExitedCal")
    ):
        stage_histogram.GetXaxis().SetBinLabel(index + 1, label)
    stage_histogram.Draw("hist")
    canvas.SaveAs(str(output_dir / "direct_daughter_furthest_stage.png"))

    print(f"CSV summary: {csv_path}")
    print(
        "Primary or same-species direct daughter reached S2: "
        f"{recovered}/{generated} "
        f"({100.0 * recovered / generated if generated else 0.0:.6f}%)"
    )
    print(f"TIMING plotting_report_s {time.time() - start:.3f}")


if __name__ == "__main__":
    main()
