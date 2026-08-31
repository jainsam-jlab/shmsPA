#!/usr/bin/env python3
"""Truth-matched pion 3/4-trigger efficiency analysis.

Usage:
    python3 CalcMissed_v2.py INPUT[.root] CENTRAL_MOMENTUM_GEV [pi+|pi-]

The pion sample is selected from Tracks (primary track: ParentID == 0), not
from detector PID. PA and Tracks are joined explicitly with EventID.
"""

import argparse
import csv
import math
import os
import sys

import ROOT


PDG_CODES = {"pi+": 211, "pi-": -211}
ENERGY_THRESHOLD = 0.5  # MeV
TIME_WINDOW = 20.0      # ns
S2Y_NPE_THRESHOLD = 100


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="input ROOT filename, with or without .root")
    parser.add_argument("central_momentum", type=float, help="central momentum in GeV")
    parser.add_argument("particle", choices=sorted(PDG_CODES), help="truth pion charge")
    return parser.parse_args()


def binomial_percent(numerator, denominator):
    if denominator == 0:
        return 0.0, 0.0
    fraction = numerator / denominator
    return 100.0 * fraction, 100.0 * math.sqrt(fraction * (1.0 - fraction) / denominator)


def trigger_multiplicity(event):
    s1x = event.S1XEnergy > ENERGY_THRESHOLD
    s1y = (event.S1YEnergy > ENERGY_THRESHOLD and
           abs(event.S1YTime - event.S1XTime) < TIME_WINDOW)
    s2x = (event.S2XEnergy > ENERGY_THRESHOLD and
           abs(event.S2XTime - event.S1XTime) < TIME_WINDOW)
    s2y = (event.S2YNPE > S2Y_NPE_THRESHOLD and
           abs(event.S2YTime - event.S1XTime) < TIME_WINDOW)
    return sum((s1x, s1y, s2x, s2y))


def main():
    args = parse_args()
    ROOT.gROOT.SetBatch(True)

    input_name = args.input if args.input.endswith(".root") else args.input + ".root"
    root_file = ROOT.TFile.Open(input_name, "READ")
    if not root_file or root_file.IsZombie():
        sys.exit(f"Error: could not open {input_name}")

    pa = root_file.Get("PA")
    tracks = root_file.Get("Tracks")
    if not pa or not tracks:
        sys.exit("Error: input must contain both PA and Tracks trees")

    # Build the truth pion event set. There is one primary track per event;
    # TrackID is deliberately not assumed to equal 1.
    wanted_pdg = PDG_CODES[args.particle]
    pion_event_ids = set()
    duplicate_primary_ids = set()
    for track in tracks:
        if int(track.ParentID) == 0 and int(track.PDG) == wanted_pdg:
            event_id = int(track.EventID)
            if event_id in pion_event_ids:
                duplicate_primary_ids.add(event_id)
            pion_event_ids.add(event_id)

    hist_multiplicity = ROOT.TH1D(
        "hodo_planes_fired_truth_pion",
        "Truth primary pions;Number of hodoscope planes fired;Events",
        5, -0.5, 4.5)

    pa_event_ids = set()
    n_trigger = 0
    for event in pa:
        event_id = int(event.EventID)
        pa_event_ids.add(event_id)
        if event_id not in pion_event_ids:
            continue

        n_fired = trigger_multiplicity(event)
        hist_multiplicity.Fill(n_fired)
        if n_fired >= 3:
            n_trigger += 1

    matched_pions = len(pion_event_ids & pa_event_ids)
    missing_pa = pion_event_ids - pa_event_ids
    trigger_pct, trigger_err = binomial_percent(n_trigger, matched_pions)
    missed = matched_pions - n_trigger
    missed_pct, missed_err = binomial_percent(missed, matched_pions)

    print("=" * 66)
    print(f"Input file:                         {input_name}")
    print(f"Truth selection:                    ParentID == 0, PDG == {wanted_pdg}")
    print(f"Truth primary pion events:          {len(pion_event_ids)}")
    print(f"Truth pions matched to PA:           {matched_pions}")
    print(f"Passed 3/4 trigger:                  {n_trigger} ({trigger_pct:.3f} +/- {trigger_err:.3f}%)")
    print(f"Missed 3/4 trigger:                  {missed} ({missed_pct:.3f} +/- {missed_err:.3f}%)")
    if missing_pa:
        print(f"WARNING: {len(missing_pa)} truth-pion EventIDs have no PA row")
    if duplicate_primary_ids:
        print(f"WARNING: {len(duplicate_primary_ids)} events contain duplicate matching primary tracks")
    print("=" * 66)

    stem = os.path.splitext(input_name)[0]
    csv_name = f"{stem}_calmissed_v2_{args.central_momentum:g}GeV.csv"
    with open(csv_name, "w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow([
            "input", "particle", "central_momentum_GeV", "truth_primary_pions",
            "matched_PA_events", "pass_3of4", "missed_3of4",
            "trigger_efficiency_percent", "trigger_efficiency_error_percent",
            "missed_trigger_percent", "missed_trigger_error_percent",
        ])
        writer.writerow([
            input_name, args.particle, args.central_momentum, len(pion_event_ids),
            matched_pions, n_trigger, missed, trigger_pct, trigger_err,
            missed_pct, missed_err,
        ])

    canvas = ROOT.TCanvas("c", "Truth pion 3/4 trigger", 700, 500)
    hist_multiplicity.SetLineColor(ROOT.kBlue)
    hist_multiplicity.SetLineWidth(2)
    hist_multiplicity.Draw("HIST")
    pdf_name = f"{stem}_calmissed_v2_{args.central_momentum:g}GeV.pdf"
    canvas.SaveAs(pdf_name)

    hist_name = f"{stem}_calmissed_v2_{args.central_momentum:g}GeV.root"
    hist_file = ROOT.TFile.Open(hist_name, "RECREATE")
    hist_multiplicity.Write()
    hist_file.Close()
    root_file.Close()
    print(f"CSV saved to:  {csv_name}")
    print(f"Plots saved to: {pdf_name}")
    print(f"Histograms:     {hist_name}")


if __name__ == "__main__":
    main()
