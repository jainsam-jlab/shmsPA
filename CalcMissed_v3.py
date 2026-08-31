#!/usr/bin/env python3
"""Streaming pi+ transport, trigger, and material-loss analysis for shmsPA.

Usage:
    ./CalcMissed_v3.py INPUT.root CENTRAL_MOMENTUM_GEV [options]

SHMS positive polarity is assumed: every event must contain exactly one
primary pi+ (ParentID=0, PDG=211).  "Any pion" always means any pi+ track,
including descendants.  The script uses the bounded-memory TreeReader from
validate_output.py and can therefore process production-size ROOT files.
"""

import argparse
import collections
import csv
import json
import math
import os
import sys

from validate_output import TreeReader


PION_PDG = 211
FLAGS = ("ReachedS1", "ReachedHGC", "ReachedAGC", "ReachedS2",
         "ReachedCal", "ExitedCal")
STAGES = (
    ("Generated", "ReachedS1"),
    ("ReachedS1", "ReachedHGC"),
    ("ReachedHGC", "ReachedAGC"),
    ("ReachedAGC", "ReachedS2"),
    ("ReachedS2", "ReachedCal"),
    ("ReachedCal", "ExitedCal"),
)
PA_BRANCHES = (
    "EventID", "PrimaryPDG",
    "S1XTime", "S1YTime", "S2XTime", "S2YTime",
    "S1XEnergy", "S1YEnergy", "S2XEnergy", "S2YNPE",
)
TRACK_BRANCHES = (
    "EventID", "TrackID", "ParentID", "PDG", "EndProcess", "EndVolume",
) + FLAGS


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Geant4 ROOT output containing PA and Tracks")
    parser.add_argument("central_momentum", type=float, help="central SHMS momentum in GeV")
    parser.add_argument("--csv-prefix", metavar="PREFIX",
                        help="output prefix (default: derived from input and momentum)")
    parser.add_argument("--json", dest="json_path", metavar="FILE",
                        help="JSON summary filename (default: <prefix>_summary.json)")
    parser.add_argument("--max-diagnostics", type=int, default=20, metavar="N",
                        help="number of process/volume rows per loss stage (default: 20)")
    parser.add_argument("--no-root-output", action="store_true",
                        help="do not create the optional ROOT histogram file")
    return parser.parse_args()


def trigger_multiplicity(row):
    s1x = row["S1XEnergy"] > 0.5
    s1y = row["S1YEnergy"] > 0.5 and abs(row["S1YTime"] - row["S1XTime"]) < 20.0
    s2x = row["S2XEnergy"] > 0.5 and abs(row["S2XTime"] - row["S1XTime"]) < 20.0
    s2y = row["S2YNPE"] > 100 and abs(row["S2YTime"] - row["S1XTime"]) < 20.0
    return int(s1x) + int(s1y) + int(s2x) + int(s2y)


def binomial(numerator, denominator):
    if denominator == 0:
        return 0.0, 0.0
    fraction = numerator / denominator
    error = math.sqrt(fraction * (1.0 - fraction) / denominator)
    return fraction, error


def metric(name, definition, numerator, denominator):
    fraction, error = binomial(numerator, denominator)
    return {
        "name": name,
        "definition": definition,
        "numerator": numerator,
        "denominator": denominator,
        "fraction": fraction,
        "percent": 100.0 * fraction,
        "error_percent": 100.0 * error,
    }


def require_schema(reader, tree_name, required):
    found = reader.branches(tree_name)
    if found is None:
        raise RuntimeError(f"input does not contain the {tree_name} tree")
    missing = sorted(set(required) - found)
    if missing:
        raise RuntimeError(f"{tree_name} is missing branches: {', '.join(missing)}")


def default_prefix(input_name, momentum):
    stem = os.path.splitext(input_name)[0]
    return f"{stem}_calmissed_v3_{momentum:g}GeV"


def analyze(args):
    reader = TreeReader(args.input)
    try:
        require_schema(reader, "PA", PA_BRANCHES)
        require_schema(reader, "Tracks", TRACK_BRANCHES)

        pa_events = set()
        duplicate_pa = 0
        pa_primary_pdg_zero = 0
        pa_primary_pdg_other = 0
        trigger_pass = set()
        trigger_counts = collections.Counter()

        for row in reader.rows("PA", PA_BRANCHES):
            event_id = int(row["EventID"])
            if event_id in pa_events:
                duplicate_pa += 1
            pa_events.add(event_id)
            pdg = int(row["PrimaryPDG"])
            if pdg == 0:
                pa_primary_pdg_zero += 1
            elif pdg != PION_PDG:
                pa_primary_pdg_other += 1
            multiplicity = trigger_multiplicity(row)
            trigger_counts[multiplicity] += 1
            if multiplicity >= 3:
                trigger_pass.add(event_id)

        track_events = set()
        primary_events = set()
        primary_flag_events = {flag: set() for flag in FLAGS}
        any_flag_events = {flag: set() for flag in FLAGS}
        primary_end = {}
        duplicate_primary_events = set()
        unexpected_primary_events = set()
        noncontiguous_events = set()

        current_event = None
        primary_count = 0
        expected_primary_count = 0

        def finish_event():
            if current_event is None:
                return
            if primary_count != 1 or expected_primary_count != 1:
                duplicate_primary_events.add(current_event)

        for row in reader.rows("Tracks", TRACK_BRANCHES):
            event_id = int(row["EventID"])
            if event_id != current_event:
                finish_event()
                if event_id in track_events:
                    noncontiguous_events.add(event_id)
                track_events.add(event_id)
                current_event = event_id
                primary_count = 0
                expected_primary_count = 0

            parent_id = int(row["ParentID"])
            pdg = int(row["PDG"])
            if parent_id == 0:
                primary_count += 1
                if pdg == PION_PDG:
                    expected_primary_count += 1
                    primary_events.add(event_id)
                    primary_end[event_id] = (str(row["EndProcess"]), str(row["EndVolume"]))
                    for flag in FLAGS:
                        if int(row[flag]):
                            primary_flag_events[flag].add(event_id)
                else:
                    unexpected_primary_events.add(event_id)

            if pdg == PION_PDG:
                for flag in FLAGS:
                    if int(row[flag]):
                        any_flag_events[flag].add(event_id)

        finish_event()

        generated = len(primary_events)
        matched_events = primary_events & pa_events
        efficiencies = []
        for definition, flag_sets in (("primary_pi+", primary_flag_events),
                                      ("any_pi+", any_flag_events)):
            for flag in FLAGS:
                numerator = len(flag_sets[flag] & primary_events)
                efficiencies.append(metric(flag, definition, numerator, generated))

        stage_metrics = []
        diagnostic_counts = collections.defaultdict(collections.Counter)
        for definition, flag_sets in (("primary_pi+", primary_flag_events),
                                      ("any_pi+", any_flag_events)):
            for upstream, downstream in STAGES:
                upstream_events = primary_events if upstream == "Generated" else flag_sets[upstream]
                downstream_events = upstream_events & flag_sets[downstream]
                result = metric(f"{downstream}|{upstream}", definition,
                                len(downstream_events), len(upstream_events))
                result["loss_fraction"] = 1.0 - result["fraction"] if result["denominator"] else 0.0
                result["loss_percent"] = 100.0 * result["loss_fraction"]
                stage_metrics.append(result)

                # End diagnostics describe the primary track even for the
                # any-pion definition; this avoids calling secondary recovery
                # itself an absorption process.
                for event_id in upstream_events - downstream_events:
                    process, volume = primary_end.get(event_id, ("Unknown", "Unknown"))
                    diagnostic_counts[(definition, f"{upstream}->{downstream}")][
                        (process, volume)] += 1

        trigger_metrics = []
        matched_trigger = trigger_pass & matched_events
        trigger_metrics.append(metric("3_of_4_trigger", "primary_pi+",
                                      len(matched_trigger), len(matched_events)))
        cross_tables = {}
        for flag in ("ReachedS2", "ReachedCal", "ExitedCal"):
            reached = primary_flag_events[flag] & matched_events
            table = {
                "trigger_pass_flag_true": len(matched_trigger & reached),
                "trigger_pass_flag_false": len(matched_trigger - reached),
                "trigger_fail_flag_true": len(reached - matched_trigger),
                "trigger_fail_flag_false": len(matched_events - matched_trigger - reached),
            }
            cross_tables[flag] = table

        issues = {
            "duplicate_PA_rows": duplicate_pa,
            "PA_PrimaryPDG_zero": pa_primary_pdg_zero,
            "PA_PrimaryPDG_other": pa_primary_pdg_other,
            "events_without_exactly_one_pi+_primary": len(duplicate_primary_events),
            "events_with_unexpected_primary_PDG": len(unexpected_primary_events),
            "noncontiguous_Tracks_event_blocks": len(noncontiguous_events),
            "PA_events_without_primary_track": len(pa_events - primary_events),
            "primary_events_without_PA_row": len(primary_events - pa_events),
        }

        ordering_issues = {}
        for definition, flag_sets in (("primary_pi+", primary_flag_events),
                                      ("any_pi+", any_flag_events)):
            for upstream, downstream in STAGES:
                if upstream == "Generated":
                    continue
                bad = flag_sets[downstream] - flag_sets[upstream]
                ordering_issues[f"{definition}:{downstream}_without_{upstream}"] = len(bad)

        diagnostics = []
        for (definition, transition), counts in diagnostic_counts.items():
            total = sum(counts.values())
            for (process, volume), count in counts.most_common(args.max_diagnostics):
                diagnostics.append({
                    "definition": definition,
                    "transition": transition,
                    "end_process": process,
                    "end_volume": volume,
                    "count": count,
                    "percent_of_transition_losses": 100.0 * count / total if total else 0.0,
                })

        return {
            "input": os.path.abspath(args.input),
            "central_momentum_GeV": args.central_momentum,
            "reader": reader.backend,
            "particle": "pi+",
            "PDG": PION_PDG,
            "counts": {
                "PA_events": len(pa_events),
                "Tracks_events": len(track_events),
                "generated_primary_pi+": generated,
                "matched_PA_primary_pi+": len(matched_events),
                "trigger_multiplicity": {str(i): trigger_counts[i] for i in range(5)},
            },
            "issues": issues,
            "ordering_issues": ordering_issues,
            "efficiencies": efficiencies,
            "stage_efficiencies": stage_metrics,
            "trigger_efficiencies": trigger_metrics,
            "trigger_cross_tables": cross_tables,
            "loss_diagnostics": diagnostics,
        }
    finally:
        reader.close()


def write_metric_csv(path, rows):
    fields = ("name", "definition", "numerator", "denominator", "fraction",
              "percent", "error_percent", "loss_fraction", "loss_percent")
    with open(path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_diagnostics_csv(path, rows):
    fields = ("definition", "transition", "end_process", "end_volume", "count",
              "percent_of_transition_losses")
    with open(path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_root_histograms(path, result):
    try:
        import ROOT
    except ImportError:
        return False
    ROOT.gROOT.SetBatch(True)
    output = ROOT.TFile.Open(path, "RECREATE")
    histogram = ROOT.TH1D("hodo_planes_fired_primary_pi_plus",
                          "Primary pi+;Hodoscope planes fired;Events", 5, -0.5, 4.5)
    for multiplicity, count in result["counts"]["trigger_multiplicity"].items():
        histogram.SetBinContent(int(multiplicity) + 1, count)
    histogram.Write()
    output.Close()
    return True


def print_summary(result):
    counts = result["counts"]
    print("=" * 78)
    print(f"Input:                         {result['input']}")
    print(f"Reader:                        {result['reader']}")
    print(f"Generated primary pi+:         {counts['generated_primary_pi+']}")
    print(f"Matched primary pi+ / PA:      {counts['matched_PA_primary_pi+']}")
    print("-" * 78)
    print(f"{'Definition':<14} {'Selection':<25} {'Count':>12} {'Percent':>12}")
    for row in result["efficiencies"] + result["trigger_efficiencies"]:
        count = f"{row['numerator']}/{row['denominator']}"
        value = f"{row['percent']:.4f} +/- {row['error_percent']:.4f}"
        print(f"{row['definition']:<14} {row['name']:<25} {count:>12} {value:>12}")
    print("-" * 78)
    nonzero_issues = {key: value for key, value in result["issues"].items() if value}
    nonzero_ordering = {key: value for key, value in result["ordering_issues"].items() if value}
    if nonzero_issues or nonzero_ordering:
        print("VALIDATION WARNINGS")
        for key, value in {**nonzero_issues, **nonzero_ordering}.items():
            print(f"  {key}: {value}")
    else:
        print("Internal validation checks: PASS")
    print("=" * 78)


def main():
    args = parse_args()
    if args.max_diagnostics < 0:
        sys.exit("Error: --max-diagnostics must be nonnegative")
    if not os.path.isfile(args.input):
        sys.exit(f"Error: input file does not exist: {args.input}")
    prefix = args.csv_prefix or default_prefix(args.input, args.central_momentum)
    json_path = args.json_path or prefix + "_summary.json"

    try:
        result = analyze(args)
    except (OSError, RuntimeError) as exc:
        sys.exit(f"Error: {exc}")

    efficiency_csv = prefix + "_efficiencies.csv"
    stage_csv = prefix + "_stages.csv"
    diagnostics_csv = prefix + "_loss_diagnostics.csv"
    write_metric_csv(efficiency_csv,
                     result["efficiencies"] + result["trigger_efficiencies"])
    write_metric_csv(stage_csv, result["stage_efficiencies"])
    write_diagnostics_csv(diagnostics_csv, result["loss_diagnostics"])
    with open(json_path, "w", encoding="utf-8") as output:
        json.dump(result, output, indent=2, sort_keys=True)
        output.write("\n")

    root_path = prefix + "_histograms.root"
    wrote_root = False if args.no_root_output else write_root_histograms(root_path, result)
    print_summary(result)
    print(f"Efficiencies CSV: {efficiency_csv}")
    print(f"Stage CSV:        {stage_csv}")
    print(f"Diagnostics CSV:  {diagnostics_csv}")
    print(f"JSON summary:     {json_path}")
    if wrote_root:
        print(f"ROOT histograms:  {root_path}")
    elif not args.no_root_output:
        print("ROOT histograms:  skipped (PyROOT unavailable)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
