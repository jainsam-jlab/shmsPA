#!/usr/bin/env python3
"""Classify secondary S2 response when the primary pi+ did not reach S2."""

import argparse
import collections
import csv
import json
import os
import sys

from validate_output import TreeReader


PION_PDG = 211
PA_BRANCHES = (
    "EventID", "S1XTime", "S1YTime", "S2XTime", "S2YTime",
    "S1XEnergy", "S1YEnergy", "S2XEnergy", "S2YNPE",
)
TRACK_BRANCHES = ("EventID", "TrackID", "ParentID", "PDG", "ReachedS2")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Geant4 ROOT output containing PA and Tracks")
    parser.add_argument("--output-prefix", default="secondary_s2",
                        help="output prefix (default: secondary_s2)")
    return parser.parse_args()


def trigger_multiplicity(row):
    s1x = row["S1XEnergy"] > 0.5
    s1y = row["S1YEnergy"] > 0.5 and abs(row["S1YTime"] - row["S1XTime"]) < 20.0
    s2x = row["S2XEnergy"] > 0.5 and abs(row["S2XTime"] - row["S1XTime"]) < 20.0
    s2y = row["S2YNPE"] > 100 and abs(row["S2YTime"] - row["S1XTime"]) < 20.0
    return int(s1x) + int(s1y) + int(s2x) + int(s2y)


def require_schema(reader, tree, required):
    found = reader.branches(tree)
    if found is None:
        raise RuntimeError(f"missing {tree} tree")
    missing = sorted(set(required) - found)
    if missing:
        raise RuntimeError(f"{tree} missing branches: {', '.join(missing)}")


def analyze(input_path):
    reader = TreeReader(input_path)
    try:
        require_schema(reader, "PA", PA_BRANCHES)
        require_schema(reader, "Tracks", TRACK_BRANCHES)

        trigger_pass = set()
        for row in reader.rows("PA", PA_BRANCHES):
            if trigger_multiplicity(row) >= 3:
                trigger_pass.add(int(row["EventID"]))

        categories = collections.Counter()
        species_event_counts = collections.Counter()
        records = []
        current_event = None
        primary_pi_found = False
        primary_reached_s2 = False
        secondary_pi_reached_s2 = False
        secondary_species = set()

        def finish_event():
            if current_event is None or current_event not in trigger_pass:
                return
            if not primary_pi_found:
                categories["missing_primary_pi+"] += 1
                return
            if primary_reached_s2:
                categories["primary_pi+_reached_S2"] += 1
                return

            if secondary_pi_reached_s2:
                category = "secondary_pi+_reached_S2"
            elif secondary_species:
                category = "other_secondary_reached_S2"
            else:
                category = "no_recorded_track_reached_S2"
            categories[category] += 1
            for pdg in secondary_species:
                species_event_counts[pdg] += 1
            records.append({
                "EventID": current_event,
                "category": category,
                "secondary_pi+_reached_S2": int(secondary_pi_reached_s2),
                "secondary_PDGs_reaching_S2": ";".join(map(str, sorted(secondary_species))),
            })

        for row in reader.rows("Tracks", TRACK_BRANCHES):
            event_id = int(row["EventID"])
            if event_id != current_event:
                finish_event()
                current_event = event_id
                primary_pi_found = False
                primary_reached_s2 = False
                secondary_pi_reached_s2 = False
                secondary_species = set()

            # Avoid further work for the roughly 1% of events that did not
            # pass the trigger, while still preserving event boundaries.
            if event_id not in trigger_pass:
                continue
            parent_id = int(row["ParentID"])
            pdg = int(row["PDG"])
            reached_s2 = bool(int(row["ReachedS2"]))
            if parent_id == 0 and pdg == PION_PDG:
                primary_pi_found = True
                primary_reached_s2 = reached_s2
            elif parent_id != 0 and reached_s2:
                secondary_species.add(pdg)
                if pdg == PION_PDG:
                    secondary_pi_reached_s2 = True
        finish_event()

        candidates = (categories["secondary_pi+_reached_S2"] +
                      categories["other_secondary_reached_S2"] +
                      categories["no_recorded_track_reached_S2"])
        return {
            "input": os.path.abspath(input_path),
            "reader": reader.backend,
            "trigger_pass_events": len(trigger_pass),
            "trigger_pass_primary_pi+_not_reached_S2": candidates,
            "categories": dict(categories),
            "secondary_species_event_counts": {
                str(pdg): count for pdg, count in species_event_counts.most_common()
            },
            "important_note": (
                "ReachedS2 classifies track presence, not which track caused each "
                "event-level trigger signal."
            ),
            "events": records,
        }
    finally:
        reader.close()


def main():
    args = parse_args()
    if not os.path.isfile(args.input):
        sys.exit(f"Error: input file does not exist: {args.input}")
    try:
        result = analyze(args.input)
    except (OSError, RuntimeError) as exc:
        sys.exit(f"Error: {exc}")

    json_path = args.output_prefix + ".json"
    csv_path = args.output_prefix + "_events.csv"
    with open(json_path, "w", encoding="utf-8") as output:
        json.dump(result, output, indent=2, sort_keys=True)
        output.write("\n")
    with open(csv_path, "w", newline="", encoding="utf-8") as output:
        fields = ("EventID", "category", "secondary_pi+_reached_S2",
                  "secondary_PDGs_reaching_S2")
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(result["events"])

    print(f"Trigger-pass events:                         {result['trigger_pass_events']}")
    print("Trigger pass, primary pi+ did not reach S2: "
          f"{result['trigger_pass_primary_pi+_not_reached_S2']}")
    for category in ("secondary_pi+_reached_S2", "other_secondary_reached_S2",
                     "no_recorded_track_reached_S2"):
        print(f"  {category:<38} {result['categories'].get(category, 0)}")
    print(f"JSON: {json_path}")
    print(f"CSV:  {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
