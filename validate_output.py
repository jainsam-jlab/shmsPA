#!/usr/bin/env python3
"""Validate the PA and Tracks trees produced by the modified shmsPA simulation.

The checks implement the machine-testable items in the Validation Strategy of
PionAbsorption_Geant4_Analysis_Note.tex.  Either PyROOT or uproot may be used to
read the input file; no output ROOT file is created.
"""

import argparse
import collections
import json
import os
import sys


PDG_CODES = {"pi+": 211, "pi-": -211}
PA_BRANCHES = (
    "S1XTime", "S1YTime", "S2XTime", "S2YTime",
    "S1XEnergy", "S1YEnergy", "S2XEnergy", "S2YEnergy",
    "S2YNPE", "CopyNo", "AGCNPE", "HGCNPE", "NGCNPE", "CalEnergy",
    "EventID", "PrimaryPDG", "PrimaryReachedCal", "PrimaryExitedCal",
)
TRACK_BRANCHES = (
    "EventID", "TrackID", "ParentID", "PDG", "CreatorProcess",
    "StartX", "StartY", "StartZ", "StartPx", "StartPy", "StartPz", "StartKE",
    "StartVolume", "EndX", "EndY", "EndZ", "EndPx", "EndPy", "EndPz", "EndKE",
    "EndVolume", "EndProcess", "TrackLength", "ReachedS1", "ReachedHGC",
    "ReachedAGC", "ReachedS2", "ReachedCal", "ExitedCal",
)
NPE_BRANCHES = ("S2YNPE", "AGCNPE", "HGCNPE", "NGCNPE")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Geant4 ROOT output file")
    parser.add_argument("--particle", choices=sorted(PDG_CODES), default="pi+",
                        help="expected generated primary (default: pi+)")
    parser.add_argument("--kinematics", metavar="FILE",
                        help="input kinematics file; valid data lines are counted")
    parser.add_argument("--beam-on", type=int, metavar="N",
                        help="requested /run/beamOn count (for exhaustion checks)")
    parser.add_argument("--inspect-events", type=int, default=0, metavar="N",
                        help="print up to N representative event genealogies")
    parser.add_argument("--max-examples", type=int, default=10, metavar="N",
                        help="maximum IDs shown for each issue (default: 10)")
    parser.add_argument("--json", dest="json_path", metavar="FILE",
                        help="also write the full validation result as JSON")
    return parser.parse_args()


def _plain(value):
    """Convert ROOT/numpy scalar and byte-string values to ordinary Python."""
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


class TreeReader:
    def __init__(self, filename):
        self.filename = filename
        self.backend = None
        self.handle = None
        self._open()

    def _open(self):
        try:
            import uproot
            self.handle = uproot.open(self.filename)
            self.backend = "uproot"
            return
        except ImportError:
            pass
        except Exception as exc:
            raise RuntimeError(f"could not open {self.filename}: {exc}") from exc

        try:
            import ROOT
        except ImportError as exc:
            raise RuntimeError(
                "neither uproot nor PyROOT is available; load the project's ROOT "
                "environment or install uproot"
            ) from exc
        ROOT.gROOT.SetBatch(True)
        self.handle = ROOT.TFile.Open(self.filename, "READ")
        if not self.handle or self.handle.IsZombie():
            raise RuntimeError(f"could not open {self.filename}")
        self.backend = "PyROOT"

    def branches(self, tree_name):
        if self.backend == "uproot":
            if tree_name not in self.handle:
                return None
            return set(self.handle[tree_name].keys())
        tree = self.handle.Get(tree_name)
        if not tree:
            return None
        return {branch.GetName() for branch in tree.GetListOfBranches()}

    def rows(self, tree_name, branches):
        if self.backend == "uproot":
            # A production Tracks tree can be tens of gigabytes.  Iterate in
            # bounded chunks rather than materializing the complete tree.
            for arrays in self.handle[tree_name].iterate(
                    list(branches), library="np", step_size="64 MB"):
                count = len(arrays[branches[0]])
                for index in range(count):
                    yield {name: _plain(arrays[name][index]) for name in branches}
            return
        tree = self.handle.Get(tree_name)
        for entry in tree:
            yield {name: _plain(getattr(entry, name)) for name in branches}

    def close(self):
        if self.backend == "PyROOT" and self.handle:
            self.handle.Close()


class Report:
    def __init__(self, max_examples):
        self.max_examples = max_examples
        self.checks = []
        self.metrics = {}

    def add(self, name, status, detail, examples=None):
        item = {"name": name, "status": status, "detail": detail}
        if examples:
            item["examples"] = list(examples)[:self.max_examples]
        self.checks.append(item)

    def outcome(self):
        if any(check["status"] == "FAIL" for check in self.checks):
            return "FAIL"
        if any(check["status"] == "WARN" for check in self.checks):
            return "WARN"
        return "PASS"


def count_kinematic_lines(path):
    valid = 0
    malformed = []
    with open(path, encoding="utf-8") as source:
        for number, line in enumerate(source, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            try:
                if len(fields) != 5:
                    raise ValueError
                tuple(float(field) for field in fields)
                valid += 1
            except ValueError:
                malformed.append(number)
    return valid, malformed


def trigger_multiplicity(row):
    s1x = row["S1XEnergy"] > 0.5
    s1y = row["S1YEnergy"] > 0.5 and abs(row["S1YTime"] - row["S1XTime"]) < 20.0
    s2x = row["S2XEnergy"] > 0.5 and abs(row["S2XTime"] - row["S1XTime"]) < 20.0
    s2y = row["S2YNPE"] > 100 and abs(row["S2YTime"] - row["S1XTime"]) < 20.0
    return int(s1x) + int(s1y) + int(s2x) + int(s2y)


def validate(args):
    report = Report(args.max_examples)
    reader = TreeReader(args.input)
    try:
        pa_found = reader.branches("PA")
        tracks_found = reader.branches("Tracks")
        missing_trees = [name for name, found in (("PA", pa_found), ("Tracks", tracks_found)) if found is None]
        if missing_trees:
            report.add("ROOT trees", "FAIL", "missing required tree(s): " + ", ".join(missing_trees))
            return report, reader.backend, None

        missing_pa = sorted(set(PA_BRANCHES) - pa_found)
        missing_tracks = sorted(set(TRACK_BRANCHES) - tracks_found)
        if missing_pa or missing_tracks:
            detail = f"missing PA branches={missing_pa}; missing Tracks branches={missing_tracks}"
            report.add("ROOT schema", "FAIL", detail)
            return report, reader.backend, None
        report.add("ROOT schema", "PASS", "PA and Tracks contain all required branches")

        pa_by_event = {}
        duplicate_pa = []
        trigger_counts = collections.Counter()
        trigger_by_event = {}
        trigger_pass = set()
        npe_nonzero = collections.Counter()
        invalid_npe = []
        for row in reader.rows("PA", PA_BRANCHES):
            event_id = int(row["EventID"])
            if event_id in pa_by_event:
                duplicate_pa.append(event_id)
            pa_by_event[event_id] = row
            multiplicity = trigger_multiplicity(row)
            trigger_counts[multiplicity] += 1
            trigger_by_event[event_id] = multiplicity
            if multiplicity >= 3:
                trigger_pass.add(event_id)
            for branch in NPE_BRANCHES:
                value = int(row[branch])
                if value < 0:
                    invalid_npe.append((event_id, branch, value))
                if value > 0:
                    npe_nonzero[branch] += 1

        track_events = set()
        primary_by_event = {}
        duplicate_track_count = 0
        duplicate_tracks = []
        unresolved_count = 0
        unresolved = []
        bad_primary_count_total = 0
        bad_primary_count = []
        noncontiguous_count = 0
        noncontiguous_events = []
        pion_reached = set()
        pion_exited = set()
        primary_reached = set()
        primary_exited = set()
        primary_reached_s2 = set()
        optical_count = 0
        optical_candidates = []
        track_exit_count = 0
        track_exit_without_reach = []
        track_entry_count = 0

        inspect_ids = set()
        if args.inspect_events:
            for multiplicity in range(5):
                candidates = [event for event, value in trigger_by_event.items()
                              if value == multiplicity]
                if candidates and len(inspect_ids) < args.inspect_events:
                    inspect_ids.add(min(candidates))
            for event_id in sorted(pa_by_event):
                if len(inspect_ids) >= args.inspect_events:
                    break
                inspect_ids.add(event_id)
        inspect_rows = collections.defaultdict(list) if inspect_ids else None
        wanted_pdg = PDG_CODES[args.particle]

        current_event = None
        current_ids = set()
        current_parents = []
        current_primaries = []

        def keep_example(items, value):
            if len(items) < args.max_examples:
                items.append(value)

        def finish_event():
            nonlocal unresolved_count, bad_primary_count_total
            if current_event is None:
                return
            for track_id, parent_id in current_parents:
                if parent_id not in current_ids:
                    unresolved_count += 1
                    keep_example(unresolved, (current_event, track_id, parent_id))
            expected = [item for item in current_primaries if item[1] == wanted_pdg]
            if len(current_primaries) != 1 or len(expected) != 1:
                bad_primary_count_total += 1
                keep_example(bad_primary_count,
                             (current_event, len(current_primaries), len(expected)))
            if len(current_primaries) == 1:
                primary_by_event[current_event] = current_primaries[0]

        for row in reader.rows("Tracks", TRACK_BRANCHES):
            event_id = int(row["EventID"])
            track_id = int(row["TrackID"])
            parent_id = int(row["ParentID"])
            pdg = int(row["PDG"])
            track_entry_count += 1

            if event_id != current_event:
                finish_event()
                if event_id in track_events:
                    noncontiguous_count += 1
                    keep_example(noncontiguous_events, event_id)
                track_events.add(event_id)
                current_event = event_id
                current_ids = set()
                current_parents = []
                current_primaries = []

            if track_id in current_ids:
                duplicate_track_count += 1
                keep_example(duplicate_tracks, (event_id, track_id))
            current_ids.add(track_id)
            if parent_id:
                current_parents.append((track_id, parent_id))
            if parent_id == 0:
                current_primaries.append(
                    (track_id, pdg, int(row["ReachedCal"]), int(row["ExitedCal"])))
            if pdg == 0:
                optical_count += 1
                keep_example(optical_candidates,
                             (event_id, track_id, str(row["CreatorProcess"])))
            if int(row["ExitedCal"]) and not int(row["ReachedCal"]):
                track_exit_count += 1
                keep_example(track_exit_without_reach, (event_id, track_id, pdg))
            if pdg == wanted_pdg:
                if int(row["ReachedCal"]):
                    pion_reached.add(event_id)
                if int(row["ExitedCal"]):
                    pion_exited.add(event_id)
                if parent_id == 0:
                    if int(row["ReachedCal"]):
                        primary_reached.add(event_id)
                    if int(row["ExitedCal"]):
                        primary_exited.add(event_id)
                    if int(row["ReachedS2"]):
                        primary_reached_s2.add(event_id)
            if inspect_rows is not None and event_id in inspect_ids:
                inspect_rows[event_id].append(row)

        finish_event()

        all_events = set(pa_by_event) | track_events
        generated = len(all_events)
        report.metrics.update({
            "pa_entries": len(pa_by_event) + len(duplicate_pa),
            "track_entries": track_entry_count,
            "generated_events": generated,
            "trigger_multiplicity": {str(key): trigger_counts[key] for key in range(5)},
            "trigger_pass": len(trigger_pass),
            "primary_reached_cal": len(primary_reached),
            "primary_exited_cal": len(primary_exited),
            "any_pion_reached_cal": len(pion_reached),
            "any_pion_exited_cal": len(pion_exited),
        })

        if duplicate_pa:
            report.add("PA EventID uniqueness", "FAIL", f"{len(duplicate_pa)} duplicate PA rows", duplicate_pa)
        else:
            report.add("PA EventID uniqueness", "PASS", f"{len(pa_by_event)} unique event rows")
        if duplicate_track_count or noncontiguous_count:
            report.add("track-key uniqueness", "FAIL",
                       f"{duplicate_track_count} duplicate keys; {noncontiguous_count} non-contiguous event blocks",
                       duplicate_tracks + noncontiguous_events)
        else:
            report.add("track-key uniqueness", "PASS", "all (EventID, TrackID) keys are unique")

        report.add("parent genealogy", "FAIL" if unresolved_count else "PASS",
                   f"{unresolved_count} nonzero ParentID values do not resolve" if unresolved_count else
                   "every nonzero ParentID resolves within its event", unresolved)

        bad_primary_pdg = []
        bad_primary_pdg_count = 0
        for event_id, primary in primary_by_event.items():
            track_id, pdg, _, _ = primary
            if pdg != wanted_pdg:
                bad_primary_pdg_count += 1
                keep_example(bad_primary_pdg, (event_id, track_id, pdg))
        report.add("expected primary", "FAIL" if bad_primary_count_total else "PASS",
                   f"{bad_primary_count_total} events lack exactly one {args.particle} primary" if bad_primary_count_total else
                   f"all {generated} events have exactly one PDG {wanted_pdg} primary", bad_primary_count)
        report.add("primary PDG", "FAIL" if bad_primary_pdg_count else "PASS",
                   f"{bad_primary_pdg_count} primary rows have unexpected PDG" if bad_primary_pdg_count else
                   f"all primary track rows have PDG {wanted_pdg}", bad_primary_pdg)

        for label, reached, exited in (
            ("primary calorimeter ordering", primary_reached, primary_exited),
            ("any-pion calorimeter ordering", pion_reached, pion_exited),
        ):
            bad = sorted(exited - reached)
            counts_ok = len(exited) <= len(reached) <= generated
            report.add(label, "FAIL" if bad or not counts_ok else "PASS",
                       f"Exited={len(exited)}, Reached={len(reached)}, Generated={generated}", bad)

        summary_mismatch = []
        summary_mismatch_count = 0
        zero_primary_pdg = []
        zero_primary_pdg_count = 0
        missing_pa_events = sorted(track_events - set(pa_by_event))
        missing_track_events = sorted(set(pa_by_event) - track_events)
        for event_id, pa_row in pa_by_event.items():
            if int(pa_row["PrimaryExitedCal"]) and not int(pa_row["PrimaryReachedCal"]):
                summary_mismatch_count += 1
                keep_example(summary_mismatch, (event_id, "PA ExitedCal without ReachedCal"))
            if int(pa_row["PrimaryPDG"]) == 0:
                zero_primary_pdg_count += 1
                keep_example(zero_primary_pdg, event_id)
            track = primary_by_event.get(event_id)
            if track is not None and track[1] == wanted_pdg:
                comparisons = (
                    ("PDG", int(pa_row["PrimaryPDG"]), track[1]),
                    ("ReachedCal", int(pa_row["PrimaryReachedCal"]), track[2]),
                    ("ExitedCal", int(pa_row["PrimaryExitedCal"]), track[3]),
                )
                for field, left, right in comparisons:
                    if left != right:
                        summary_mismatch_count += 1
                        keep_example(summary_mismatch, (event_id, field, left, right))
        report.add("PA/primary-track agreement", "FAIL" if summary_mismatch_count else "PASS",
                   f"{summary_mismatch_count} summary discrepancies" if summary_mismatch_count else
                   "event summaries agree with primary track rows", summary_mismatch)
        report.add("PA/Tracks event join", "FAIL" if missing_pa_events or missing_track_events else "PASS",
                   f"missing PA={len(missing_pa_events)}, missing Tracks={len(missing_track_events)}",
                   [("missing PA", event) for event in missing_pa_events] +
                   [("missing Tracks", event) for event in missing_track_events])
        report.add("PrimaryPDG zero rows", "FAIL" if zero_primary_pdg_count else "PASS",
                   f"{zero_primary_pdg_count} PA rows have PrimaryPDG=0" if zero_primary_pdg_count else
                   "no PA rows have PrimaryPDG=0", zero_primary_pdg)

        report.add("optical-photon filtering", "FAIL" if optical_count else "PASS",
                   f"{optical_count} PDG=0 track rows may be optical photons" if optical_count else
                   "Tracks contains no PDG=0 optical-photon candidates", optical_candidates)
        report.add("per-track calorimeter flags", "FAIL" if track_exit_count else "PASS",
                   f"{track_exit_count} tracks exited without reaching Cal" if track_exit_count else
                   "every ExitedCal track also has ReachedCal=1", track_exit_without_reach)
        if invalid_npe:
            report.add("NPE summaries", "FAIL", f"{len(invalid_npe)} negative NPE values", invalid_npe)
        elif not sum(npe_nonzero.values()):
            report.add("NPE summaries", "WARN", "all four NPE branches are zero; verify the sample should produce photons")
        else:
            populated = ", ".join(f"{key}={value}" for key, value in npe_nonzero.items())
            report.add("NPE summaries", "PASS", "events with nonzero NPE: " + populated)

        secondary_trigger = sorted(trigger_pass - primary_reached_s2)
        report.metrics["trigger_pass_without_primary_reached_s2"] = len(secondary_trigger)
        report.add("trigger versus primary ReachedS2", "WARN" if secondary_trigger else "PASS",
                   f"{len(secondary_trigger)} trigger-passing events have primary ReachedS2=0" if secondary_trigger else
                   "every trigger-passing event has primary ReachedS2=1", secondary_trigger)

        if args.kinematics:
            line_count, malformed = count_kinematic_lines(args.kinematics)
            report.metrics["kinematic_lines"] = line_count
            if malformed:
                report.add("kinematics format", "FAIL", f"{len(malformed)} malformed data lines", malformed)
            else:
                report.add("kinematics format", "PASS", f"{line_count} valid five-column data lines")
            requested = args.beam_on if args.beam_on is not None else len(pa_by_event)
            status = "PASS" if len(pa_by_event) == requested else "FAIL"
            report.add("input/event accounting", status,
                       f"PA={len(pa_by_event)}, requested={requested}, input lines={line_count}")
            if requested > line_count:
                report.add("input exhaustion", "FAIL",
                           "requested events exceed input rows; the current generator reuses the final row")
            elif requested == line_count:
                report.add("input exhaustion", "PASS", "run consumed exactly all available input rows")
            else:
                report.add("input exhaustion", "PASS",
                           f"run stopped {line_count - requested} rows before input exhaustion")

        inspection = (inspect_rows, trigger_by_event) if inspect_rows is not None else None
        return report, reader.backend, inspection
    finally:
        reader.close()


def choose_inspection_events(rows_by_event, trigger_by_event, limit):
    chosen = []
    # First include one event from each available trigger-multiplicity class.
    for multiplicity in range(5):
        candidates = [event for event, value in trigger_by_event.items()
                      if value == multiplicity and event in rows_by_event]
        if candidates and len(chosen) < limit:
            chosen.append(min(candidates))
    # Then prefer complex genealogies and fill deterministically by EventID.
    ordered = sorted(rows_by_event, key=lambda event: (-len(rows_by_event[event]), event))
    for event_id in ordered:
        if len(chosen) >= limit:
            break
        if event_id not in chosen:
            chosen.append(event_id)
    return chosen


def print_genealogies(inspection, limit):
    if not inspection or limit <= 0:
        return
    rows_by_event, trigger_by_event = inspection
    print("\nRepresentative event genealogies")
    for event_id in choose_inspection_events(rows_by_event, trigger_by_event, limit):
        rows = sorted(rows_by_event[event_id], key=lambda row: int(row["TrackID"]))
        print(f"  Event {event_id} (trigger multiplicity {trigger_by_event.get(event_id, 'n/a')}, "
              f"{len(rows)} tracks)")
        for row in rows:
            print("    T{track} <- P{parent}: PDG={pdg}, creator={creator}, end={end}, "
                  "S2={s2}, Cal={cal}, Exit={exit}".format(
                      track=int(row["TrackID"]), parent=int(row["ParentID"]), pdg=int(row["PDG"]),
                      creator=row["CreatorProcess"], end=row["EndProcess"],
                      s2=int(row["ReachedS2"]), cal=int(row["ReachedCal"]),
                      exit=int(row["ExitedCal"])))


def main():
    args = parse_args()
    if args.max_examples < 0 or args.inspect_events < 0:
        sys.exit("Error: --max-examples and --inspect-events must be nonnegative")
    if not os.path.isfile(args.input):
        sys.exit(f"Error: input file does not exist: {args.input}")
    try:
        report, backend, inspection_rows = validate(args)
    except (OSError, RuntimeError) as exc:
        sys.exit(f"Error: {exc}")

    print(f"Input:   {args.input}")
    print(f"Reader:  {backend}")
    for check in report.checks:
        print(f"[{check['status']:<4}] {check['name']}: {check['detail']}")
        if check.get("examples"):
            print("       examples: " + ", ".join(map(str, check["examples"])))
    print(f"Overall: {report.outcome()}")

    print_genealogies(inspection_rows, args.inspect_events)

    payload = {
        "input": os.path.abspath(args.input),
        "particle": args.particle,
        "reader": backend,
        "outcome": report.outcome(),
        "metrics": report.metrics,
        "checks": report.checks,
    }
    if args.json_path:
        with open(args.json_path, "w", encoding="utf-8") as output:
            json.dump(payload, output, indent=2, sort_keys=True)
            output.write("\n")
        print(f"JSON:    {args.json_path}")

    return 1 if report.outcome() == "FAIL" else 0


if __name__ == "__main__":
    sys.exit(main())
