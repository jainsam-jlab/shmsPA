#!/usr/bin/env python3
"""Event-by-event comparison of old CalcMissed and track-truth selections.

This is a read-only diagnostic.  It reproduces the trigger and sequential PID
logic from CalcMissed.py exactly, joins PA and Tracks by EventID, and compares
those event-response selections with primary-pi+ and any-pi+ ReachedS2 truth.
"""

import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[2]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from validate_output import TreeReader  # noqa: E402


PION_PDG = 211
ENERGY_THRESHOLD_MEV = 0.5
TIME_WINDOW_NS = 20.0
S2Y_NPE_THRESHOLD = 100

PA_BRANCHES = (
    "EventID", "PrimaryPDG",
    "S1XTime", "S1YTime", "S2XTime", "S2YTime",
    "S1XEnergy", "S1YEnergy", "S2XEnergy", "S2YEnergy",
    "S2YNPE", "AGCNPE", "HGCNPE", "NGCNPE", "CalEnergy",
)
TRACK_BRANCHES = (
    "EventID", "TrackID", "ParentID", "PDG",
    "ReachedS1", "ReachedHGC", "ReachedAGC", "ReachedS2", "ReachedCal",
    "EndVolume", "EndProcess",
)
PID_CATEGORIES = ("pion-like", "proton-like", "kaon-like",
                  "positron-like", "contamination", "trigger-fail")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="ROOT file containing PA and Tracks")
    parser.add_argument("central_momentum", type=float,
                        help="central momentum in GeV, used by old e+ PID")
    parser.add_argument("--output-dir", default="v1_v3_diagnostics",
                        help="diagnostic output directory")
    parser.add_argument("--no-plots", action="store_true",
                        help="skip optional matplotlib diagnostic plots")
    return parser.parse_args()


def require_schema(reader, tree_name, required):
    found = reader.branches(tree_name)
    if found is None:
        raise RuntimeError(f"input does not contain the {tree_name} tree")
    missing = sorted(set(required) - found)
    if missing:
        raise RuntimeError(f"{tree_name} missing branches: {', '.join(missing)}")


def trigger_multiplicity(row):
    """Exact 3/4 trigger logic from CalcMissed.py."""
    s1x = row["S1XEnergy"] > ENERGY_THRESHOLD_MEV
    s1y = (row["S1YEnergy"] > ENERGY_THRESHOLD_MEV and
           abs(row["S1YTime"] - row["S1XTime"]) < TIME_WINDOW_NS)
    s2x = (row["S2XEnergy"] > ENERGY_THRESHOLD_MEV and
           abs(row["S2XTime"] - row["S1XTime"]) < TIME_WINDOW_NS)
    s2y = (row["S2YNPE"] > S2Y_NPE_THRESHOLD and
           abs(row["S2YTime"] - row["S1XTime"]) < TIME_WINDOW_NS)
    return int(s1x) + int(s1y) + int(s2x) + int(s2y)


def old_pid_classification(row, central_momentum_gev):
    """Reproduce the sequential v1 classification, including its ordering."""
    if trigger_multiplicity(row) < 3:
        return "trigger-fail"

    agc = int(row["AGCNPE"])
    hgc = int(row["HGCNPE"])
    ngc = int(row["NGCNPE"])
    momentum_mev = central_momentum_gev * 1000.0
    cal_ratio = row["CalEnergy"] / momentum_mev

    if agc < 20 and hgc < 50:
        return "proton-like"
    if agc > 20 and hgc < 50:
        return "kaon-like"
    if agc > 20 and hgc > 50 and ngc > 5 and cal_ratio > 0.7:
        return "positron-like"
    if agc > 20 and hgc > 50:
        return "pion-like"
    return "contamination"


def percent(count, denominator):
    return 100.0 * count / denominator if denominator else 0.0


def particle_group(pdg):
    """Group a reached-S2 secondary by PDG; this does not assign causality."""
    if pdg == 211:
        return "secondary_pi+"
    if pdg == 2212:
        return "secondary_proton"
    if abs(pdg) == 11:
        return "secondary_electron_or_positron"
    if abs(pdg) == 13:
        return "secondary_muon"
    if abs(pdg) in (211, 321) or pdg > 1000000000:
        return "other_charged_secondary"
    return "other_or_neutral_secondary"


def compact_track(row):
    return {
        "TrackID": int(row["TrackID"]),
        "ParentID": int(row["ParentID"]),
        "PDG": int(row["PDG"]),
        "ReachedS1": int(row["ReachedS1"]),
        "ReachedHGC": int(row["ReachedHGC"]),
        "ReachedAGC": int(row["ReachedAGC"]),
        "ReachedS2": int(row["ReachedS2"]),
        "ReachedCal": int(row["ReachedCal"]),
        "EndVolume": str(row["EndVolume"]),
        "EndProcess": str(row["EndProcess"]),
    }


def read_input(input_path, central_momentum):
    reader = TreeReader(input_path)
    try:
        require_schema(reader, "PA", PA_BRANCHES)
        require_schema(reader, "Tracks", TRACK_BRANCHES)

        pa_by_event = {}
        duplicate_pa = []
        for row in reader.rows("PA", PA_BRANCHES):
            event_id = int(row["EventID"])
            if event_id in pa_by_event:
                duplicate_pa.append(event_id)
            plain = {name: row[name] for name in PA_BRANCHES}
            plain["EventID"] = event_id
            plain["nTrig"] = trigger_multiplicity(plain)
            plain["oldPIDClassification"] = old_pid_classification(
                plain, central_momentum)
            pa_by_event[event_id] = plain

        primary_by_event = {}
        secondary_s2_by_event = {}
        any_pi_s2 = set()
        track_event_ids = set()
        noncontiguous_events = []
        malformed = []

        current_event = None
        current_parent_zero = []
        current_primary_pi = []
        current_secondary_s2 = []

        def finish_event():
            if current_event is None:
                return
            if len(current_parent_zero) != 1 or len(current_primary_pi) != 1:
                malformed.append({
                    "EventID": current_event,
                    "ParentID0_tracks": len(current_parent_zero),
                    "ParentID0_PDG211_tracks": len(current_primary_pi),
                })
            if len(current_primary_pi) == 1:
                primary_by_event[current_event] = current_primary_pi[0]
            # Secondary details are needed only for the reverse-case sample.
            # Keeping them for every event would defeat bounded-memory
            # streaming on the 581-million-row production Tracks tree.
            pa_row = pa_by_event.get(current_event)
            primary_missed_s2 = (len(current_primary_pi) == 1 and
                                 not current_primary_pi[0]["ReachedS2"])
            if (current_secondary_s2 and pa_row is not None and
                    pa_row["nTrig"] >= 3 and primary_missed_s2):
                secondary_s2_by_event[current_event] = list(current_secondary_s2)

        for row in reader.rows("Tracks", TRACK_BRANCHES):
            event_id = int(row["EventID"])
            if event_id != current_event:
                finish_event()
                if event_id in track_event_ids:
                    noncontiguous_events.append(event_id)
                track_event_ids.add(event_id)
                current_event = event_id
                current_parent_zero = []
                current_primary_pi = []
                current_secondary_s2 = []

            pdg = int(row["PDG"])
            parent_id = int(row["ParentID"])
            reached_s2 = int(row["ReachedS2"])
            if pdg == PION_PDG and reached_s2:
                any_pi_s2.add(event_id)
            if parent_id == 0:
                track = compact_track(row)
                current_parent_zero.append(track)
                if pdg == PION_PDG:
                    current_primary_pi.append(track)
            elif reached_s2:
                current_secondary_s2.append({
                    "TrackID": int(row["TrackID"]),
                    "ParentID": parent_id,
                    "PDG": pdg,
                    "group": particle_group(pdg),
                })
        finish_event()

        return {
            "reader": reader.backend,
            "pa_by_event": pa_by_event,
            "duplicate_pa": duplicate_pa,
            "primary_by_event": primary_by_event,
            "secondary_s2_by_event": secondary_s2_by_event,
            "any_pi_s2": any_pi_s2,
            "track_event_ids": track_event_ids,
            "noncontiguous_events": noncontiguous_events,
            "malformed": malformed,
        }
    finally:
        reader.close()


def classify_secondary_recovery(tracks):
    groups = collections.Counter(track["group"] for track in tracks)
    pdgs = collections.Counter(track["PDG"] for track in tracks)
    if not tracks:
        exclusive = "no_recorded_secondary_reached_S2"
    elif len(tracks) > 1 or len(groups) > 1:
        exclusive = "multiple_secondaries_reached_S2"
    else:
        exclusive = next(iter(groups))
    return exclusive, groups, pdgs


def build_results(data, input_path, central_momentum):
    pa = data["pa_by_event"]
    primary = data["primary_by_event"]
    pa_ids = set(pa)
    primary_ids = set(primary)
    generated_ids = pa_ids & primary_ids
    primary_s2 = {event for event in generated_ids if primary[event]["ReachedS2"]}
    any_pi_s2 = data["any_pi_s2"] & generated_ids
    trigger_pass = {event for event in generated_ids if pa[event]["nTrig"] >= 3}
    pion_pid = {event for event in generated_ids
                if pa[event]["oldPIDClassification"] == "pion-like"}

    global_pid = collections.Counter(pa[event]["oldPIDClassification"]
                                     for event in generated_ids)
    primary_s2_pid = collections.Counter(pa[event]["oldPIDClassification"]
                                         for event in primary_s2)
    any_s2_pid = collections.Counter(pa[event]["oldPIDClassification"]
                                     for event in any_pi_s2)

    primary_table = {
        "trigger_pass_reached_yes": len(trigger_pass & primary_s2),
        "trigger_pass_reached_no": len(trigger_pass - primary_s2),
        "trigger_fail_reached_yes": len(primary_s2 - trigger_pass),
        "trigger_fail_reached_no": len(generated_ids - trigger_pass - primary_s2),
    }
    any_table = {
        "trigger_pass_reached_yes": len(trigger_pass & any_pi_s2),
        "trigger_pass_reached_no": len(trigger_pass - any_pi_s2),
        "trigger_fail_reached_yes": len(any_pi_s2 - trigger_pass),
        "trigger_fail_reached_no": len(generated_ids - trigger_pass - any_pi_s2),
    }

    recovery_events = sorted(trigger_pass - primary_s2)
    recovery_exclusive = collections.Counter()
    recovery_group_presence = collections.Counter()
    recovery_pdg_presence = collections.Counter()
    recovery_details = {}
    for event in recovery_events:
        tracks = data["secondary_s2_by_event"].get(event, [])
        exclusive, groups, pdgs = classify_secondary_recovery(tracks)
        recovery_exclusive[exclusive] += 1
        for group in groups:
            recovery_group_presence[group] += 1
        for pdg in pdgs:
            recovery_pdg_presence[pdg] += 1
        recovery_details[event] = {
            "exclusive_category": exclusive,
            "tracks": tracks,
        }

    primary_overlap = {
        "reachedS2_and_v1Pion": len(primary_s2 & pion_pid),
        "reachedS2_not_v1Pion": len(primary_s2 - pion_pid),
        "v1Pion_not_reachedS2": len(pion_pid - primary_s2),
        "neither": len(generated_ids - primary_s2 - pion_pid),
    }
    any_overlap = {
        "reachedS2_and_v1Pion": len(any_pi_s2 & pion_pid),
        "reachedS2_not_v1Pion": len(any_pi_s2 - pion_pid),
        "v1Pion_not_reachedS2": len(pion_pid - any_pi_s2),
        "neither": len(generated_ids - any_pi_s2 - pion_pid),
    }

    generated_identity = sum(global_pid.values())
    primary_identity = sum(primary_s2_pid.values())
    any_identity = sum(any_s2_pid.values())
    assertions = {
        "one_PA_row_per_event": len(data["duplicate_pa"]) == 0,
        "PA_and_Tracks_event_ranges_align": pa_ids == data["track_event_ids"],
        "exactly_one_primary_pi+_per_event": len(data["malformed"]) == 0,
        "no_noncontiguous_track_event_blocks": len(data["noncontiguous_events"]) == 0,
        "generated_accounting_closes": generated_identity == len(generated_ids),
        "primary_S2_accounting_closes": primary_identity == len(primary_s2),
        "any_pi_S2_accounting_closes": any_identity == len(any_pi_s2),
        "primary_trigger_table_closes": sum(primary_table.values()) == len(generated_ids),
        "any_trigger_table_closes": sum(any_table.values()) == len(generated_ids),
    }

    result = {
        "input": os.path.abspath(input_path),
        "reader": data["reader"],
        "central_momentum_GeV": central_momentum,
        "definitions": {
            "primary_pi+": "ParentID == 0 and PDG == 211",
            "primary_ReachedS2": "primary track pre/post logical volume is S2XLogical or S2YLogical",
            "any_pi+_ReachedS2": "any PDG == 211 track has ReachedS2 == 1",
            "trigger": "at least 3 of the four exact CalcMissed.py hodoscope conditions",
            "v1_pion": "trigger pass, not earlier e+ classification, AGCNPE > 20 and HGCNPE > 50",
        },
        "validation": {
            "PA_rows": len(pa),
            "PA_unique_EventIDs": len(pa_ids),
            "Tracks_unique_EventIDs": len(data["track_event_ids"]),
            "EventID_min_PA": min(pa_ids) if pa_ids else None,
            "EventID_max_PA": max(pa_ids) if pa_ids else None,
            "EventID_min_Tracks": min(data["track_event_ids"]) if data["track_event_ids"] else None,
            "EventID_max_Tracks": max(data["track_event_ids"]) if data["track_event_ids"] else None,
            "duplicate_PA_EventIDs": data["duplicate_pa"],
            "malformed_primary_events": data["malformed"],
            "noncontiguous_track_event_blocks": data["noncontiguous_events"],
            "assertions": assertions,
        },
        "counts": {
            "generated_primary_pi+_events": len(generated_ids),
            "trigger_pass": len(trigger_pass),
            "trigger_fail": len(generated_ids - trigger_pass),
            "v1_pion_like": len(pion_pid),
            "primary_pi+_ReachedS2": len(primary_s2),
            "any_pi+_ReachedS2": len(any_pi_s2),
        },
        "table1_primary_S2_vs_trigger": primary_table,
        "table1b_any_pi_S2_vs_trigger": any_table,
        "table2_primary_reached_S2": {
            "N_reachedS2": len(primary_s2),
            "N_reachedS2_triggerPass": len(primary_s2 & trigger_pass),
            "N_reachedS2_triggerFail": len(primary_s2 - trigger_pass),
            "N_reachedS2_triggerPass_pionPIDPass": len(primary_s2 & pion_pid),
            "N_reachedS2_triggerPass_pionPIDFail": len(primary_s2 & trigger_pass - pion_pid),
        },
        "table3_secondary_recovery": {
            "events": len(recovery_events),
            "exclusive_categories": dict(recovery_exclusive),
            "group_presence_event_counts": dict(recovery_group_presence),
            "PDG_presence_event_counts": {str(k): v for k, v in recovery_pdg_presence.items()},
        },
        "table4_v1_classification_all_events": dict(global_pid),
        "table4_v1_classification_primary_reached_S2": dict(primary_s2_pid),
        "v1_classification_any_pi_reached_S2": dict(any_s2_pid),
        "primary_S2_vs_v1_pion_overlap": primary_overlap,
        "any_pi_S2_vs_v1_pion_overlap": any_overlap,
        "sets": {
            "generated": generated_ids,
            "primary_s2": primary_s2,
            "any_pi_s2": any_pi_s2,
            "trigger_pass": trigger_pass,
            "pion_pid": pion_pid,
            "recovery_events": set(recovery_events),
        },
        "recovery_details": recovery_details,
    }
    return result


EVENT_FIELDS = (
    "EventID", "PrimaryTrackID", "PrimaryPDG", "ReachedS1", "ReachedHGC",
    "ReachedAGC", "ReachedS2", "ReachedCal", "EndVolume", "EndProcess",
    "S1XEnergy", "S1YEnergy", "S2XEnergy", "S2YNPE",
    "S1XTime", "S1YTime", "S2XTime", "S2YTime",
    "AGCNPE", "HGCNPE", "NGCNPE", "CalEnergy", "nTrig",
    "oldPIDClassification", "secondaryTracksReachingS2",
)


def event_output_row(event_id, data, result):
    pa = data["pa_by_event"][event_id]
    primary = data["primary_by_event"][event_id]
    secondaries = result["recovery_details"].get(event_id, {}).get("tracks", [])
    secondary_text = ";".join(
        f"T{x['TrackID']}<-P{x['ParentID']}:PDG{x['PDG']}:{x['group']}"
        for x in secondaries)
    return {
        "EventID": event_id,
        "PrimaryTrackID": primary["TrackID"],
        "PrimaryPDG": primary["PDG"],
        "ReachedS1": primary["ReachedS1"],
        "ReachedHGC": primary["ReachedHGC"],
        "ReachedAGC": primary["ReachedAGC"],
        "ReachedS2": primary["ReachedS2"],
        "ReachedCal": primary["ReachedCal"],
        "EndVolume": primary["EndVolume"],
        "EndProcess": primary["EndProcess"],
        "S1XEnergy": pa["S1XEnergy"],
        "S1YEnergy": pa["S1YEnergy"],
        "S2XEnergy": pa["S2XEnergy"],
        "S2YNPE": pa["S2YNPE"],
        "S1XTime": pa["S1XTime"],
        "S1YTime": pa["S1YTime"],
        "S2XTime": pa["S2XTime"],
        "S2YTime": pa["S2YTime"],
        "AGCNPE": pa["AGCNPE"],
        "HGCNPE": pa["HGCNPE"],
        "NGCNPE": pa["NGCNPE"],
        "CalEnergy": pa["CalEnergy"],
        "nTrig": pa["nTrig"],
        "oldPIDClassification": pa["oldPIDClassification"],
        "secondaryTracksReachingS2": secondary_text,
    }


def write_event_csv(path, event_ids, data, result):
    with open(path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=EVENT_FIELDS)
        writer.writeheader()
        for event_id in sorted(event_ids):
            writer.writerow(event_output_row(event_id, data, result))


def make_plots(output_dir, data, result, central_momentum):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return ["matplotlib unavailable; plots skipped"]

    pa = data["pa_by_event"]
    sets = result["sets"]
    primary_s2 = sets["primary_s2"]
    trigger_pass = sets["trigger_pass"]
    pion_pid = sets["pion_pid"]
    warnings = []

    variables = (
        ("S1XEnergy", "S1X energy [MeV]"),
        ("S1YEnergy", "S1Y energy [MeV]"),
        ("S2XEnergy", "S2X energy [MeV]"),
        ("S2YNPE", "S2Y NPE"),
        ("dtS1Y", "tS1Y - tS1X [ns]"),
        ("dtS2X", "tS2X - tS1X [ns]"),
        ("dtS2Y", "tS2Y - tS1X [ns]"),
    )

    def values(events, name):
        if name == "dtS1Y":
            return [pa[e]["S1YTime"] - pa[e]["S1XTime"] for e in events]
        if name == "dtS2X":
            return [pa[e]["S2XTime"] - pa[e]["S1XTime"] for e in events]
        if name == "dtS2Y":
            return [pa[e]["S2YTime"] - pa[e]["S1XTime"] for e in events]
        return [pa[e][name] for e in events]

    pass_events = primary_s2 & trigger_pass
    fail_events = primary_s2 - trigger_pass
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    for ax, (name, label) in zip(axes.flat, variables):
        ax.hist(values(pass_events, name), bins=60, histtype="step", density=True,
                label=f"trigger pass ({len(pass_events)})")
        if fail_events:
            ax.hist(values(fail_events, name), bins=60, histtype="step", density=True,
                    label=f"trigger fail ({len(fail_events)})")
        else:
            ax.text(0.5, 0.8, "No primary-ReachedS2 trigger failures",
                    transform=ax.transAxes, ha="center")
        ax.set_xlabel(label)
        ax.set_ylabel("normalized events")
        ax.legend(fontsize=8)
    for ax in axes.flat[len(variables):]:
        ax.axis("off")
    fig.suptitle("Primary pi+ ReachedS2: trigger-pass versus trigger-fail")
    fig.tight_layout()
    fig.savefig(output_dir / "primary_reachedS2_trigger_comparison.png", dpi=160)
    plt.close(fig)

    pid_pass = primary_s2 & pion_pid
    pid_fail = primary_s2 & trigger_pass - pion_pid
    pid_vars = (
        ("AGCNPE", "AGC NPE"), ("HGCNPE", "HGC NPE"),
        ("NGCNPE", "NGC NPE"), ("CalRatio", "CalEnergy / P"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, (name, label) in zip(axes.flat, pid_vars):
        def pid_values(events):
            if name == "CalRatio":
                return [pa[e]["CalEnergy"] / (central_momentum * 1000.0) for e in events]
            return [pa[e][name] for e in events]
        ax.hist(pid_values(pid_pass), bins=70, histtype="step", density=True,
                label=f"v1 pion-like ({len(pid_pass)})")
        ax.hist(pid_values(pid_fail), bins=70, histtype="step", density=True,
                label=f"other v1 class ({len(pid_fail)})")
        ax.set_xlabel(label)
        ax.set_ylabel("normalized events")
        ax.legend(fontsize=8)
    fig.suptitle("Primary pi+ ReachedS2 and trigger pass: old PID pass/fail")
    fig.tight_layout()
    fig.savefig(output_dir / "primary_reachedS2_pid_comparison.png", dpi=160)
    plt.close(fig)

    groups = result["table3_secondary_recovery"]["group_presence_event_counts"]
    fig, ax = plt.subplots(figsize=(11, 6))
    labels, counts = zip(*sorted(groups.items(), key=lambda item: item[1], reverse=True))
    ax.bar(labels, counts)
    ax.set_ylabel("events (categories may overlap)")
    ax.set_title("Secondary species groups reaching S2 when primary missed S2 but trigger passed")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(output_dir / "secondary_recovery_species.png", dpi=160)
    plt.close(fig)
    return warnings


def sanitize_for_json(result):
    clean = dict(result)
    clean.pop("sets", None)
    clean.pop("recovery_details", None)
    return clean


def write_report(path, result):
    c = result["counts"]
    t1 = result["table1_primary_S2_vs_trigger"]
    t2 = result["table2_primary_reached_S2"]
    t3 = result["table3_secondary_recovery"]
    primary_overlap = result["primary_S2_vs_v1_pion_overlap"]
    any_overlap = result["any_pi_S2_vs_v1_pion_overlap"]
    global_pid = result["table4_v1_classification_all_events"]
    reached_pid = result["table4_v1_classification_primary_reached_S2"]
    n = c["generated_primary_pi+_events"]

    def pid_rows(counter):
        return "\n".join(
            f"| {category} | {counter.get(category, 0):,} | "
            f"{percent(counter.get(category, 0), sum(counter.values())):.3f}% |"
            for category in PID_CATEGORIES)

    text = f"""# V1–V3 Event-Selection Discrepancy

## Executive finding

The premise that approximately 97,760 **primary** pi+ tracks reached S2 is not
consistent with the current v3 implementation or saved output.  The exact
counts from this file are:

- primary pi+ ReachedS2: **{c['primary_pi+_ReachedS2']:,}**;
- any pi+ ReachedS2: **{c['any_pi+_ReachedS2']:,}**;
- old v1 pion-like selection: **{c['v1_pion_like']:,}**.

Thus the approximately 1,068-event numerical difference compares the v1
pion-like sample with the **any-pi+** S2 definition, not the primary-pi+
definition.  The samples are not nested, so their overlap—not their simple
difference—is the relevant accounting.

## Verified definitions

- Primary pi+: `ParentID == 0 && PDG == 211`; `TrackID == 1` is not assumed.
- Primary ReachedS2: that primary track has either its pre-step or post-step
  logical volume equal to `S2XLogical` or `S2YLogical`.  No energy deposition,
  timing cut, NPE cut, or requirement to traverse both planes is imposed.
- Any-pi+ ReachedS2: at least one `PDG == 211` track in the event has that flag.
- Old trigger: exact 0.5 MeV, 20 ns, and S2Y NPE > 100 cuts from
  `CalcMissed.py`, with at least three of four conditions required.
- Old pion-like classification: trigger pass, then sequential proton, kaon,
  positron, pion classification.  The positron test precedes the pion test,
  so an event satisfying AGC > 20 and HGC > 50 is classified positron-like
  when NGC > 5 and CalEnergy/P > 0.7; only the remainder is pion-like.

## EventID and primary validation

- PA rows: {result['validation']['PA_rows']:,}
- unique PA EventIDs: {result['validation']['PA_unique_EventIDs']:,}
- unique Tracks EventIDs: {result['validation']['Tracks_unique_EventIDs']:,}
- PA EventID range: {result['validation']['EventID_min_PA']}–{result['validation']['EventID_max_PA']}
- Tracks EventID range: {result['validation']['EventID_min_Tracks']}–{result['validation']['EventID_max_Tracks']}
- malformed primary events: {len(result['validation']['malformed_primary_events'])}
- all accounting assertions pass: {all(result['validation']['assertions'].values())}

## Table 1: primary S2 reach versus 3/4 trigger

| | Primary reached S2 | Primary missed S2 |
|---|---:|---:|
| Trigger pass | {t1['trigger_pass_reached_yes']:,} | {t1['trigger_pass_reached_no']:,} |
| Trigger fail | {t1['trigger_fail_reached_yes']:,} | {t1['trigger_fail_reached_no']:,} |

Primary pi+ tracks that reached S2 but failed the trigger:
**{t1['trigger_fail_reached_yes']:,}**.

## Table 2: old trigger and PID within primary-ReachedS2

- primary ReachedS2: {t2['N_reachedS2']:,};
- trigger pass: {t2['N_reachedS2_triggerPass']:,};
- trigger fail: {t2['N_reachedS2_triggerFail']:,};
- trigger pass and old pion-like PID: {t2['N_reachedS2_triggerPass_pionPIDPass']:,};
- trigger pass but old pion-like PID fail: {t2['N_reachedS2_triggerPass_pionPIDFail']:,}.

## Table 3: secondary recovery when primary missed S2

There are {t3['events']:,} events in which the primary missed S2 but the 3/4
trigger passed.  Presence of a secondary at S2 does not prove that it uniquely
caused the aggregate PA trigger response.

Exclusive categories:

```text
{json.dumps(t3['exclusive_categories'], indent=2, sort_keys=True)}
```

Species-group presence counts, which may overlap within an event:

```text
{json.dumps(t3['group_presence_event_counts'], indent=2, sort_keys=True)}
```

## Table 4: v1 classification for all generated events

| v1 category | Count | Percent of {n:,} |
|---|---:|---:|
{pid_rows(global_pid)}

Accounting closes exactly: {sum(global_pid.values()):,} = {n:,}.

## V1 classification restricted to primary-ReachedS2

| v1 category | Count | Percent of primary-ReachedS2 |
|---|---:|---:|
{pid_rows(reached_pid)}

Accounting closes exactly: {sum(reached_pid.values()):,} =
{c['primary_pi+_ReachedS2']:,}.

## Correct overlap accounting

### Primary pi+ ReachedS2 versus v1 pion-like

```text
intersection                       = {primary_overlap['reachedS2_and_v1Pion']:,}
primary ReachedS2, not v1 pion     = {primary_overlap['reachedS2_not_v1Pion']:,}
v1 pion, primary did not reach S2  = {primary_overlap['v1Pion_not_reachedS2']:,}
neither                            = {primary_overlap['neither']:,}
```

### Any pi+ ReachedS2 versus v1 pion-like

```text
intersection                   = {any_overlap['reachedS2_and_v1Pion']:,}
any pi+ ReachedS2, not v1 pion = {any_overlap['reachedS2_not_v1Pion']:,}
v1 pion, no pi+ reached S2      = {any_overlap['v1Pion_not_reachedS2']:,}
neither                        = {any_overlap['neither']:,}
```

The count difference obeys the set identity:

```text
N(any pi+ ReachedS2) - N(v1 pion)
  = N(any S2 but not v1 pion) - N(v1 pion but not any S2)
  = {any_overlap['reachedS2_not_v1Pion']:,} - {any_overlap['v1Pion_not_reachedS2']:,}
  = {c['any_pi+_ReachedS2'] - c['v1_pion_like']:,}.
```

## Physical interpretation

`1 - N_primaryReachedS2/N_generated` is a primary transport reach loss, not a
microscopic absorption probability: scattering, charge exchange, inelastic
production with descendants, world/acceptance loss, and other processes may
all remove the original track from the reach numerator.

`1 - N_v1Pion/N_generated` combines trigger response with sequential
detector-response PID classification.  It therefore includes trigger and PID
failures and cannot be called microscopic hadronic absorption.

## Diagnostic files

- `reachedS2_triggerFail.csv`
- `reachedS2_triggerPass_pidFail.csv`
- `missedS2_triggerPass.csv`
- `reachedS2_pionPIDPass.csv`
- `full_event_classification.csv`
- `primary_reachedS2_trigger_comparison.png`
- `primary_reachedS2_pid_comparison.png`
- `secondary_recovery_species.png`
"""
    path.write_text(text, encoding="utf-8")


def print_summary(result):
    c = result["counts"]
    t = result["table1_primary_S2_vs_trigger"]
    print("=" * 78)
    print(f"Generated primary pi+ events:   {c['generated_primary_pi+_events']:>10,}")
    print(f"Old v1 pion-like:               {c['v1_pion_like']:>10,}")
    print(f"Primary pi+ ReachedS2:          {c['primary_pi+_ReachedS2']:>10,}")
    print(f"Any pi+ ReachedS2:              {c['any_pi+_ReachedS2']:>10,}")
    print(f"3/4 trigger pass:               {c['trigger_pass']:>10,}")
    print(f"3/4 trigger fail:               {c['trigger_fail']:>10,}")
    print("-" * 78)
    print("Primary S2 versus trigger")
    print(f"  trigger pass, reached S2:     {t['trigger_pass_reached_yes']:>10,}")
    print(f"  trigger fail, reached S2:     {t['trigger_fail_reached_yes']:>10,}")
    print(f"  trigger pass, missed S2:      {t['trigger_pass_reached_no']:>10,}")
    print(f"  trigger fail, missed S2:      {t['trigger_fail_reached_no']:>10,}")
    print("-" * 78)
    for name, passed in result["validation"]["assertions"].items():
        print(f"[{'PASS' if passed else 'FAIL'}] {name}")
    print("=" * 78)


def main():
    args = parse_args()
    input_path = Path(args.input)
    if not input_path.is_file():
        sys.exit(f"Error: input file does not exist: {input_path}")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        data = read_input(str(input_path), args.central_momentum)
        result = build_results(data, str(input_path), args.central_momentum)
    except (OSError, RuntimeError) as exc:
        sys.exit(f"Error: {exc}")

    failed_assertions = [name for name, passed in
                         result["validation"]["assertions"].items() if not passed]
    if failed_assertions:
        raise RuntimeError("accounting/validation failures: " + ", ".join(failed_assertions))

    sets = result["sets"]
    write_event_csv(output_dir / "reachedS2_triggerFail.csv",
                    sets["primary_s2"] - sets["trigger_pass"], data, result)
    write_event_csv(output_dir / "reachedS2_triggerPass_pidFail.csv",
                    sets["primary_s2"] & sets["trigger_pass"] - sets["pion_pid"],
                    data, result)
    write_event_csv(output_dir / "missedS2_triggerPass.csv",
                    sets["trigger_pass"] - sets["primary_s2"], data, result)
    write_event_csv(output_dir / "reachedS2_pionPIDPass.csv",
                    sets["primary_s2"] & sets["pion_pid"], data, result)
    write_event_csv(output_dir / "full_event_classification.csv",
                    sets["generated"], data, result)

    plot_warnings = [] if args.no_plots else make_plots(
        output_dir, data, result, args.central_momentum)
    clean = sanitize_for_json(result)
    clean["plot_warnings"] = plot_warnings
    with open(output_dir / "selection_comparison_summary.json", "w", encoding="utf-8") as output:
        json.dump(clean, output, indent=2, sort_keys=True)
        output.write("\n")
    write_report(output_dir / "V1_V3_SELECTION_DISCREPANCY.md", result)
    print_summary(result)
    print(f"Results: {output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
