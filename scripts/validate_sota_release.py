#!/usr/bin/env python3
"""Hard validation for data/spear_sota."""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "spear_sota"
AA = set("ACDEFGHIKLMNPQRSTVWY")
errors = []


def rows(name):
    return list(csv.DictReader((DATA / name).open(newline="", encoding="utf-8"), delimiter="\t"))


def fail(message):
    errors.append(message)


wetlab = rows("df_wetlab_validated_sota.tsv")
if list(wetlab[0]) != ["sequence", "wetlab_experiment", "mic", "paper_title", "paper_doi"]:
    fail("five-column mentor table schema changed")
if len(wetlab) != 2449 or len({r["sequence"] for r in wetlab}) != 596:
    fail("unexpected wet-lab counts")
if len({r["paper_doi"] for r in wetlab}) != 22:
    fail("not all 22 primary papers are present")
for i, row in enumerate(wetlab, 2):
    if not row["sequence"] or not set(row["sequence"]) <= AA:
        fail(f"wet-lab row {i}: invalid sequence")
    if not row["wetlab_experiment"] or not row["paper_title"] or not row["paper_doi"]:
        fail(f"wet-lab row {i}: missing required evidence field")

structured = rows("df_wetlab_assays_structured.tsv")
if len(structured) != len(wetlab):
    fail("structured assay table is not row-complete")
if len({r["assay_id"] for r in structured}) != len(structured):
    fail("assay_id is not unique")

decisions = rows("df_wetlab_compound_decisions.tsv")
if len(decisions) != 601:
    fail("unexpected compound-paper decision count")
if sum(r["segmentation_decision"] == "include" for r in decisions) != 270:
    fail("unexpected strict inclusion count")

gold = rows("df_spear_segmentation_gold.tsv")
if len(gold) != 270 or len({r["gold_label_id"] for r in gold}) != len(gold):
    fail("unexpected or duplicate strict-gold labels")
for i, row in enumerate(gold, 2):
    a, b = int(row["peptide_start_0based"]), int(row["peptide_end_0based_exclusive"])
    if row["source_peptide"][a:b] != row["encrypted_peptide"]:
        fail(f"gold row {i}: coordinate mismatch")
    if len(row["source_peptide"]) <= len(row["encrypted_peptide"]):
        fail(f"gold row {i}: parent is not longer")
    if row["recommended_split"] != "wetlab_holdout":
        fail(f"gold row {i}: not assigned to holdout")

split_rows = rows("df_spear_homology_splits.tsv")
cluster_splits = {}
source_splits = {}
for row in split_rows:
    previous = cluster_splits.setdefault(row["homology_cluster"], row["split"])
    if previous != row["split"]:
        fail(f"homology cluster {row['homology_cluster']} leaks across splits")
    if row["source_id"] in source_splits:
        fail(f"source {row['source_id']} appears more than once")
    source_splits[row["source_id"]] = row["split"]

json_sources = set()
for split in ("train", "val", "test", "wetlab_holdout"):
    for line_number, line in enumerate((DATA / "segmentation_splits" / f"{split}.jsonl").open(), 1):
        item = json.loads(line)
        sid = item["source_id"]
        if sid in json_sources:
            fail(f"source {sid} appears in multiple JSONL records")
        json_sources.add(sid)
        if item["split"] != split or source_splits.get(sid) != split:
            fail(f"{split}:{line_number}: split metadata disagreement")
        if len(item["cleavage_sites"]) != len(item["encrypted_peptides"]):
            fail(f"{split}:{line_number}: site/peptide count mismatch")
        for (a, b), peptide in zip(item["cleavage_sites"], item["encrypted_peptides"]):
            if item["source_sequence"][a:b] != peptide:
                fail(f"{split}:{line_number}: coordinate mismatch")
if json_sources != set(source_splits):
    fail("JSONL source set differs from split manifest")

challenge = rows("df_amp_challenge_similarity.tsv")
if len(challenge) != 596:
    fail("challenge audit does not cover all sequences")
if sum(r["exact_reference_match"] == "True" for r in challenge) != 23:
    fail("unexpected challenge exact-match count")
if sum(float(r["max_levenshtein_ratio"]) > 0.8 for r in challenge) != 93:
    fail("unexpected challenge >0.8 similarity count")

broad_view = rows("SPEAR_596_wetlab_validated.tsv")
if len(broad_view) != 596 or len({r["sequence"] for r in broad_view}) != 596:
    fail("broad handoff view must contain exactly 596 unique sequences")
if any(r["wetlab_validated"] != "True" for r in broad_view):
    fail("broad handoff view contains a row without positive wet-lab evidence")

strict_view = rows("SPEAR_270_strict_segmentation.tsv")
if strict_view != gold:
    fail("strict handoff view differs from the validated segmentation gold table")

activity_peptides = rows("SPEAR_596_activity_peptides.tsv")
if len(activity_peptides) != 596 or len({r["sequence"] for r in activity_peptides}) != 596:
    fail("activity peptide index must contain exactly 596 unique sequences")
if {r["activity_split"] for r in activity_peptides} != {"train", "val", "test"}:
    fail("activity peptide index must contain train, val, and test splits")
activity_cluster_splits = {}
for row in activity_peptides:
    cluster = row["activity_homology_cluster"]
    split = row["activity_split"]
    if cluster in activity_cluster_splits and activity_cluster_splits[cluster] != split:
        fail(f"activity cluster crosses splits: {cluster}")
    activity_cluster_splits[cluster] = split

activity_measurements = rows("SPEAR_2449_activity_measurements.tsv")
if len(activity_measurements) != 2449:
    fail("activity measurement table must contain exactly 2,449 assay rows")
peptide_lookup = {r["peptide_id"]: r for r in activity_peptides}
for row in activity_measurements:
    peptide = peptide_lookup.get(row["peptide_id"])
    if not peptide or peptide["sequence"] != row["sequence"]:
        fail(f"activity measurement has invalid peptide link: {row['assay_id']}")
        continue
    if peptide["activity_homology_cluster"] != row["activity_homology_cluster"] or peptide["activity_split"] != row["activity_split"]:
        fail(f"activity measurement crosses its peptide split: {row['assay_id']}")
if sum(r["regression_eligible"] == "True" for r in activity_measurements) != 2352:
    fail("unexpected activity MIC regression row count")
if sum(r["primary_uM_regression_eligible"] == "True" for r in activity_measurements) != 1957:
    fail("unexpected primary µM regression row count")
if sum(r["regression_task"] == "MIC_ug_per_mL" for r in activity_measurements) != 395:
    fail("unexpected secondary mass-concentration regression row count")

if errors:
    print("FAIL")
    for error in errors[:100]:
        print("-", error)
    sys.exit(1)
print("PASS")
print(f"{len(wetlab)} assays; 596 sequences; 270 strict-gold labels")
print(f"{len(split_rows)} unique parent sequences; {len(cluster_splits)} homology clusters; no leakage")
