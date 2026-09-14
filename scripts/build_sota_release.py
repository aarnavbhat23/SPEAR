#!/usr/bin/env python3
"""Build the chemistry-aware SPEAR wet-lab and segmentation release.

The five-column table remains compatible with the schema requested by Vinay.
The structured and decision tables prevent modified/synthetic compounds from
silently becoming protein-cleavage labels.  MMseqs2 clusters parent proteins so
related parents cannot leak between model splits.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


AA = set("ACDEFGHIKLMNPQRSTVWY")
AMP_DOI = "10.1016/j.cell.2024.05.013"
AMPSPHERE_MIC = {
    "synechocucin-1": 8,
    "proteobacticin-1": 16,
    "actynomycin-1": 64,
    "lachnospirin-1": 2,
    "enterococcin-1": 1,
    "alphaprotecin-1": 1,
    "oscillospirin": 8,
    "ampspherin-4": 8,
    "methylocellin-1": 2,
    "reyranin-1": 16,
}


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], columns: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: row.get(k, "") for k in columns} for row in rows)


def stable_id(prefix: str, *parts: str, n: int = 16):
    digest = hashlib.sha256("\x1f".join(parts).encode()).hexdigest()[:n]
    return f"{prefix}_{digest}"


def bracket_locator(text: str):
    matches = re.findall(r"\[([^\]]+)\]", text)
    return "; ".join(matches)


def peptide_name(text: str):
    loc = bracket_locator(text)
    if not loc:
        return ""
    parts = [x.strip() for x in loc.split(";")]
    candidate = parts[-1]
    if re.search(r"figure|table|results|supplement", candidate, re.I):
        return ""
    return candidate


def assay_target(text: str):
    core = text.split(" [", 1)[0]
    match = re.search(r"\bagainst\s+(.+?)(?:;|$)", core, re.I)
    return match.group(1).strip() if match else "at least one assayed organism"


def chemical_form(text: str):
    lower = text.lower()
    is_d = "all-d" in lower or "d-amino" in lower or "retro-inverso" in lower
    is_amide = "c-terminal amide" in lower or "c-terminally amidated" in lower
    if is_d and is_amide:
        return "all_D; C_terminal_amide"
    if is_d:
        return "all_D_or_retro_inverso"
    if is_amide:
        return "C_terminal_amide"
    return "canonical_L_form_or_not_reported"


def parse_mic(value: str):
    value = value.strip()
    if not value:
        return "", "", ""
    match = re.fullmatch(r"(<=|>=|<|>)?\s*([0-9.]+)\s*(.+)", value)
    if not match:
        return "", value, ""
    relation = match.group(1) or "="
    return relation, match.group(2), match.group(3).replace("μ", "µ")


def fasta_sequences(path: Path):
    seqs, buf = [], []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            if buf:
                seqs.append("".join(buf).strip().upper())
                buf = []
        else:
            buf.append(raw.strip())
    if buf:
        seqs.append("".join(buf).strip().upper())
    return [s for s in seqs if s]


def challenge_similarity(sequences: list[str], fasta: Path):
    try:
        import Levenshtein
    except ImportError as exc:
        raise SystemExit("python-Levenshtein is required when --challenge-fasta is used") from exc
    references = fasta_sequences(fasta)
    refset = set(references)
    result = {}
    for i, seq in enumerate(sorted(sequences), 1):
        nearest, score = "", -1.0
        for ref in references:
            current = Levenshtein.ratio(seq, ref)
            if current > score:
                nearest, score = ref, current
                if score == 1.0:
                    break
        result[seq] = {
            "exact_reference_match": seq in refset,
            "nearest_reference_sequence": nearest,
            "max_levenshtein_ratio": round(score, 6),
            "passes_full_library_no_exact": seq not in refset,
            "passes_top100_similarity": score <= 0.8,
        }
        if i % 50 == 0:
            print(f"challenge similarity: {i}/{len(sequences)}")
    return result


def build_source_records(spec_rows, provenance_rows):
    records = {}
    row_mappings = defaultdict(list)
    for spec, prov in zip(spec_rows, provenance_rows):
        sequence = spec["source_peptide"]
        peptide = spec["encrypted_peptide"]
        if not sequence or peptide not in sequence:
            continue
        sid = stable_id("SRC", sequence)
        rec = records.setdefault(sid, {
            "source_id": sid,
            "source_sequence": sequence,
            "source_accessions": set(),
            "source_proteomes": set(),
            "dois": set(),
            "sites": set(),
            "has_existing_wetlab_label": False,
        })
        rec["source_accessions"].update(x for x in prov["parent_accessions"].split(";") if x)
        rec["source_proteomes"].add(spec["source_proteome"])
        rec["dois"].add(spec["doi_paper"])
        start = sequence.find(peptide)
        rec["sites"].add((start, start + len(peptide), peptide))
        if spec["wet_lab"] not in ("", "None", "not_reported"):
            rec["has_existing_wetlab_label"] = True
        row_mappings[(peptide, spec["doi_paper"])].append((spec, prov, sid, start))
    return records, row_mappings


def mmseqs_clusters(records: dict, executable: str):
    if not shutil.which(executable):
        raise SystemExit(f"MMseqs2 executable not found: {executable}")
    with tempfile.TemporaryDirectory(prefix="spear_sota_mmseqs_") as tmp:
        tmp = Path(tmp)
        fasta = tmp / "sources.fasta"
        with fasta.open("w", encoding="utf-8") as fh:
            for sid, rec in sorted(records.items()):
                fh.write(f">{sid}\n{rec['source_sequence']}\n")
        prefix = tmp / "clustered"
        work = tmp / "work"
        subprocess.run([
            executable, "easy-cluster", str(fasta), str(prefix), str(work),
            "--min-seq-id", "0.30", "-c", "0.80", "--cov-mode", "1",
            "--threads", "4",
        ], check=True)
        membership = {}
        with (tmp / "clustered_cluster.tsv").open(encoding="utf-8") as fh:
            for line in fh:
                representative, member = line.rstrip("\n").split("\t")
                membership[member] = stable_id("HC", representative, n=12)
        missing = set(records) - set(membership)
        if missing:
            raise SystemExit(f"MMseqs2 did not assign {len(missing)} sources")
        return membership


def assign_splits(records: dict, membership: dict, gold_source_ids: set[str], seed: int):
    clusters = defaultdict(list)
    for sid, cid in membership.items():
        clusters[cid].append(sid)
    forced = set()
    for cid, members in clusters.items():
        if any(records[s]["has_existing_wetlab_label"] or s in gold_source_ids for s in members):
            forced.add(cid)
    assign = {cid: "wetlab_holdout" for cid in forced}
    strata = defaultdict(list)
    for cid, members in clusters.items():
        if cid in forced:
            continue
        labels = sorted(p for s in members for p in records[s]["source_proteomes"])
        strata[labels[0] if labels else "unknown"].append(cid)
    rng = random.Random(seed)
    for label, cids in sorted(strata.items()):
        cids = sorted(cids)
        rng.shuffle(cids)
        n_val = round(len(cids) * 0.10)
        n_test = round(len(cids) * 0.10)
        for cid in cids[:n_val]:
            assign[cid] = "val"
        for cid in cids[n_val:n_val + n_test]:
            assign[cid] = "test"
        for cid in cids[n_val + n_test:]:
            assign[cid] = "train"
    return assign


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--challenge-fasta", type=Path, required=True)
    parser.add_argument("--challenge-commit", default="unknown")
    parser.add_argument("--mmseqs", default="mmseqs")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    data = args.repo / "data"
    out = data / "spear_sota"
    wetlab = read_tsv(data / "df_wetlab_validated.tsv")
    spec = read_tsv(data / "df_amp_encryptome.tsv")
    provenance = read_tsv(data / "df_amp_provenance.tsv")
    engineered = read_tsv(data / "df_amp_engineered.tsv")
    murine = read_tsv(data / "df_cesar_murine_validated.tsv")
    if len(spec) != len(provenance):
        raise SystemExit("encryptome and provenance row counts differ")

    # Source-read improvements that do not alter row identity.
    for row in wetlab:
        if row["paper_doi"] == AMP_DOI and not row["mic"]:
            lower = row["wetlab_experiment"].lower()
            matches = [name for name in AMPSPHERE_MIC if name in lower]
            if len(matches) != 1:
                raise SystemExit(f"cannot resolve AMPSphere name: {row}")
            row["mic"] = f"{AMPSPHERE_MIC[matches[0]]} µM"
            row["wetlab_experiment"] = row["wetlab_experiment"].replace(
                "in-vitro antimicrobial susceptibility;",
                "broth microdilution MIC against A. baumannii;",
            )
        if row["paper_doi"] == "10.1016/j.xcrp.2023.101459":
            row["wetlab_experiment"] += "; active against at least one tested strain (manual SI confirmation)"

    five_cols = ["sequence", "wetlab_experiment", "mic", "paper_title", "paper_doi"]
    wetlab.sort(key=lambda r: tuple(r[x] for x in ("paper_doi", "sequence", "wetlab_experiment", "mic")))
    write_tsv(out / "df_wetlab_validated_sota.tsv", wetlab, five_cols)

    similarity = challenge_similarity(sorted({r["sequence"] for r in wetlab}), args.challenge_fasta)
    challenge_rows = [{"sequence": seq, **similarity[seq]} for seq in sorted(similarity)]
    write_tsv(out / "df_amp_challenge_similarity.tsv", challenge_rows, [
        "sequence", "exact_reference_match", "nearest_reference_sequence",
        "max_levenshtein_ratio", "passes_full_library_no_exact", "passes_top100_similarity",
    ])

    murine_keys = {(r["sequence"], r["paper_doi"]) for r in murine if r["sequence"]}
    engineered_by_native = defaultdict(list)
    for row in engineered:
        if row["native_fragment"]:
            engineered_by_native[(row["native_fragment"], row["doi_paper"])].append(row)

    structured = []
    by_key = defaultdict(list)
    for row in wetlab:
        key = (row["sequence"], row["paper_doi"])
        relation, value, unit = parse_mic(row["mic"])
        form = chemical_form(row["wetlab_experiment"])
        eng = engineered_by_native.get(key, [])
        if eng:
            form = "engineered_construct; " + "; ".join(sorted({e["modification"] for e in eng}))
        record = {
            "assay_id": stable_id("ASSAY", row["sequence"], row["paper_doi"], row["wetlab_experiment"], row["mic"]),
            "compound_id": stable_id("CMP", row["sequence"], row["paper_doi"]),
            "sequence": row["sequence"],
            "peptide_name": peptide_name(row["wetlab_experiment"]),
            "chemical_form": form,
            "tested_sequence_note": "canonical string is not the complete tested chemistry" if form != "canonical_L_form_or_not_reported" else "",
            "assay_type": row["wetlab_experiment"].split(" [", 1)[0],
            "assay_target": assay_target(row["wetlab_experiment"]),
            "mic_relation": relation,
            "mic_value": value,
            "mic_unit": unit,
            "mic_source_text": row["mic"],
            "source_locator": bracket_locator(row["wetlab_experiment"]),
            "evidence_level": "exact_source_reported_MIC" if row["mic"] else "source_figure_positive_no_exact_MIC",
            "murine_validated": key in murine_keys,
            "paper_title": row["paper_title"],
            "paper_doi": row["paper_doi"],
            **similarity[row["sequence"]],
        }
        structured.append(record)
        by_key[key].append(record)

    struct_cols = [
        "assay_id", "compound_id", "sequence", "peptide_name", "chemical_form",
        "tested_sequence_note", "assay_type", "assay_target", "mic_relation", "mic_value",
        "mic_unit", "mic_source_text", "source_locator", "evidence_level", "murine_validated",
        "paper_title", "paper_doi", "exact_reference_match", "nearest_reference_sequence",
        "max_levenshtein_ratio", "passes_full_library_no_exact", "passes_top100_similarity",
    ]
    write_tsv(out / "df_wetlab_assays_structured.tsv", structured, struct_cols)

    records, row_mappings = build_source_records(spec, provenance)
    strict_gold = []
    decisions = []
    gold_source_ids = set()
    for key in sorted(by_key):
        assays = by_key[key]
        seq, doi = key
        mappings = row_mappings.get(key, [])
        forms = sorted({a["chemical_form"] for a in assays})
        is_engineered = key in engineered_by_native
        modified_only = all(f != "canonical_L_form_or_not_reported" for f in forms)
        eligible = []
        for spec_row, prov, sid, start in mappings:
            if is_engineered or modified_only:
                continue
            if len(spec_row["source_peptide"]) <= len(seq):
                continue
            if prov["distinct_parent_sequence_count"] != "1":
                continue
            eligible.append((spec_row, prov, sid, start))
        if eligible:
            decision, reason = "include", "exact natural fragment; unique longer parent; direct antimicrobial evidence"
        elif is_engineered:
            decision, reason = "exclude", "tested evidence is for an engineered construct, not the unmodified natural fragment"
        elif modified_only:
            decision, reason = "exclude", "activity is supported only for modified chemistry/stereochemistry"
        elif not mappings:
            decision, reason = "exclude", "no resolved biological parent mapping; commonly synthetic-design evidence"
        elif all(len(x[0]["source_peptide"]) <= len(seq) for x in mappings):
            decision, reason = "exclude", "peptide equals the complete source sequence; no cleavage boundary target"
        elif all(x[1]["distinct_parent_sequence_count"] != "1" for x in mappings):
            decision, reason = "exclude", "parent sequence is ambiguous"
        else:
            decision, reason = "exclude", "failed strict segmentation eligibility rules"

        decisions.append({
            "compound_id": assays[0]["compound_id"], "sequence": seq, "paper_doi": doi,
            "chemical_forms": " | ".join(forms), "assay_rows": len(assays),
            "murine_validated": any(a["murine_validated"] for a in assays),
            "segmentation_decision": decision, "decision_reason": reason,
            "resolved_parent_candidates": len(mappings), "eligible_parent_mappings": len(eligible),
            "challenge_exact_match": similarity[seq]["exact_reference_match"],
            "challenge_max_similarity": similarity[seq]["max_levenshtein_ratio"],
        })
        for spec_row, prov, sid, start in eligible:
            gold_source_ids.add(sid)
            numeric = [float(a["mic_value"]) for a in assays if a["mic_value"] and a["mic_unit"] == "µM"]
            strict_gold.append({
                "gold_label_id": stable_id("GOLD", seq, doi, sid, str(start)),
                "encrypted_peptide": seq,
                "source_id": sid,
                "source_accessions": prov["parent_accessions"],
                "source_peptide": spec_row["source_peptide"],
                "source_proteome": spec_row["source_proteome"],
                "peptide_start_0based": start,
                "peptide_end_0based_exclusive": start + len(seq),
                "parent_match_method": prov["parent_match_method"],
                "paper_title": assays[0]["paper_title"],
                "paper_doi": doi,
                "active_assay_rows": len(assays),
                "tested_targets": len({a["assay_target"] for a in assays}),
                "best_exact_mic_uM": min(numeric) if numeric else "",
                "murine_validated": any(a["murine_validated"] for a in assays),
                "evidence_tier": "strict_gold_exact_natural_fragment",
                "challenge_exact_match": similarity[seq]["exact_reference_match"],
                "challenge_max_similarity": similarity[seq]["max_levenshtein_ratio"],
                "challenge_top100_eligible": similarity[seq]["passes_top100_similarity"],
            })

    write_tsv(out / "df_wetlab_compound_decisions.tsv", decisions, [
        "compound_id", "sequence", "paper_doi", "chemical_forms", "assay_rows",
        "murine_validated", "segmentation_decision", "decision_reason",
        "resolved_parent_candidates", "eligible_parent_mappings", "challenge_exact_match",
        "challenge_max_similarity",
    ])

    membership = mmseqs_clusters(records, args.mmseqs)
    cluster_split = assign_splits(records, membership, gold_source_ids, args.seed)
    gold_cols = [
        "gold_label_id", "encrypted_peptide", "source_id", "source_accessions", "source_peptide",
        "source_proteome", "peptide_start_0based", "peptide_end_0based_exclusive",
        "parent_match_method", "paper_title", "paper_doi", "active_assay_rows", "tested_targets",
        "best_exact_mic_uM", "murine_validated", "evidence_tier", "homology_cluster",
        "recommended_split", "challenge_exact_match", "challenge_max_similarity",
        "challenge_top100_eligible",
    ]
    for row in strict_gold:
        row["homology_cluster"] = membership[row["source_id"]]
        row["recommended_split"] = "wetlab_holdout"
    strict_gold.sort(key=lambda r: (r["source_id"], int(r["peptide_start_0based"]), r["encrypted_peptide"]))
    write_tsv(out / "df_spear_segmentation_gold.tsv", strict_gold, gold_cols)

    split_rows = []
    split_dir = out / "segmentation_splits"
    split_dir.mkdir(parents=True, exist_ok=True)
    handles = {name: (split_dir / f"{name}.jsonl").open("w", encoding="utf-8")
               for name in ("train", "val", "test", "wetlab_holdout")}
    for sid, rec in sorted(records.items()):
        cid = membership[sid]
        split = cluster_split[cid]
        sites = sorted(rec["sites"])
        payload = {
            "source_id": sid,
            "source_accessions": sorted(rec["source_accessions"]),
            "source_sequence": rec["source_sequence"],
            "source_proteomes": sorted(rec["source_proteomes"]),
            "dois": sorted(rec["dois"]),
            "cleavage_sites": [[a, b] for a, b, _ in sites],
            "encrypted_peptides": [p for _, _, p in sites],
            "homology_cluster": cid,
            "split": split,
        }
        handles[split].write(json.dumps(payload, separators=(",", ":")) + "\n")
        split_rows.append({
            "source_id": sid,
            "source_accessions": ";".join(sorted(rec["source_accessions"])),
            "source_proteomes": ";".join(sorted(rec["source_proteomes"])),
            "sequence_length": len(rec["source_sequence"]),
            "peptide_sites": len(sites),
            "homology_cluster": cid,
            "split": split,
            "contains_strict_gold": sid in gold_source_ids,
            "contains_any_wetlab_label": rec["has_existing_wetlab_label"],
        })
    for fh in handles.values():
        fh.close()
    write_tsv(out / "df_spear_homology_splits.tsv", split_rows, [
        "source_id", "source_accessions", "source_proteomes", "sequence_length", "peptide_sites",
        "homology_cluster", "split", "contains_strict_gold", "contains_any_wetlab_label",
    ])

    gold_json = defaultdict(lambda: {"sites": set(), "peptides": set(), "dois": set()})
    for row in strict_gold:
        item = gold_json[row["source_id"]]
        item["sites"].add((int(row["peptide_start_0based"]), int(row["peptide_end_0based_exclusive"])))
        item["peptides"].add(row["encrypted_peptide"])
        item["dois"].add(row["paper_doi"])
    with (out / "spear_segmentation_gold.jsonl").open("w", encoding="utf-8") as fh:
        for sid, item in sorted(gold_json.items()):
            rec = records[sid]
            ordered = sorted((a, b, rec["source_sequence"][a:b]) for a, b in item["sites"])
            fh.write(json.dumps({
                "source_id": sid,
                "source_accessions": sorted(rec["source_accessions"]),
                "source_sequence": rec["source_sequence"],
                "cleavage_sites": [[a, b] for a, b, _ in ordered],
                "encrypted_peptides": [p for _, _, p in ordered],
                "paper_dois": sorted(item["dois"]),
                "homology_cluster": membership[sid],
                "split": "wetlab_holdout",
            }, separators=(",", ":")) + "\n")

    counts = {
        "assay_rows": len(structured),
        "unique_sequences": len({r["sequence"] for r in wetlab}),
        "compound_paper_records": len(decisions),
        "primary_papers": len({r["paper_doi"] for r in wetlab}),
        "exact_mic_rows": sum(bool(r["mic"]) for r in wetlab),
        "blank_mic_rows": sum(not bool(r["mic"]) for r in wetlab),
        "strict_gold_labels": len(strict_gold),
        "strict_gold_peptides": len({r["encrypted_peptide"] for r in strict_gold}),
        "strict_gold_sources": len(gold_source_ids),
        "source_sequences": len(records),
        "homology_clusters": len(set(membership.values())),
        "challenge_exact_matches": sum(r["exact_reference_match"] for r in challenge_rows),
        "challenge_similarity_over_0_8": sum(float(r["max_levenshtein_ratio"]) > 0.8 for r in challenge_rows),
        "split_source_counts": dict(Counter(r["split"] for r in split_rows)),
    }
    manifest = {
        "release": "SPEAR SOTA wet-lab and segmentation dataset",
        "schema_version": "1.0.0",
        "split_seed": args.seed,
        "homology_rule": "MMseqs2 min_seq_id=0.30, coverage=0.80, cov_mode=1; clusters kept intact",
        "challenge_reference_commit": args.challenge_commit,
        "challenge_reference_sha256": hashlib.sha256(args.challenge_fasta.read_bytes()).hexdigest(),
        "counts": counts,
    }
    (out / "release_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
