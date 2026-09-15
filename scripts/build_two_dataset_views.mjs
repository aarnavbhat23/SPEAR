import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import crypto from "node:crypto";
import { execFileSync } from "node:child_process";
import { Workbook } from "/Users/aarnav/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";

const base = "data/spear_sota";

function parseTsv(text) {
  const matrix = text.trimEnd().split("\n").map((line) => line.split("\t"));
  const header = matrix[0];
  return matrix.slice(1).map((values) => Object.fromEntries(header.map((key, i) => [key, values[i] ?? ""])));
}

function clean(value) {
  return String(value ?? "").replace(/[\t\r\n]+/g, " ").trim();
}

function uniqueJoined(values) {
  return [...new Set(values.map(clean).filter(Boolean))].sort().join(" | ");
}

function stableId(prefix, ...parts) {
  const digest = crypto.createHash("sha256").update(parts.join("\x1f")).digest("hex").slice(0, 12);
  return `${prefix}_${digest}`;
}

async function peptideClusters(sequences) {
  const tmp = await fs.mkdtemp(path.join(os.tmpdir(), "spear_activity_mmseqs_"));
  try {
    const fasta = path.join(tmp, "peptides.fasta");
    const prefix = path.join(tmp, "clustered");
    const work = path.join(tmp, "work");
    await fs.writeFile(fasta, sequences.map((sequence) => `>${stableId("PEP", sequence)}\n${sequence}\n`).join(""));
    execFileSync("mmseqs", [
      "easy-cluster", fasta, prefix, work,
      "--min-seq-id", "0.60", "-c", "0.80", "--cov-mode", "1", "--threads", "4",
    ], { stdio: "ignore" });
    const idToSequence = new Map(sequences.map((sequence) => [stableId("PEP", sequence), sequence]));
    const membership = new Map();
    const lines = (await fs.readFile(`${prefix}_cluster.tsv`, "utf8")).trim().split("\n");
    for (const line of lines) {
      const [representative, member] = line.split("\t");
      membership.set(idToSequence.get(member), stableId("AHC", representative));
    }
    if (membership.size !== sequences.length) throw new Error("MMseqs2 did not assign every peptide");
    return membership;
  } finally {
    await fs.rm(tmp, { recursive: true, force: true });
  }
}

function assignClusterSplits(membership) {
  const sizes = new Map();
  for (const cluster of membership.values()) sizes.set(cluster, (sizes.get(cluster) ?? 0) + 1);
  const ordered = [...sizes].sort(([a], [b]) => stableId("ORDER", a).localeCompare(stableId("ORDER", b)));
  const targets = { train: 0.80 * membership.size, val: 0.10 * membership.size };
  const assignment = new Map();
  let assigned = 0;
  for (const [cluster, size] of ordered) {
    const split = assigned < targets.train ? "train" : assigned < targets.train + targets.val ? "val" : "test";
    assignment.set(cluster, split);
    assigned += size;
  }
  return assignment;
}

function normalizedTarget(value) {
  return clean(value)
    .replace(/\s+\([+-]\)(?=\s+-|$)/g, "")
    .replace(/ATCC\s*(\d)/g, "ATCC $1")
    .replace(/\(ATCC\s+([^)]+)\)/g, "ATCC $1")
    .replace(/\s+/g, " ")
    .trim();
}

const assays = parseTsv(await fs.readFile(`${base}/df_wetlab_assays_structured.tsv`, "utf8"));
const bySequence = new Map();
for (const row of assays) {
  if (!bySequence.has(row.sequence)) bySequence.set(row.sequence, []);
  bySequence.get(row.sequence).push(row);
}

const broadHeader = [
  "sequence", "wetlab_validated", "wetlab_experiments", "mic_measurements",
  "paper_titles", "paper_dois", "assay_row_count", "exact_mic_row_count",
  "chemical_forms", "murine_validated",
];
const broadRows = [...bySequence.entries()].sort(([a], [b]) => a.localeCompare(b)).map(([sequence, rows]) => [
  sequence,
  "True",
  uniqueJoined(rows.map((row) => row.assay_type)),
  uniqueJoined(rows.map((row) => row.mic_source_text)),
  uniqueJoined(rows.map((row) => row.paper_title)),
  uniqueJoined(rows.map((row) => row.paper_doi)),
  rows.length,
  rows.filter((row) => row.mic_source_text).length,
  uniqueJoined(rows.map((row) => row.chemical_form)),
  rows.some((row) => row.murine_validated === "True") ? "True" : "False",
]);
if (broadRows.length !== 596) throw new Error(`Expected 596 broad peptides, found ${broadRows.length}`);

const strictText = await fs.readFile(`${base}/df_spear_segmentation_gold.tsv`, "utf8");
const strictMatrix = strictText.trimEnd().split("\n").map((line) => line.split("\t"));
if (strictMatrix.length - 1 !== 270) throw new Error(`Expected 270 strict labels, found ${strictMatrix.length - 1}`);

const workbook = Workbook.create();
const broadSheet = workbook.worksheets.add("596 wet-lab peptides");
broadSheet.getRange("A1").write([broadHeader, ...broadRows]);
const strictSheet = workbook.worksheets.add("270 strict labels");
strictSheet.getRange("A1").write(strictMatrix);
workbook.recalculate();
const broadCheck = await workbook.inspect({
  kind: "table", range: "596 wet-lab peptides!A1:J6", include: "values,formulas",
  tableMaxRows: 6, tableMaxCols: 10,
});
const strictCheck = await workbook.inspect({
  kind: "table", range: "270 strict labels!A1:F6", include: "values,formulas",
  tableMaxRows: 6, tableMaxCols: 6,
});
console.log(broadCheck.ndjson);
console.log(strictCheck.ndjson);

const sequences = [...bySequence.keys()].sort((a, b) => a.localeCompare(b));
const membership = await peptideClusters(sequences);
const splitByCluster = assignClusterSplits(membership);
const peptideHeader = [
  "peptide_id", "sequence", "sequence_length", "chemical_forms", "assay_row_count",
  "distinct_target_count", "exact_mic_row_count", "exact_mic_uM_count",
  "exact_mic_ug_per_mL_count", "mic_label_status", "murine_validated", "paper_dois",
  "activity_homology_cluster", "activity_split", "challenge_exact_match",
  "challenge_max_similarity", "challenge_top100_eligible",
];
const peptideRows = sequences.map((sequence) => {
  const rows = bySequence.get(sequence);
  const cluster = membership.get(sequence);
  return [
    stableId("PEP", sequence), sequence, sequence.length,
    uniqueJoined(rows.map((row) => row.chemical_form)), rows.length,
    new Set(rows.map((row) => normalizedTarget(row.assay_target))).size,
    rows.filter((row) => row.mic_value && row.mic_relation === "=").length,
    rows.filter((row) => row.mic_value && row.mic_relation === "=" && row.mic_unit === "µM").length,
    rows.filter((row) => row.mic_value && row.mic_relation === "=" && row.mic_unit === "µg/mL").length,
    rows.some((row) => row.mic_value && row.mic_unit === "µM") ? "primary_uM_regression" :
      rows.some((row) => row.mic_value && row.mic_unit === "µg/mL") ? "secondary_mass_concentration_regression" :
        "positive_wetlab_no_numeric_MIC",
    rows.some((row) => row.murine_validated === "True") ? "True" : "False",
    uniqueJoined(rows.map((row) => row.paper_doi)), cluster, splitByCluster.get(cluster),
    rows[0].exact_reference_match, rows[0].max_levenshtein_ratio, rows[0].passes_top100_similarity,
  ];
});

const measurementHeader = [
  "assay_id", "peptide_id", "sequence", "chemical_form", "assay_type",
  "assay_target_raw", "assay_target_normalized", "target_id", "mic_relation",
  "mic_value", "mic_unit", "mic_log2_value", "mic_source_text", "regression_task",
  "regression_eligible", "primary_uM_regression_eligible", "sequence_task_sample_weight",
  "evidence_level", "murine_validated", "paper_title",
  "paper_doi", "activity_homology_cluster", "activity_split",
];
const exactMicCounts = new Map();
for (const sequence of sequences) {
  for (const unit of ["µM", "µg/mL"]) {
    exactMicCounts.set(`${sequence}\x1f${unit}`, bySequence.get(sequence).filter(
      (row) => row.mic_value && row.mic_relation === "=" && row.mic_unit === unit,
    ).length);
  }
}
const measurementRows = assays.map((row) => {
  const target = normalizedTarget(row.assay_target);
  const cluster = membership.get(row.sequence);
  const eligible = Boolean(row.mic_value && row.mic_relation === "=");
  const primaryEligible = eligible && row.mic_unit === "µM";
  const denominator = exactMicCounts.get(`${row.sequence}\x1f${row.mic_unit}`);
  const numericMic = eligible ? Number(row.mic_value) : NaN;
  const task = primaryEligible ? "MIC_uM" : eligible && row.mic_unit === "µg/mL" ? "MIC_ug_per_mL" : "positive_no_numeric_MIC";
  return [
    row.assay_id, stableId("PEP", row.sequence), row.sequence, row.chemical_form,
    row.assay_type, row.assay_target, target, stableId("TGT", target),
    row.mic_relation === "=" ? "equal" : row.mic_relation,
    row.mic_value, row.mic_unit, eligible ? Math.log2(numericMic).toFixed(8) : "",
    row.mic_source_text, task, eligible ? "True" : "False", primaryEligible ? "True" : "False",
    eligible ? (1 / denominator).toFixed(8) : "", row.evidence_level,
    row.murine_validated, row.paper_title, row.paper_doi, cluster, splitByCluster.get(cluster),
  ];
});

const peptideSheet = workbook.worksheets.add("596 activity peptides");
peptideSheet.getRange("A1").write([peptideHeader, ...peptideRows]);
const measurementSheet = workbook.worksheets.add("2449 activity assays");
measurementSheet.getRange("A1").write([measurementHeader, ...measurementRows]);
workbook.recalculate();
for (const range of ["596 activity peptides!A1:Q6", "2449 activity assays!A1:X6"]) {
  const check = await workbook.inspect({ kind: "table", range, include: "values,formulas", tableMaxRows: 6, tableMaxCols: 20 });
  console.log(check.ndjson);
}

const broadText = [broadHeader, ...broadRows].map((row) => row.map(clean).join("\t")).join("\n") + "\n";
await fs.writeFile(`${base}/SPEAR_596_wetlab_validated.tsv`, broadText);
await fs.writeFile(`${base}/SPEAR_270_strict_segmentation.tsv`, strictText);
await fs.writeFile(`${base}/SPEAR_596_activity_peptides.tsv`, [peptideHeader, ...peptideRows].map((row) => row.map(clean).join("\t")).join("\n") + "\n");
await fs.writeFile(`${base}/SPEAR_2449_activity_measurements.tsv`, [measurementHeader, ...measurementRows].map((row) => row.map(clean).join("\t")).join("\n") + "\n");
