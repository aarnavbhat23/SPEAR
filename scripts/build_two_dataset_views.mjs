import fs from "node:fs/promises";
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

const broadText = [broadHeader, ...broadRows].map((row) => row.map(clean).join("\t")).join("\n") + "\n";
await fs.writeFile(`${base}/SPEAR_596_wetlab_validated.tsv`, broadText);
await fs.writeFile(`${base}/SPEAR_270_strict_segmentation.tsv`, strictText);

