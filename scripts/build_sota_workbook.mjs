import fs from "node:fs/promises";
import { Workbook } from "/Users/aarnav/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";

const base = "data/spear_sota";
const workbook = Workbook.create();

async function matrix(name) {
  const text = await fs.readFile(`${base}/${name}`, "utf8");
  return text.trimEnd().split("\n").map((line) => line.split("\t"));
}

function styleSheet(sheet, data, widths) {
  const cols = data[0].length;
  const end = String.fromCharCode(64 + Math.min(cols, 26));
  sheet.getRange(`A1:${end}1`).format = {
    fill: "#17365D",
    font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" },
  };
  sheet.getRange(`A1:${end}${data.length}`).format.verticalAlignment = "center";
  sheet.getRange(`A2:${end}${data.length}`).format.font = { name: "Arial", size: 9 };
  widths.forEach((width, index) => {
    sheet.getRange(`${String.fromCharCode(65 + index)}:${String.fromCharCode(65 + index)}`).format.columnWidth = width;
  });
  sheet.freezePanes.freezeRows(1);
}

const overview = workbook.worksheets.add("README");
const summary = [
  ["SPEAR SOTA release", "Use the strict gold set as a wet-lab holdout; train on the homology-safe JSONL split."],
  ["Assay rows", 2449],
  ["Unique canonical sequences", 596],
  ["Primary papers", 22],
  ["Strict segmentation labels", 270],
  ["Strict gold parent sequences", 217],
  ["All parent sequences", 12278],
  ["Homology clusters", 3811],
  ["Challenge exact matches", 23],
  ["Challenge similarity > 0.80", 93],
  ["Important", "The sequence-only field cannot encode amidation, D stereochemistry, or engineered termini. Use the structured assay table."],
];
overview.getRange("A1:B11").write(summary);
overview.getRange("A1:B1").format = { fill: "#17365D", font: { bold: true, color: "#FFFFFF", size: 12 } };
overview.getRange("A:A").format.columnWidth = 34;
overview.getRange("B:B").format.columnWidth = 96;
overview.getRange("A1:B11").format.wrapText = true;

for (const [file, sheetName, widths] of [
  ["df_wetlab_validated_sota.tsv", "Vinay 5-column", [36, 80, 14, 62, 30]],
  ["df_spear_segmentation_gold.tsv", "Segmentation gold", [24, 34, 24, 34, 72, 22, 14, 14, 30, 62, 30]],
  ["df_wetlab_compound_decisions.tsv", "Inclusion decisions", [24, 34, 30, 50, 12, 16, 18, 72, 16, 16, 18, 18]],
]) {
  const data = await matrix(file);
  const sheet = workbook.worksheets.add(sheetName);
  sheet.getRange("A1").write(data);
  styleSheet(sheet, data, widths);
}

workbook.recalculate();
const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!",
  options: { useRegex: true, maxResults: 20 },
  summary: "formula error scan",
});
console.log(errors.ndjson);
const preview = await workbook.render({ sheetName: "README", range: "A1:B11", scale: 1.3, format: "png" });
await fs.writeFile(`${base}/SPEAR_SOTA_preview.png`, new Uint8Array(await preview.arrayBuffer()));
const xlsx = await workbook.export({ format: "xlsx" });
await fs.writeFile(`${base}/SPEAR_SOTA_release.xlsx`, new Uint8Array(await xlsx.arrayBuffer()));
