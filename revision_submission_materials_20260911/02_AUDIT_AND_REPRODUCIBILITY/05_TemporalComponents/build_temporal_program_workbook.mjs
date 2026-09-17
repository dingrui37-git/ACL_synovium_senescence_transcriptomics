import fs from "node:fs/promises";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const outputDir = "D:/workspace/ACLsenescence2_temporal_program_summary";
const csvDir = outputDir;
const outXlsx = `${outputDir}/Supplementary_Table_temporal_programs_genes_and_Reactome_KEGG.xlsx`;

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"') {
        if (text[i + 1] === '"') { field += '"'; i++; }
        else quoted = false;
      } else field += ch;
    } else if (ch === '"') quoted = true;
    else if (ch === ',') { row.push(field); field = ""; }
    else if (ch === "\n") { row.push(field.replace(/\r$/, "")); rows.push(row); row = []; field = ""; }
    else field += ch;
  }
  if (field.length || row.length) { row.push(field); rows.push(row); }
  const header = rows.shift();
  return rows.filter(r => r.length > 1 || (r.length === 1 && r[0] !== "")).map(r => {
    const out = {};
    header.forEach((h, i) => {
      // Remove XML-invalid control characters that can occur in legacy
      // provenance strings while retaining ordinary Unicode text.
      out[h] = (r[i] ?? "").replace(/[\u0000-\u0008\u000B\u000C\u000E-\u001F]/g, "");
    });
    return out;
  });
}

async function readCsv(name) {
  const text = await fs.readFile(`${csvDir}/${name}`, "utf8");
  return parseCsv(text);
}

const summary = await readCsv("Supplementary_Table_temporal_programs_summary.csv");
const genes = await readCsv("Supplementary_Table_temporal_programs_representative_genes.csv");
const top = await readCsv("Supplementary_Table_temporal_programs_top_pathways.csv");
const allOra = await readCsv("Supplementary_Table_temporal_programs_all_ORA_results.csv");
const qa = await readCsv("Supplementary_Table_temporal_programs_QA.csv");

const n = v => {
  if (v === "" || v === null || v === undefined) return null;
  const x = Number(v);
  return Number.isFinite(x) ? x : v;
};
const fdrSupported = top.filter(r => r.display_basis === "FDR-supported");
const names = (tc, dr, db, limit = 3) => fdrSupported
  .filter(r => r.temporal_class === tc && r.direction === dr && r.database === db)
  .sort((a, b) => Number(a.display_rank) - Number(b.display_rank))
  .slice(0, limit).map(r => r.pathway_name).join("; ");
const supportText = (tc, dr) => {
  const z = allOra.filter(r => r.temporal_class === tc && r.direction === dr && r.FDR_status === "FDR_supported");
  return z.length ? `${z.length} FDR-supported tested pathway rows` : "No FDR-supported pathway; nominal-only rows shown";
};

const summaryRows = summary.map(r => {
  const tc = r.temporal_class;
  const pattern = tc === "persistent"
    ? "Up: extracellular-matrix/collagen remodeling; down: lipid and amino-acid metabolism"
    : tc === "1-week-only"
      ? "Up: neutrophil/innate immune processes; down: mitochondrial respiration and oxidative phosphorylation"
      : "No FDR-supported pathway; only nominal signals in the smaller 4-week-only component";
  return [
    tc, n(r.strict_gene_count), n(r.up_count), n(r.down_count),
    r.representative_up_genes, r.representative_down_genes,
    names(tc, "Up", "Reactome"), names(tc, "Up", "KEGG"),
    names(tc, "Down", "Reactome"), names(tc, "Down", "KEGG"),
    supportText(tc, "Up"), supportText(tc, "Down"), pattern
  ];
});

const geneRows = genes.map(r => [
  r.temporal_class, r.direction, n(r.gene_rank), n(r.gene_id), r.gene_symbol,
  r.effect_timepoint, n(r.logFC_1W), n(r.FDR_1W), n(r.logFC_4W), n(r.FDR_4W), n(r.ranking_score)
]);

const topRows = top.map(r => [
  r.temporal_class, r.direction, r.database, n(r.display_rank), r.display_basis,
  r.pathway_id, r.pathway_name, n(r.pathway_size_in_universe),
  n(r.query_size_in_universe), n(r.overlap_size), n(r.expected_overlap),
  n(r.enrichment_ratio), n(r.p_value), n(r.FDR_BH_joint_Reactome_KEGG),
  r.overlap_gene_symbols
]);

const allRows = allOra.map(r => [
  r.temporal_class, r.direction, r.database, r.database_release, r.pathway_id,
  r.pathway_name, n(r.universe_size), n(r.query_size_in_universe),
  n(r.pathway_size_in_universe), n(r.overlap_size), n(r.expected_overlap),
  n(r.enrichment_ratio), n(r.p_value), n(r.FDR_BH_joint_Reactome_KEGG),
  r.FDR_status, r.overlap_gene_symbols
]);

const qaRows = qa.map(r => [r.metric, r.value]);

const wb = Workbook.create();
const summarySheet = wb.worksheets.add("Summary");
const geneSheet = wb.worksheets.add("Representative genes");
const topSheet = wb.worksheets.add("Top pathways");
const allSheet = wb.worksheets.add("All ORA results");
const qaSheet = wb.worksheets.add("Methods and QA");

const font = "Arial";
const navy = "#1F4E78";
const blue = "#D9EAF7";
const pale = "#F5F8FB";
const grey = "#5B6573";
const green = "#E2F0D9";
const amber = "#FFF2CC";

function styleSheet(sheet) {
  sheet.showGridLines = false;
  sheet.getUsedRange()?.format?.font;
}
function title(sheet, text, lastCol) {
  sheet.getRange(`A1:${lastCol}1`).merge();
  sheet.getRange("A1").values = [[text]];
  sheet.getRange(`A1:${lastCol}1`).format = {
    font: { name: font, size: 14, bold: true, color: navy },
    verticalAlignment: "center"
  };
  sheet.getRange(`A1:${lastCol}1`).format.rowHeight = 24;
}
function note(sheet, text, lastCol) {
  sheet.getRange(`A2:${lastCol}2`).merge();
  sheet.getRange("A2").values = [[text]];
  sheet.getRange(`A2:${lastCol}2`).format = {
    font: { name: font, size: 10, italic: true, color: grey },
    wrapText: true, verticalAlignment: "center"
  };
  sheet.getRange(`A2:${lastCol}2`).format.rowHeight = 34;
}
function header(sheet, range) {
  sheet.getRange(range).format = {
    fill: navy,
    font: { name: font, size: 10, bold: true, color: "#FFFFFF" },
    horizontalAlignment: "center", verticalAlignment: "center", wrapText: true,
    borders: { preset: "all", style: "thin", color: "#FFFFFF" }
  };
}
function body(sheet, range) {
  sheet.getRange(range).format = {
    font: { name: font, size: 10, color: "#1F2937" },
    verticalAlignment: "center", wrapText: true,
    borders: { preset: "insideHorizontal", style: "thin", color: "#D9E2F3" }
  };
}

// Summary
styleSheet(summarySheet);
title(summarySheet, "Supplementary Table: temporal mouse programs and pathway context", "M");
note(summarySheet, "Strict DEG threshold: FDR < 0.05 and |logFC| > 1. Persistent genes are significant at both 1W and 4W with consistent direction. ORA is run separately for up/down genes using the common two-timepoint eligibility universe; BH is pooled across 937 Reactome + 350 KEGG pathways within each temporal-class × direction list.", "M");
const sumHeaders = ["Temporal class", "Strict genes", "Up", "Down", "Representative up genes", "Representative down genes", "Top Reactome up", "Top KEGG up", "Top Reactome down", "Top KEGG down", "Up pathway support", "Down pathway support", "Interpretive pattern"];
summarySheet.getRange("A4:M4").values = [sumHeaders];
summarySheet.getRange("A5:M7").values = summaryRows;
header(summarySheet, "A4:M4"); body(summarySheet, "A5:M7");
summarySheet.getRange("A5:A7").format.font = { name: font, size: 10, bold: true, color: navy };
summarySheet.getRange("B5:D7").format.numberFormat = "#,##0";
summarySheet.getRange("A5:M7").conditionalFormats.add("expression", { formula: "=$A5=\"4-week-only\"", format: { fill: amber } });
summarySheet.getRange("A5:M7").conditionalFormats.add("expression", { formula: "=$A5=\"persistent\"", format: { fill: green } });
summarySheet.freezePanes.freezeRows(4);
summarySheet.getRange("A:A").format.columnWidth = 18;
summarySheet.getRange("B:D").format.columnWidth = 11;
summarySheet.getRange("E:F").format.columnWidth = 34;
summarySheet.getRange("G:J").format.columnWidth = 29;
summarySheet.getRange("K:L").format.columnWidth = 25;
summarySheet.getRange("M:M").format.columnWidth = 47;
summarySheet.getRange("A5:M7").format.rowHeight = 72;
summarySheet.tables.add("A4:M7", true, "TemporalProgramSummary");

// Representative genes
styleSheet(geneSheet);
title(geneSheet, "Representative genes by temporal class and direction", "K");
note(geneSheet, "Up to ten representatives per class × direction, ranked by −log10(FDR) × |logFC|; for persistent genes the ranking score is averaged across 1W and 4W. Values are shown for both time points to retain the derivation context.", "K");
const geneHeaders = ["Temporal class", "Direction", "Rank", "Entrez ID", "Gene symbol", "Effect time point", "logFC 1W", "FDR 1W", "logFC 4W", "FDR 4W", "Ranking score"];
geneSheet.getRange(`A4:K${4 + geneRows.length}`).values = [geneHeaders, ...geneRows];
header(geneSheet, "A4:K4"); body(geneSheet, `A5:K${4 + geneRows.length}`);
geneSheet.getRange(`C5:C${4 + geneRows.length}`).format.numberFormat = "0";
geneSheet.getRange(`D5:D${4 + geneRows.length}`).format.numberFormat = "0";
geneSheet.getRange(`G5:G${4 + geneRows.length}`).format.numberFormat = "0.000";
geneSheet.getRange(`H5:H${4 + geneRows.length}`).format.numberFormat = "0.000E+00";
geneSheet.getRange(`I5:I${4 + geneRows.length}`).format.numberFormat = "0.000";
geneSheet.getRange(`J5:J${4 + geneRows.length}`).format.numberFormat = "0.000E+00";
geneSheet.getRange(`K5:K${4 + geneRows.length}`).format.numberFormat = "0.0";
geneSheet.freezePanes.freezeRows(4);
geneSheet.getRange("A:A").format.columnWidth = 18;
geneSheet.getRange("B:B").format.columnWidth = 11;
geneSheet.getRange("C:D").format.columnWidth = 11;
geneSheet.getRange("E:E").format.columnWidth = 17;
geneSheet.getRange("F:F").format.columnWidth = 17;
geneSheet.getRange("G:K").format.columnWidth = 13;
geneSheet.tables.add(`A4:K${4 + geneRows.length}`, true, "RepresentativeGeneTable");

// Top pathway display
styleSheet(topSheet);
title(topSheet, "Top Reactome and KEGG pathways by temporal class and direction", "O");
note(topSheet, "Up to five pathways are shown per database × direction. FDR-supported rows are ranked first. Where no pathway survived pooled BH correction, the table displays top nominal rows and labels them as nominal-only.", "O");
const topHeaders = ["Temporal class", "Direction", "Database", "Display rank", "Display basis", "Pathway ID", "Pathway", "Pathway size", "Query size", "Overlap", "Expected overlap", "Enrichment ratio", "Nominal P", "Pooled BH-FDR", "Overlap genes"];
topSheet.getRange(`A4:O${4 + topRows.length}`).values = [topHeaders, ...topRows];
header(topSheet, "A4:O4"); body(topSheet, `A5:O${4 + topRows.length}`);
topSheet.getRange(`H5:K${4 + topRows.length}`).format.numberFormat = "0.0";
topSheet.getRange(`L5:L${4 + topRows.length}`).format.numberFormat = "0.00";
topSheet.getRange(`M5:N${4 + topRows.length}`).format.numberFormat = "0.000E+00";
topSheet.freezePanes.freezeRows(4);
topSheet.getRange("A:A").format.columnWidth = 18;
topSheet.getRange("B:C").format.columnWidth = 11;
topSheet.getRange("D:D").format.columnWidth = 11;
topSheet.getRange("E:E").format.columnWidth = 42;
topSheet.getRange("F:F").format.columnWidth = 14;
topSheet.getRange("G:G").format.columnWidth = 37;
topSheet.getRange("H:N").format.columnWidth = 14;
topSheet.getRange("O:O").format.columnWidth = 38;
topSheet.getRange(`A5:O${4 + topRows.length}`).conditionalFormats.add("expression", { formula: "=LEFT($E5,3)=\"FDR\"", format: { fill: green } });
topSheet.getRange(`A5:O${4 + topRows.length}`).conditionalFormats.add("expression", { formula: "=LEFT($E5,7)=\"top nom\"", format: { fill: amber } });
topSheet.tables.add(`A4:O${4 + topRows.length}`, true, "TopPathwayTable");

// Full tested ORA table
styleSheet(allSheet);
title(allSheet, "All tested Reactome and KEGG ORA results", "P");
note(allSheet, "Complete tested result set: 7,722 rows (6 temporal-class × direction lists × 1,287 pathways). Pooled BH-FDR is calculated separately within each of the six lists across Reactome and KEGG together.", "P");
const allHeaders = ["Temporal class", "Direction", "Database", "Database release", "Pathway ID", "Pathway", "Universe", "Query", "Pathway size", "Overlap", "Expected overlap", "Enrichment ratio", "Nominal P", "Pooled BH-FDR", "FDR status", "Overlap genes"];
allSheet.getRange(`A4:P${4 + allRows.length}`).values = [allHeaders, ...allRows];
header(allSheet, "A4:P4"); body(allSheet, `A5:P${4 + allRows.length}`);
allSheet.getRange(`G5:K${4 + allRows.length}`).format.numberFormat = "0.0";
allSheet.getRange(`L5:L${4 + allRows.length}`).format.numberFormat = "0.00";
allSheet.getRange(`M5:N${4 + allRows.length}`).format.numberFormat = "0.000E+00";
allSheet.freezePanes.freezeRows(4);
allSheet.getRange("A:A").format.columnWidth = 18;
allSheet.getRange("B:C").format.columnWidth = 11;
allSheet.getRange("D:D").format.columnWidth = 14;
allSheet.getRange("E:E").format.columnWidth = 14;
allSheet.getRange("F:F").format.columnWidth = 37;
allSheet.getRange("G:O").format.columnWidth = 14;
allSheet.getRange("P:P").format.columnWidth = 36;
allSheet.tables.add(`A4:P${4 + allRows.length}`, true, "AllOraResults");

// Methods and QA
styleSheet(qaSheet);
title(qaSheet, "Definitions, data sources and quality-control record", "B");
note(qaSheet, "This sheet documents the analysis contract used to produce the displayed tables. It is intended to keep the supplementary result auditable without replacing the manuscript Methods section.", "B");
qaSheet.getRange(`A4:B${4 + qaRows.length}`).values = [["Item", "Specification"], ...qaRows];
header(qaSheet, "A4:B4"); body(qaSheet, `A5:B${4 + qaRows.length}`);
qaSheet.getRange("A:A").format.columnWidth = 42;
qaSheet.getRange("B:B").format.columnWidth = 95;
qaSheet.getRange(`A5:B${4 + qaRows.length}`).format.rowHeight = 28;
qaSheet.tables.add(`A4:B${4 + qaRows.length}`, true, "MethodsQATable");

wb.recalculate();
const summaryInspect = await wb.inspect({ kind: "table", sheetId: "Summary", range: "A4:M7", include: "values,formulas", tableMaxRows: 8, tableMaxCols: 13, maxChars: 6000 });
await fs.writeFile(`${outputDir}/workbook_summary_inspect.ndjson`, summaryInspect.ndjson ?? String(summaryInspect), "utf8");
const preview = await wb.render({ sheetName: "Summary", autoCrop: "all", scale: 1, format: "png" });
await fs.writeFile(`${outputDir}/Supplementary_Table_temporal_programs_summary_preview.png`, new Uint8Array(await preview.arrayBuffer()));
const xlsx = await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(outXlsx);
console.log(`Saved ${outXlsx}`);
