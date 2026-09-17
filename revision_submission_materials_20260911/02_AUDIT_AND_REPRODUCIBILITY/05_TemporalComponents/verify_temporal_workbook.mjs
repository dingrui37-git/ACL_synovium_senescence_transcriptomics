import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";
const path = "D:/workspace/ACLsenescence2_temporal_program_summary/Supplementary_Table_temporal_programs_genes_and_Reactome_KEGG.xlsx";
const wb = await SpreadsheetFile.importXlsx(await FileBlob.load(path));
const sheets = await wb.inspect({ kind: "sheet", include: "id,name", maxChars: 2000 });
console.log(sheets.ndjson ?? String(sheets));
const check = await wb.inspect({ kind: "table", sheetId: "Summary", range: "A4:M7", include: "values,formulas", tableMaxRows: 8, tableMaxCols: 13, maxChars: 6000 });
console.log(check.ndjson ?? String(check));
const errors = await wb.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!", options: { useRegex: true, maxResults: 50 }, summary: "formula error scan" });
console.log(errors.ndjson ?? String(errors));
