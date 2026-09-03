// Optional companion workbook exporter; does not run simulations.
// Requires @oai/artifact-tool. Run with Node and the report directory argument.
import fs from 'node:fs/promises';
import path from 'node:path';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';
process.on('uncaughtException', error => {
  console.error(error.message);
  console.error(error.stack.split('\n').slice(-5).join('\n'));
  process.exit(1);
});

const root = path.resolve(process.argv[2] || '.');
const outputDir = path.join(root, 'outputs', '7beec1facbc5');
const qaDir = path.join(root, 'qa');
await fs.mkdir(outputDir, { recursive: true });
await fs.mkdir(qaDir, { recursive: true });
const wb = Workbook.create();
const inputs = wb.worksheets.add('Model');
inputs.showGridLines = false;
inputs.getRange('A1:C15').values = [
  ['Quantity', 'Value', 'Meaning / units'],
  ['Population mean', 100, 'IQ points; assumed'],
  ['Common SD', 15, 'IQ points; assumed'],
  ['Every pair correlation', 0.6, 'Population parameter; assumed'],
  ['Number of tests', 5, 'Five-test closed form below'],
  ['Expected maximum of 5 iid N(0,1)', null, 'Exact inverse-trigonometric formula'],
  ['Expected absolute pair gap', null, 'IQ points; jointly normal model'],
  ['Expected five-test range', null, 'IQ points; jointly normal model'],
  ['Difference SD / RMS', null, 'IQ points; equal means'],
  ['Source', 'Report: absolute gap and five-test range formulas', 'See accompanying PDF and R code'],
  ['Data status', 'Synthetic / theoretical', 'No measured IQ scores in this workbook'],
  ['Other sheets', 'Fixed R exports', 'They preserve the report run and do not update with inputs'],
  ['Recompute exports', 'Rscript code/iq_report_extensions.R', 'Then rerun this exporter'],
  ['Range source', 'https://arxiv.org/pdf/1111.4976', 'Finch, section 1'],
  ['CDF source', 'https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Tukey.html', 'R ptukey / qtukey; df=Inf']
];
inputs.getRange('B6:B9').formulas = [
  ['=5/(4*SQRT(PI()))*(1+6/PI()*ASIN(1/3))'],
  ['=2*B3*SQRT((1-B4)/PI())'],
  ['=IF(B5=5,2*B3*SQRT(1-B4)*B6,"Requires k=5")'],
  ['=B3*SQRT(2*(1-B4))']
];
inputs.getRange('A1:C15').format.font.name = 'Arial';
inputs.getRange('A1:C15').format.font.size = 11;
inputs.getRange('A1:C15').format.wrapText = true;
inputs.getRange('A1:A15').format.columnWidth = 38;
inputs.getRange('B1:B15').format.columnWidth = 48;
inputs.getRange('C1:C15').format.columnWidth = 52;
inputs.getRange('A1:C15').format.rowHeight = 38;
inputs.getRange('B2:B9').setNumberFormat('0.0000');
inputs.getRange('B5').setNumberFormat('0');
inputs.getRange('C15').clear({ applyTo: 'contents' });
inputs.getRange('B15:C15').merge();
inputs.getRange('B2:B5').format.font.color = '#176B88';
inputs.getRange('B6:B9').format.fill = '#E7F3F2';
inputs.getRange('A1:C1').format = {
  fill: '#183B56', font: { bold: true, color: '#FFFFFF' }, rowHeight: 30
};
inputs.freezePanes.freezeRows(1);

const exports = [
  ['probabilities', 'Threshold probabilities', 'Percent exceeding threshold; IQ-point thresholds'],
  ['quantiles', 'Quantiles', 'Probability; IQ-point quantiles'],
  ['sensitivity', 'Correlation sensitivity', 'Correlation; expected IQ-point gaps'],
  ['test_counts', 'Number of tests', 'Test count; standard-normal maximum; IQ-point range'],
  ['conditional', 'Conditional scores', 'All score and gap values in IQ points'],
  ['models', 'Alternative models', 'Expected gaps in IQ points'],
  ['smooth_simulation', 'Smooth mixture run', 'IQ-point means and MC SE; dimensionless z_mc'],
  ['figure_correlation', 'Correlation curve', 'Exact values plotted in correlation.pdf'],
  ['figure_distributions', 'Distribution curves', 'IQ-point gaps and CDF probabilities'],
  ['figure_conditional', 'Conditional curves', 'IQ-point first scores and expected gaps']
];
for (const [file, sheetName, units] of exports) {
  const csv = await fs.readFile(path.join(root, 'tables', file + '.csv'), 'utf8');
  // Import into a fresh workbook, then transfer typed values. Importing CSV
  // into an already populated collaborative workbook is not supported.
  const imported = await Workbook.fromCSV(csv, { sheetName });
  const sheet = wb.worksheets.add(sheetName);
  sheet.getRange('A1').writeValues(imported.worksheets.getItemAt(0).getUsedRange().values);
  sheet.showGridLines = false;
  const used = sheet.getUsedRange();
  used.format.font.name = 'Arial';
  used.format.font.size = 11;
  used.format.columnWidth = file === 'models' ? 35 : 27;
  used.format.rowHeight = 24;
  used.setNumberFormat('0.000000');
  if (file === 'test_counts' || file === 'probabilities') {
    sheet.getRange('A2:A' + csv.trim().split('\n').length).setNumberFormat('0');
  }
  const count = csv.trim().split('\n')[0].split(',').length;
  sheet.getRangeByIndexes(0, 0, 1, count).format = {
    fill: '#183B56', font: { bold: true, color: '#FFFFFF' },
    wrapText: true, rowHeight: 48
  };
  const rowCount = csv.trim().split('\n').length;
  sheet.getCell(rowCount + 1, 0).values = [['Units / provenance']];
  sheet.getCell(rowCount + 2, 0).values = [[units]];
  sheet.getRangeByIndexes(rowCount + 2, 0, 1, count).merge();
  sheet.getRangeByIndexes(rowCount + 2, 0, 1, count).format.wrapText = true;
  sheet.getRangeByIndexes(rowCount + 2, 0, 1, count).format.rowHeight = 34;
  sheet.getCell(rowCount + 3, 0).values = [['Source: code/iq_report_extensions.R; model calculations, 2026-09-03']];
  sheet.getRangeByIndexes(rowCount + 3, 0, 1, count).merge();
  sheet.getRangeByIndexes(rowCount + 3, 0, 1, count).format.wrapText = true;
  sheet.getRangeByIndexes(rowCount + 3, 0, 1, count).format.rowHeight = 34;
  sheet.freezePanes.freezeRows(1);
}
console.log((await wb.inspect({ kind: 'table', range: 'Model!A1:C9',
  include: 'values,formulas', tableMaxRows: 9, tableMaxCols: 3,
  maxChars: 3000 })).ndjson);
console.log((await wb.inspect({ kind: 'match',
  searchTerm: '#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A',
  options: { useRegex: true, maxResults: 30 }, maxChars: 1000,
  summary: 'Formula error scan' })).ndjson);
for (const name of ['Model', ...exports.map(x => x[1])]) {
  const sheet = wb.worksheets.getItem(name);
  const last = name === 'Model' ? 'C15' : name === 'Smooth mixture run' ? 'E8' : 'D8';
  const blob = await wb.render({ sheetName: name, range: 'A1:' + last, scale: 1.5 });
  await fs.writeFile(path.join(qaDir, 'workbook_' + name.replaceAll(' ', '_') + '.png'),
    new Uint8Array(await blob.arrayBuffer()));
}
await (await SpreadsheetFile.exportXlsx(wb)).save(path.join(outputDir, 'iq_report_data.xlsx'));
console.log('Saved iq_report_data.xlsx');
