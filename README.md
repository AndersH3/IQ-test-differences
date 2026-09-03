# IQ test differences

Exact mathematical results, reproducible R simulations, and a research report
about differences between correlated IQ scores.

**AI disclosure:** ChatGPT/OpenAI generated the research synthesis, mathematical
exposition, source code, figures, and report at the request of Anders Hellström.
All simulations are synthetic. No observed IQ-test dataset was analysed and no
independent human peer review is claimed. The report contains a prominent opening
disclosure and a complete AI-use appendix.

## Read the report

- [Full report (PDF, 38 pages)](iq_report/iq_test_differences_report.pdf)
- [Companion workbook (XLSX)](iq_report/outputs/7beec1facbc5/iq_report_data.xlsx)
- [Original mathematical proof](iq_difference_proof.md)
- [Report reproduction instructions](iq_report/README.md)
- [Source and evidence ledger](source-ledger.md)

Assuming equal means, standard deviation 15, joint normality, and Pearson
correlation 0.6 for every distinct pair of tests:

| Quantity | Expected value, IQ points |
|---|---:|
| Signed difference between two tests | 0 |
| Absolute difference between two tests | 10.7047447 |
| Highest minus lowest of five tests | 22.0656994 |

Correlation alone does not determine the absolute difference or range. The
report proves these results and gives counterexamples, conditional calculations,
sensitivity analyses, and limits on their application to actual tests.

## Contents

| Path | Contents |
|---|---|
| `iq_difference_simulation.R` | Original self-contained R simulation |
| `iq_difference_proof.md` | Original proof, assumptions, and limitations |
| `iq_simulation_results.txt` | Original recorded five-million-person run |
| `iq_report/code/` | Both R programs, PDF build script, and optional workbook exporter |
| `iq_report/*.tex`, `iq_report/references.bib` | XeLaTeX/KOMA-Script source and 21 references |
| `iq_report/tables/` | CSV data, exact plotted values, and LaTeX table fragments |
| `iq_report/figures/` | Three generated PDF figures |
| `iq_report/results/` | Recorded simulation and extension outputs |
| `iq_report/outputs/` | Workbook and workbook inspection output |
| `iq_report/qa/` | Retained report/workbook previews and validation records |
| `validation/edge_cases/` | Recorded checks for identical, independent, two-test, and ten-test cases |
| `report-source.md`, `source-ledger.md` | Original assembly snapshot and research provenance notes |
| `iq_test_differences_source_package.zip` | Source package matching the current report |
| `PROJECT_MANIFEST.sha256` | Checksums of the complete project snapshot |

The report PDF, editable sources, and source package are kept in sync; earlier
versions remain available in Git history. Some files intentionally duplicate
earlier deliverables to preserve those snapshots. `report-source.md` contains the original TeX body
despite its filename. The ledger retains working retrieval references; the
bibliography supplies the source URLs.

Temporary TeX compilation files and downloaded operating-system installation
packages are excluded. The report PDF, figures, data, run logs, workbook, and
validation records are included.

## Reproduce the calculations

Run these commands from the repository root. The R programs use only packages
supplied with R; the recorded report runs used R 4.3.3 on 64-bit Linux.

```bash
cd iq_report
Rscript code/iq_difference_simulation.R --out=results/main_run_new.txt
Rscript code/iq_report_extensions.R
```

The original simulation refuses to overwrite an existing result file, so use a
new output filename for each saved run. The extension program is independent of
the original simulation: it reads no input files and creates or replaces its
tables, figures, and `results/extensions.txt` in the working directory.

The extension script uses `cairo_pdf()` and requests DejaVu Sans. Check R's
graphics support with:

```bash
Rscript -e 'capabilities("cairo")'
```

This should report `TRUE`. No separate CRAN Cairo package is used.

## Rebuild the PDF

From the repository root:

```bash
cd iq_report
bash code/build_report.sh
```

Requires XeLaTeX, `latexmk`, BibTeX, KOMA-Script, the packages in the TeX preamble,
Latin Modern text/math fonts, and DejaVu Sans Mono. Supplied tables and figures
allow PDF compilation without rerunning the simulations.

## Optional workbook export

From the repository root, in an environment that provides Node.js and
`@oai/artifact-tool`:

```bash
node iq_report/code/export_workbook.mjs iq_report
```

This exporter depends on the original OpenAI artifact tooling; its availability
in an ordinary Node installation is not assumed. The supplied workbook can be
opened directly, and all numeric source data are also supplied as CSV. Its R
export sheets are fixed snapshots. To update those sheets, rerun the extension
program and then the exporter.

## Check snapshot integrity

From the repository root:

```bash
sha256sum -c PROJECT_MANIFEST.sha256
```

The report package also has its own `iq_report/MANIFEST.sha256`.
Checksums describe the saved snapshot; intentional edits or regenerated output
will change them.

## Git repository

The project is hosted at
[AndersH3/IQ-test-differences](https://github.com/AndersH3/IQ-test-differences),
with `main` as the default branch. Clone the complete project with:

```bash
git clone https://github.com/AndersH3/IQ-test-differences.git
cd IQ-test-differences
```

## Licence

The repository includes the [CC0 1.0 Universal dedication](LICENSE) selected
when the GitHub repository was created.
