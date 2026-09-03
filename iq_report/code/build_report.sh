#!/usr/bin/env bash
set -euo pipefail
# Run from the extracted report directory.
# Uses supplied figures and tables; does not overwrite simulation logs.
latexmk -xelatex -interaction=nonstopmode -halt-on-error iq_test_differences_report.tex
# To regenerate derived figures/tables first:
# Rscript code/iq_report_extensions.R
# To repeat the main simulation, choose a NEW output filename:
# Rscript code/iq_difference_simulation.R --out=results/main_run_new.txt
