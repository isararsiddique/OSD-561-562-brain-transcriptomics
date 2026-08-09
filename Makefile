# OSD-561/562 CNS spaceflight re-analysis - reproducible pipeline
# Usage:  make setup && make all
PY := ./.venv/bin/python

.PHONY: all setup data analyze panels tables summary clean clean-results

## Full pipeline (assumes `make setup` already run)
all: data analyze panels tables summary

## Create the virtual environment and install pinned dependencies
setup:
	python3 -m venv .venv
	$(PY) -m pip install --upgrade pip
	$(PY) -m pip install -r requirements.txt

## Download the GeneLab processed tables from NASA OSDR (~196 MB)
data:
	$(PY) scripts/download_data.py

## Core + full-data analysis -> individual figures in results/figures/
analyze:
	$(PY) scripts/analyze.py
	$(PY) scripts/analyze_full.py

## Compose the multi-panel figures -> results/figures/panels/
panels:
	$(PY) scripts/make_panels.py

## Export every result table (CSV + Excel) -> results/tables/
tables:
	$(PY) scripts/export_results.py

## Dataset summary (markdown + json) -> results/tables/
summary:
	$(PY) scripts/dataset_summary.py

## Remove generated results (keeps downloaded data)
clean-results:
	rm -rf results/figures results/tables

## Remove results AND downloaded data
clean: clean-results
	rm -rf data
