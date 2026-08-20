# OSD-682/685/698/699 GeoMx DSP brain re-analysis - reproducible pipeline
# Usage:  make setup && make all
PY := ./.venv/bin/python

.PHONY: all setup data conventions targets analyze panels tables summary clean clean-results

## Full pipeline (assumes `make setup` already run)
all: data conventions analyze panels tables summary

## Create the virtual environment and install pinned dependencies.
## Prefers uv (fast, and can fetch Python 3.12); falls back to the stdlib venv.
setup:
	@if command -v uv >/dev/null 2>&1; then \
		uv venv --python 3.12 .venv && \
		uv pip install --python $(PY) -r requirements.txt; \
	else \
		python3 -m venv .venv && \
		$(PY) -m pip install --upgrade pip && \
		$(PY) -m pip install -r requirements.txt; \
	fi

## Download the processed tables from GEO + ISA metadata from OSDR (~12 MB)
data:
	$(PY) scripts/download_data.py

## Resolve the deposited log2 sign convention from the data (writes config/)
conventions:
	$(PY) scripts/check_conventions.py

## Map the CNS target spreadsheet onto mouse GeoMx targets (needs network once)
targets:
	$(PY) scripts/map_targets.py

## QC, global structure, differential expression, target panel, recurrence
analyze: targets
	$(PY) scripts/analyze_qc.py
	$(PY) scripts/analyze_expression.py
	$(PY) scripts/analyze_de.py
	$(PY) scripts/analyze_targets.py
	$(PY) scripts/analyze_recurrence.py

## Compose the numbered main figures -> results/figures/panels/
panels:
	$(PY) scripts/make_panels.py

## Export every result table (CSV + Excel) -> results/tables/
tables:
	$(PY) scripts/export_results.py

## Dataset + results summary (markdown + json) -> results/tables/
summary:
	$(PY) scripts/dataset_summary.py

## Remove generated results (keeps downloaded data)
clean-results:
	rm -rf results/figures results/tables config/sign_convention.json

## Remove results AND downloaded data
clean: clean-results
	rm -rf data
