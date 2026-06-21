# Makefile -- Linux/Mac equivalents of run.ps1 / setup.ps1.
# Thin wrappers around run_pipeline.py; the venv interpreter is used if present.
PY := $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)

.PHONY: setup check qc integration all test menu clean-smoke

setup:
	python3 -m venv .venv
	.venv/bin/python -m pip install --upgrade pip
	.venv/bin/python -m pip install -r requirements.txt
	@echo "Setup done. Verify with: make check"

check:
	$(PY) run_pipeline.py --check

qc:
	$(PY) run_pipeline.py --stage qc $(FLAGS)

integration:
	$(PY) run_pipeline.py --stage integration $(FLAGS)

annotation:
	$(PY) run_pipeline.py --stage annotation $(FLAGS)

composition:
	$(PY) run_pipeline.py --stage composition $(FLAGS)

de:
	$(PY) run_pipeline.py --stage de $(FLAGS)

size:
	$(PY) run_pipeline.py --stage size $(FLAGS)

# Run every registered stage in order. Keep in sync with STAGE_REGISTRY.
all: qc integration annotation composition de size

test:
	$(PY) -m pytest

menu:
	$(PY) run_pipeline.py

# Convenience: SMOKE=1 DEBUG=1 make qc  ->  passes the flags through.
ifeq ($(SMOKE),1)
FLAGS += --smoke-test
endif
ifeq ($(DEBUG),1)
FLAGS += --debug
endif

clean-smoke:
	rm -f data/processed/*_smoke.h5ad
