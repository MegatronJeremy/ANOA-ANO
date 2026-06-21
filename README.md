# Single-Cell Analysis of Immune Response to Nanoplastic Particles

scRNA-seq analysis of human PBMCs exposed to carboxylated polystyrene nanoparticles (PSNPs),
investigating how particle size shapes the immune response.

Genomics Informatics project — ETF Belgrade.

## The question

Four samples from one donor:

| Label | Exposure |
|---|---|
| `PSNP_40nm` | 40 nm particles |
| `PSNP_200nm` | 200 nm particles |
| `PSNP_mixture` | 40 nm + 200 nm |
| `control` | none (reference) |

Do small vs. large nanoplastics provoke different responses in immune cells, and does the
mixture do something neither size does alone?

## Setup

No conda on the dev machine — built and verified with a plain venv. `requirements.txt` is the
environment file.

**Windows (recommended) — one-shot setup script:**

```powershell
.\setup.ps1                       # finds Python 3.10+, makes .venv, installs deps, checks raw data
.\run.ps1 check                   # environment doctor: verify everything is ready
```

**Manual / Linux / Mac:**

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt   # or: make setup
```

Data: the 4 `.h5ad` files (and optional `*_CoDi_KLD.csv` reference files) go in `data/raw/`.
See `src/config.py` for the exact filename → sample mapping.

## Running the pipeline

Everything goes through `run_pipeline.py`. Each stage reads the previous stage's checkpoint
from `data/processed/` and writes its own. On Windows, `run.ps1` is a thin wrapper around it
(new stages appear automatically as they are registered):

```powershell
.\run.ps1                          # interactive menu
.\run.ps1 check                    # environment doctor (= run_pipeline.py --check)
.\run.ps1 qc -Smoke -Debug         # run a stage on the smoke subsample, verbose
.\run.ps1 all                      # every registered stage in order
.\run.ps1 test                     # pytest
```

The headless `run_pipeline.py` calls below are the canonical, scriptable entry point
(`make qc` / `make check` are the Linux/Mac equivalents):

```bash
# headless, scriptable -- the canonical entry point reproducibility/grading depends on
python run_pipeline.py --stage qc --smoke-test --debug     # smoke test: ~500 cells/sample, seconds
python run_pipeline.py --stage qc --debug                    # full run on real data
python run_pipeline.py --stage qc                              # quiet mode

python run_pipeline.py --stage integration --smoke-test --debug
python run_pipeline.py --stage integration --debug             # ~7-8 min on the full dataset

# interactive launcher (optional convenience, not a separate code path):
python run_pipeline.py                                          # menu, only in a real terminal;
                                                                  # piped/CI invocations just print help
```

`--debug` turns on rich per-stage logging: timing banners, AnnData shape before/after, structured
tables (QC filter summary/percentiles, cluster sizes, per-sample counts), and terminal
histograms/scatter plots in a real terminal session (QC threshold distributions, before/after
integration UMAP, cluster sizes). `--subsample N` gives a custom smoke-test size instead of the
default. The interactive menu shows per-stage status (input ready? already run?) and calls the
exact same stage functions as the `--stage` flag.

## Status

| Stage | Status |
|---|---|
| 1. QC & preprocessing | done (`--stage qc`) |
| 2. Integration & clustering | done (`--stage integration`) |
| 3. Cell-type annotation | done (`--stage annotation`) — celltypist + lineage, cross-checked vs Azimuth (92.7%) & CoDi (93.1%) |
| 4. Composition analysis | in progress |
| 5. Differential expression + pathway enrichment | in progress |
| 6. Size-specific effects | in progress |

## Repo layout

```
.
├── requirements.txt      # pip environment (venv-based)
├── README.md
├── run_pipeline.py        # driver: stage registry, --stage flags, interactive menu
├── src/
│   ├── config.py          # paths, sample names, thresholds, markers, seed
│   ├── io.py               # load samples / save+load checkpoints
│   ├── logging_utils.py    # --debug logging, rich tables/panels, terminal plots, guard checks
│   ├── qc.py                # Stage 1: QC & preprocessing
│   └── integration.py        # Stage 2: Harmony integration, UMAP, Leiden clustering
├── data/                   # raw/ and processed/ (git-ignored)
├── results/                 # figures/ and tables/ (git-ignored)
└── legacy/                  # old notebook-based pipeline, local reference only (git-ignored)
```

A full README (setup detail, per-stage explanations, bonus analyses, deliverables checklist)
will be written once all 6 stages are built and verified.
