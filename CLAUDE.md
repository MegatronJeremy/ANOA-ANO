# CLAUDE.md — radni okvir za ovaj projekat

scRNA-seq pipeline: imuni odgovor ljudskih PBMC ćelija na nanoplastiku (PSNP) u funkciji
veličine čestice. Driver je `run_pipeline.py` (staged, `--stage`/meni), faze u `src/`.
Detaljan razvojni roadmap je u `PROMPTOVI.md`.

## Standardni kvalitetski okvir (poštuj na svakoj izmeni)

1. **Jedan inkrement po koraku.** Svaka faza se završava u pokretljivom stanju + checkpoint.
2. **Jedan code path.** Svaka faza je jedan red `StageSpec` u `STAGE_REGISTRY`; `--stage` i
   meni zovu istu `run_fn`. Bez paralelnih implementacija.
3. **Smoke-test paritet.** Svaka faza prolazi end-to-end na `--smoke-test` (~500 ćelija/uzorak)
   za sekunde i piše `*_smoke` izlaz. Smoke i full su isti kod.
4. **`--debug` priča istinu o podacima.** Timing banneri, AnnData shape pre/posle, rich tabele
   i plotext terminalni plotovi. Pragovi se opravdavaju odštampanim brojevima, ne nagađaju.
   Koristi helpere iz `src/logging_utils.py` (`log_stage_banner`, `log_table`, `describe_adata`,
   `log_percentile_table`, `term_hist/bar/scatter`, `progress_spinner`, `require`).
5. **Reproduktivnost.** `cfg.RANDOM_SEED` svuda (numpy/scanpy/harmony/leiden/umap); headless
   `--stage` je kanonski put; checkpoint-ovi u `data/processed/`.
6. **Testovi namere (pytest), offline.** Testiraj logiku (set-operacije, normalizaciju,
   pragove, harmonizaciju labela) na sitnim sintetičkim AnnData objektima — brzo, determinisno,
   bez teških podataka i bez mreže.
7. **Differential verification.** Gde POSTOJI nezavisna referenca, ukrsti: anotacija (Stage 3)
   vs `obs['predicted.celltype']` (Azimuth/R) **i** vs `*_CoDi_KLD.csv` (CoDi metod). Za
   DE/sastav nema eksterne reference — validacija je interna. Ne izmišljaj referencu koje nema.
8. **Biološka interpretacija.** Svaki rezultat dobija figuru/tabelu u `results/` + 2–4 rečenice
   „šta ovo biološki znači".
9. **Jezik:** kod, komentari i izlaz na **engleskom**; komunikacija sa Vukom na srpskom.
10. **Verify-before-claim.** Pre nego što predaš artefakt ili tvrdnju o kodu/podacima, proveri
    svaku netrivijalnu činjenicu naspram **stvarnog fajla/komande** (otvori, `grep`, `head`,
    pokreni), ne iz imena ili pretpostavke. Što ne možeš potvrditi — označi ASSUMED. Bar jedan
    ground-truth prolaz uvek ide PRE isporuke.

## Ograničenja specifična za ovaj dataset (ne zaboravi)

- **Jedan donor, jedan uzorak po uslovu — NEMA bioloških replika.** Zato DE = cell-level
  Wilcoxon (`rank_genes_groups`), NE pseudobulk/DESeq2; sastav = deskriptivni udeli +
  fold-change, svaki p-value eksploratoran uz caveat o pseudoreplikaciji.
- **`.h5ad` su Seurat→AnnData konverzije.** `.X` je Seurat LogNormalize (scale.factor=1000),
  `.layers['counts']` su sirovi count-ovi (Stage 1 resetuje `.X` na counts i sve recomputuje
  scanpy-jem). `.obs['predicted.celltype']` je Azimuth iz R-a — referenca, ne ulaz.
- **`*_CoDi_KLD.csv` NIJE gen-divergencija** (ime „KLD" zavarava) — to je per-ćelija anotacija
  tipa ćelije iz CoDi metoda (`cell_id, CoDi, CoDi_dist, CoDi_contrastive, ...`). Druga
  referenca za anotaciju (Stage 3), NE referenca za DE.
- Faze 3–5 logike (legacy) su u `legacy/notebooks/*.ipynb` — portuju se u čiste `src/` module.
  Faze 4–6 proizvode tabele/figure, NE nove 2.7 GB h5ad checkpoint-e.

## Pokretanje

```
.\run.ps1                 # interaktivni meni
.\run.ps1 check           # environment doctor (--check)
.\run.ps1 qc -Smoke -Debug
python run_pipeline.py --stage qc --smoke-test --debug   # ekvivalent, headless
.\run.ps1 test            # pytest
```
