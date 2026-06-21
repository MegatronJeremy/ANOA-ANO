# GI projekat — promptovi (razvojni roadmap)

Isti pristup kao na PSZ-u (Sea-Of-Sorrow): projekat se gradi kroz **niz fokusiranih
faznih prompt-ova**, gde je svaki jedan proverljiv inkrement koji se završava u stanju
koje možeš da pokreneš i pogledaš. Nije „napiši ceo projekat" — nego lanac koraka koji
svaki put ostavlja zelen build, smoke-test koji prolazi i checkpoint za sledeću fazu.

Svaki prompt ispod je **paste-ready** — kopiraš ga u Claude Code kao jednu poruku.
Redom. Ne preskači verifikaciju između koraka.

---

## Postavka projekta (rezime na koji se promptovi oslanjaju)

**Tema:** Single-cell (scRNA-seq) analiza imunog odgovora ljudskih PBMC ćelija na
nanoplastiku (polistirenske čestice, PSNP), u funkciji veličine čestice.

**Podaci** (`data/raw/`, AnnData `.h5ad`, Zenodo DOI 10.5281/zenodo.15866724) — 4 uzorka
jednog donora: `PSNP_40nm` (S1), `PSNP_200nm` (S2), `PSNP_mixture` (S3), `control` (S4).
Bitno (već zabeleženo u `src/config.py`): `.h5ad` su Seurat→AnnData konverzije; `.X` je
Seurat LogNormalize (scale.factor=1000), `.layers['counts']` su sirovi count-ovi, a
`.obs['predicted.celltype']` je Azimuth anotacija urađena u R-u.

**Šta je zapravo `*_CoDi_KLD.csv` (provereno otvaranjem fajla):** per-ćelija tabela
(`cell_id, CoDi, CoDi_confidence, CoDi_dist, ...`) sa vrednostima tipa „CD4+ T cell" /
„Cytotoxic T cell" — to je **anotacija tipa ćelije iz nezavisnog metoda (CoDi)**, NIJE
gen-level divergencija. Dakle to je **drugi referentni cross-check za Stage 3 (anotaciju)**,
zajedno sa `predicted.celltype` — a NE referenca za DE (Stage 5). (Ime „KLD" zavarava.)

**6 zadataka (60 „internih" poena kvaliteta):** 1) QC & preprocessing (opravdati pragove)
· 2) Integracija + klasterovanje (batch correction po izboru, UMAP, klasteri) · 3) Anotacija
tipova ćelija (marker geni + Azimuth PBMC referenca) · 4) Composition analiza (udeli tipova
po uzorku) · 5) Diferencijalna ekspresija po tipu ćelije (svaki izloženi vs kontrola) +
pathway enrichment (GO/Reactome/KEGG) · 6) Size-specific efekti (unique 40nm / unique 200nm
/ shared / emergent u smeši — biološka interpretacija).

**4 deliverable-a (65 poena):** GitHub repo + environment fajl + README za reprodukciju
(20) · 3–5 dodatnih analiza/uvida (10) · PowerPoint slajdovi sa vizuelizacijama, na GitHub-u
(10) · video prezentacija 5–10 min na YouTube (20).

**Trenutno stanje koda** (provereno): `STAGE_REGISTRY` u `run_pipeline.py` ima žičano samo
Stage 1 (`qc`) i Stage 2 (`integration`) — i samo `src/qc.py` i `src/integration.py` postoje
kao čisti moduli. Faze 3–6 su zakomentarisan template. `data/processed/` i `results/` ipak
sadrže `03_annotated`, composition i DE izlaze — ali ta logika živi u `legacy/notebooks/`
(`03_annotation.ipynb`, `04_composition.ipynb`, `05_differential_expression.ipynb`), NE u
`src/`. Znači zadatak za faze 3–5 nije „vrati u registar" (nema šta u `src/` da se vrati) nego
**portuj logiku iz legacy notebook-a u čist `src/<stage>.py` modul i registruj ga**. Stage 6,
bonus analize i sva 4 deliverable-a (sem početnog README) ne postoje uopšte.

**Arhitektonska napomena (bitno za faze 4–6):** `StageSpec` pretpostavlja JEDAN `.h5ad`
ulaz/izlaz po fazi i `stage_status` proverava postojanje tog `.h5ad`. To lepo pasuje za 1–3
(svaka piše checkpoint), ali faze 4–6 proizvode **tabele/figure, ne novi AnnData** — i ne
treba da dupliraju ~2.7 GB h5ad bez razloga (disk!). Zato malo generalizuj registar: dozvoli
da `output_checkpoint` bude reprezentativni izlazni fajl (npr. glavni CSV) za status-proveru,
a faze 4–6 čitaju `03_annotated` (DAG, ne strogi lanac — 4 i 5 obe granaju iz anotacije,
6 čita DE tabele).

**Ključno ograničenje dizajna (boji faze 4–6):** 4 uzorka su od **jednog donora, jedan
uzorak po uslovu — nema bioloških replika.** Posledica: metode koje ocenjuju varijansu iz
replika (DESeq2/pseudobulk DE, propeller za sastav) ovde **nisu validne** u klasičnom smislu.
Za DE koristimo cell-level test (Wilcoxon / `rank_genes_groups`), a sastav prijavljujemo
deskriptivno (udeli + fold-change); svaki p-value je eksploratoran, uz eksplicitan caveat o
pseudoreplikaciji (ćelije nisu nezavisne jedinice). **Za DE nema eksternog gen-level
ground-truth-a** (CoDi_KLD je anotacija, ne DE referenca) — pa se DE validira *interno*:
konzistentnost preko tipova ćelija, biološka plauzibilnost top gena, i koherentnost pathway
enrichment-a. Budi iskren o ovome, ne izmišljaj eksternu potvrdu koje nema.

---

## Standardni kvalitetski okvir (drži ga u CLAUDE.md projekta)

Ovo su pravila koja je PSZ nosio implicitno; ovde ih zapiši jednom da svaki prompt ne mora
da ih ponavlja. **Prompt 0 ih ubacuje u `CLAUDE.md`.**

1. **Jedan inkrement po koraku.** Svaka faza se završava u pokretljivom stanju + checkpoint.
2. **Jedan code path.** Svaka nova faza je **jedan red** `StageSpec` u `STAGE_REGISTRY`;
   meni i `--stage` zovu istu `run_fn`. Bez paralelnih implementacija.
3. **Smoke-test paritet.** Svaka faza mora da prođe end-to-end na `--smoke-test`
   (~500 ćelija/uzorak) za sekunde i da upiše `*_smoke` checkpoint. Smoke i full su isti kod.
4. **`--debug` priča istinu o podacima.** Timing banneri, AnnData shape pre/posle, rich
   tabele (percentili, veličine klastera, broj po uzorku) i terminalni plotovi (plotext)
   u pravom terminalu. Pragovi se **opravdavaju odštampanim brojevima**, ne nagađaju.
5. **Reproduktivnost.** `RANDOM_SEED` svuda (numpy/scanpy/harmony/leiden/umap); headless
   `--stage` je kanonski put od kog zavisi ocenjivanje; checkpoint-ovi u `data/processed/`.
6. **Testovi namere (pytest), bez teških podataka/mreže.** Kao PSZ-ovih 122 testa — testiraj
   logiku (set-operacije, normalizaciju, pragove, parsiranje labela) na sitnim sintetičkim
   AnnData objektima. Brzi, deterministički, rade offline.
7. **Differential verification.** Gde POSTOJI nezavisna referenca, ukrsti se: **anotacija**
   (Stage 3) vs `obs['predicted.celltype']` (Azimuth/R) **i** vs `*_CoDi_KLD.csv` (CoDi metod) —
   troslojno slaganje. Za DE/sastav **nema** eksterne reference (jedan donor, bez replika),
   pa validacija je interna (konzistentnost + biološka plauzibilnost). Eksplicitno prijavi
   slaganje/neslaganje; ne izmišljaj referencu koje nema.
8. **Biološka interpretacija.** Svaki rezultat dobija figure/tabelu u `results/` + 2–4
   rečenice „šta ovo biološki znači" — ne samo brojeve.
9. **Jezik:** kod, komentari i izlaz na **engleskom** (per postavci); ovi promptovi i tvoja
   komunikacija sa mnom na srpskom.
10. **Verify-before-claim.** Pre nego što mi predaš artefakt ili tvrdnju o kodu/podacima,
    proveri svaku netrivijalnu činjenicu naspram **stvarnog fajla/komande** (otvori, `grep`,
    `head`), ne iz imena ili pretpostavke. Što ne možeš da potvrdiš — označi ASSUMED. Bar
    jedan ground-truth prolaz uvek ide PRE isporuke, ne posle mog podsećanja.

---

## Prompt 0 — Postavka, environment, skripte i orijentacija (sanity)

**Cilj:** reproduktivno okruženje, PSZ-stil skripte (`setup.ps1` / `run.ps1` / `Makefile`),
„doktor" provera, i potvrda da postojeće faze (1–2) rade — pre nego što diramo bilo šta novo.

```
Pročitaj README.md, run_pipeline.py i ceo src/ ovog projekta (GI/ANOA-ANO) da razumeš
postojeći staged pipeline. Zatim:

1. Napiši/ažuriraj CLAUDE.md u korenu projekta sa "Standardnim kvalitetskim okvirom"
   iz PROMPTOVI.md (10 pravila: jedan inkrement, jedan code path preko STAGE_REGISTRY,
   smoke-test paritet, --debug štampa prave brojeve, RANDOM_SEED svuda, testovi namere
   offline, differential verification vs predicted.celltype i CoDi_KLD, biološka
   interpretacija, kod na engleskom, verify-before-claim). Dodaj i ograničenje dizajna:
   jedan donor, bez bioloških replika -> DE je cell-level (Wilcoxon), sastav deskriptivan,
   p-value eksploratoran uz caveat.
2. Proveri da requirements.txt potpuno opisuje okruženje (scanpy, anndata, harmonypy,
   leidenalg, celltypist, gseapy/decoupler, plotext, rich, questionary, pytest, ...);
   ako nešto fali za kasnije faze, dodaj sa fiksiranim verzijama. (Provereno 2026-06-21:
   `scikit-misc` za seurat_v3 HVG je VEĆ tu; realno su falili `pytest`, `upsetplot`,
   `python-pptx` -- dodati.) Napomena: celltypist (model) i gseapy (Enrichr) traže internet
   na prvom pozivu -- to je u redu za same faze, ali testovi namere moraju ostati offline
   (sintetički AnnData).
3. Napravi PSZ-stil skripte za Windows, kao tanak omotač oko run_pipeline.py (NE drugi
   code path):
   - setup.ps1: pronađe Python 3.10+, napravi .venv, instalira requirements.txt, i
     proveri da su sva 4 raw .h5ad fajla u data/raw/ (ispiše šta fali i odakle se skida).
   - run.ps1: interaktivni meni + direktne komande koje samo prosleđuju run_pipeline.py:
     `run.ps1 qc|integration|annotation|composition|de|size` (-Smoke -Debug prekidači
     prosleđuju --smoke-test/--debug), `run.ps1 all` (sve faze redom), `run.ps1 test`
     (pytest), `run.ps1 check` (doctor), `run.ps1 menu`. Pošto samo forward-uje --stage,
     nove faze se pojavljuju automatski kad se registruju.
   - Makefile sa istim ciljevima za Linux/Mac (make setup/qc/all/test/check).
4. Dodaj lagani `--check` / "doctor" režim u run_pipeline.py (u stilu PSZ debug-a):
   ispiše Python verziju, da li su sva 4 raw .h5ad fajla prisutna, koji checkpoint-ovi
   postoje u data/processed/, i da li su ključne biblioteke importabilne — bez pokretanja
   teške analize.
5. Pokreni smoke-test za Stage 1 i Stage 2 (--smoke-test --debug) i potvrdi da prolaze i
   pišu *_smoke checkpoint-ove. Ako pukne, popravi pre nego što kreneš dalje.

Ne diraj logiku faza 1–2 osim ako smoke-test ne pukne. Ažuriraj README sa novim skriptama.
Na kraju mi javi tačno stanje svake faze (registrovana? full done? smoke done?) i šta ide
u sledeći korak.
```
**Definition of done:** `CLAUDE.md`, `setup.ps1`/`run.ps1`/`Makefile` i `--check` rade, smoke faza 1–2 zelen, jasan izveštaj stanja.

---

## Prompt 1 — Zadatak 1: QC & preprocessing (opravdaj pragove + figure)

**Cilj:** QC iz „starting point" pragova prebaciti u **opravdane** pragove, sa pre/posle
figurama i testovima. Postavka eksplicitno traži „Justify thresholds".

```
Radi na Stage 1 (qc). Postavka traži da opravdamo QC pragove, ne da ih nagađamo.

1. U --debug modu, qc.run mora da odštampa stvarnu raspodelu za ovaj dataset PRE filtriranja:
   percentili (1/5/25/50/75/95/99) za n_genes_by_counts, total_counts i pct_counts_mt,
   po uzorku i ukupno, kao rich tabela + plotext histogrami u pravom terminalu.
2. Na osnovu tih brojeva preispitaj QC_* pragove u config.py (min/max gena, max mito%).
   Ostavi ih konfigurabilnim, ali u komentaru config.py opravdaj svaki prag konkretnim
   brojem iz raspodele (npr. "6000 ~ p99, iznad su verovatni dubleti").
3. Sačuvaj QC figure u results/figures/: violin (genes/counts/mito) pre i posle filtriranja,
   i scatter total_counts vs pct_mt. I za full i za smoke (sufiks _smoke), isti kod.
4. U --debug ispiši "QC filter summary" tabelu: koliko ćelija/gena je palo na koji prag,
   po uzorku.
5. Dodaj pytest testove namere (offline, sitni sintetički AnnData): da filtriranje izbacuje
   tačno ćelije van pragova, da normalizacija daje target_sum po ćeliji, da log1p i HVG
   selekcija rade na counts sloju (ne na već-normalizovanom .X), i da Stage 1 resetuje .X
   na .layers['counts'] pre obrade.

Drži smoke-test paritet i RANDOM_SEED. Na kraju pokreni full Stage 1 i pokaži mi summary.
```
**Definition of done:** opravdani pragovi sa brojevima u config-u, QC figure (full+smoke), QC summary tabela, zeleni testovi.

---

## Prompt 2 — Zadatak 2: Integracija + klasterovanje (proveri batch-mixing)

**Cilj:** potvrditi da Harmony stvarno meša uzorke (a ne da klasteri = uzorci), UMAP figure,
sanity nad veličinama klastera.

```
Radi na Stage 2 (integration). Već radi PCA+Harmony+UMAP+Leiden; treba da DOKAŽEMO da
integracija radi i da je vizuelizujemo.

1. Sačuvaj u results/figures/ UMAP obojen po (a) uzorku i (b) Leiden klasteru, i to PRE i
   POSLE Harmony korekcije (pre = PCA-only UMAP, posle = Harmony UMAP), da se vidi efekat
   batch correction-a. Full i smoke.
2. U --debug dodaj kvantitativnu meru mešanja batch-eva (npr. LISI/iLISI ili prosti
   per-cluster sastav uzoraka kao rich tabela): za svaki klaster udeo svakog uzorka.
   Ako bi klaster bio ~100% jedan uzorak, to je crveni flag i mora da se istakne.
3. Iskoristi postojeće CLUSTER_TINY_FRACTION / CLUSTER_DOMINANT_FRACTION upozorenja:
   u --debug ispiši veličine klastera + warninge, plotext bar veličina klastera.
4. Pytest testovi namere (offline, sintetički): da je BATCH_KEY postavljen na svim ćelijama,
   da broj klastera nije degenerisan (1 ili = broj ćelija), da je seed determinisan
   (dva pokretanja istog smoke daju iste labele).

Smoke paritet + seed. Pokreni full Stage 2 i pokaži mi UMAP-ove i tabelu sastava klastera.
```
**Definition of done:** pre/posle UMAP figure, kvantitativna mera mešanja, sanity nad klasterima, zeleni testovi.

---

## Prompt 3 — Zadatak 3: Anotacija tipova ćelija (portuj iz legacy + dvostruki cross-check)

**Cilj:** portovati anotaciju iz `legacy/notebooks/03_annotation.ipynb` u čist `src/annotation.py`
+ `STAGE_REGISTRY`, anotirati tipove ćelija (marker geni + celltypist), i **ukrstiti sa OBE**
nezavisne reference — `predicted.celltype` (Azimuth) i CoDi (differential verification).

```
Logika anotacije postoji SAMO kao legacy/notebooks/03_annotation.ipynb (i izlaz
03_annotated.h5ad); NEMA je u src/ ni u STAGE_REGISTRY. Portuj je u čist modul.

1. Pročitaj legacy/notebooks/03_annotation.ipynb kao referencu, pa napiši src/annotation.py
   i dodaj run_annotation_stage + jedan red StageSpec u STAGE_REGISTRY (key="annotation",
   input_checkpoint="02_clustered", output="03_annotated") -- isti code path za --stage i meni.
   Smoke varijanta čita 02_clustered_smoke. Ne kopiraj naslepo iz notebook-a; preuzmi logiku,
   ali ispoštuj kvalitetski okvir (debug, smoke paritet, testovi).
2. Anotacija: per-klaster dodela imunog tipa (T CD4/CD8, B, NK, Mono CD14/FCGR3A, DC,
   Platelet...) kombinujući (a) celltypist (PBMC model) i (b) scoring MARKER_GENES iz
   config.py. Upiši finalnu labelu u obs i obrazloži po klasteru u --debug rich tabeli
   (top marker geni + score po klasteru).
3. DIFFERENTIAL VERIFICATION (obavezno) -- imamo DVE nezavisne reference, ukrsti sa obe:
   (a) obs['predicted.celltype'] (Azimuth pokrenut u R-u; postavka traži Azimuth PBMC
   referencu, a ona je već tu -- logična pretpostavka koju postavka dozvoljava), i
   (b) *_CoDi_KLD.csv kolone CoDi/CoDi_dist/CoDi_contrastive (per-ćelija anotacija iz CoDi
   metoda; merge po cell_id). Ispiši troslojnu agreement tabelu (naša vs Azimuth vs CoDi),
   % slaganja po tipu, i istakni ćelije/klastere gde se sve tri ne slažu -- to su biološki
   najzanimljiviji / najnesigurniji. Ne blokiraj se na pravom RDS-u (opciono, traži R most).
   Pazi na mapiranje imena tipova (naše "T cell (CD4)" vs Azimuth "CD4 T" vs CoDi "CD4+ T cell")
   -- napravi eksplicitan rečnik harmonizacije labela i pokrij ga testom.
4. Figure u results/figures/: UMAP obojen po tipu ćelije, i marker dotplot/matrixplot
   (geni x tipovi). Full i smoke.
5. Pytest testovi namere (offline): mapiranje klaster->tip je deterministički i pokriva sve
   klastere, da labela nije "Unknown" za klaster sa jasnim markerom, da scoring funkcija
   radi na sintetičkom AnnData.

Smoke paritet + seed. Pokreni full Stage 3, pokaži UMAP po tipu i troslojnu agreement tabelu
(naša vs Azimuth vs CoDi).
```
**Definition of done:** src/annotation.py portovan iz legacy + StageSpec u registru, troslojni agreement (naša vs Azimuth vs CoDi) sa rečnikom harmonizacije labela, UMAP+dotplot, zeleni testovi.

---

## Prompt 4 — Zadatak 4: Composition analiza (deskriptivni udeli + trend po veličini)

**Cilj:** udeli tipova ćelija po uzorku i relativno na kontrolu kao **primarni deskriptivni
rezultat** (bez replika nema validnog testa nad uslovima); eksploratorni test samo uz caveat.

```
Logika je u legacy/notebooks/04_composition.ipynb (izlazi u results/tables/), nije u src/.
Portuj je u src/composition.py i registruj StageSpec (input_checkpoint="03_annotated").
Ova faza proizvodi tabele/figure, ne novi h5ad -- za status-proveru koristi reprezentativni
CSV (npr. 04_composition_proportions.csv), ne prazni 2.7 GB checkpoint.

1. Po uzorku: broj i udeo svakog tipa ćelije; tabela udela i tabela relativno-na-kontrolu
   (log2 fold-change udela vs control). CSV-ovi u results/tables/, stacked bar + grouped bar
   u results/figures/. Full i smoke.
2. PRIMARNI rezultat su deskriptivni udeli + log2 fold-change vs control (jer je jedan
   uzorak po uslovu -- nema replika, pa propeller/ANOVA nad uslovima nije validan).
   Eksploratorno smeš dodati chi-square/Fisher na sirovim brojevima ćelija po tipu
   (izložen vs kontrola), ALI obavezno označi kao pseudoreplikaciju (ćelije nisu nezavisne)
   i ne prodaj ga kao čvrst dokaz. U --debug rich tabela: tip ćelije, udeo po uzorku,
   fold-change vs control, (eksploratorni) p/FDR sa caveat oznakom.
3. Biološka interpretacija (2-4 rečenice po izloženom uzorku): koji tipovi rastu/padaju i
   ima li doza/veličina-zavisnog trenda (40nm vs 200nm vs smeša).
4. Pytest testovi namere (offline): da udeli po uzorku sumiraju na 1, da relativno-na-kontrolu
   za kontrolu daje 0 (log2 1), da test vraća validne p-vrednosti na sintetičkim brojevima.

Smoke paritet + seed. Pokreni full Stage 4 i pokaži tabelu udela + fold-change vs control.
```
**Definition of done:** src/composition.py portovan iz legacy + StageSpec, udeli + relativno-na-kontrolu (primarno) + eksploratorni test sa caveat-om, figure/CSV, interpretacija, zeleni testovi.

---

## Prompt 5 — Zadatak 5: Diferencijalna ekspresija + pathway enrichment

**Cilj:** po tipu ćelije, svaki izloženi vs kontrola — DE (cell-level Wilcoxon, jer nema
replika) + volcano + GO/Reactome/KEGG enrichment, sa internom validacijom (nema eksterne
gen-level reference).

```
Logika je u legacy/notebooks/05_differential_expression.ipynb (delimični CSV izlazi);
portuj je u src/differential_expression.py i registruj StageSpec (input_checkpoint="03_annotated").
Tabele/figure, ne novi h5ad -- status preko reprezentativnog CSV-a.

1. Za svaki MAJOR tip ćelije i svaki izloženi uzorak vs control, izračunaj DE. VAŽNO:
   jedan donor, jedan uzorak po uslovu -> NEMA bioloških replika, pa pseudobulk+DESeq2/
   pydeseq2 NIJE validan (ocenjuje varijansu iz replika kojih nema). Default je cell-level
   Wilcoxon (scanpy rank_genes_groups, metoda 'wilcoxon') na log-normalizovanim podacima,
   sa eksplicitnim caveat-om o pseudoreplikaciji u komentaru i izlazu. Sačuvaj per-(tip,uzorak)
   CSV (gen, logFC, p, FDR) u results/tables/ (uskladi imena sa postojećim 05_DE_*).
2. Volcano plot po (tip ćelije, uzorak) u results/figures/, sa obeleženim top genima.
   Full i smoke (smoke sme da bude prazan/oskudan -- bitno je da kod prođe).
3. Pathway enrichment (gseapy/decoupler: GO, Reactome, KEGG) nad signifikantnim genima po
   grupi; tabela top puteva (naziv, NES/odds, FDR) + bar/dot figure.
4. VALIDACIJA (nema eksterne gen-level reference -- CoDi_KLD je anotacija, ne DE): validiraj
   INTERNO. (a) Konzistentnost: da li se isti stres/inflamacija geni javljaju preko više
   tipova ćelija; (b) biološka plauzibilnost top gena (poznati inflamatorni/oksidativni
   markeri?); (c) koherentnost: da li pathway enrichment iz koraka 3 podržava top gene.
   NE izmišljaj eksternu potvrdu. (Sanity nad ekspresijom konkretnih markera možeš ukrstiti
   sa CoDi anotacijom tipova iz Stage 3, ali to nije DE referenca.)
5. U --debug: broj sig. gena po (tip,uzorak) rich tabela, upozori gde je tip premali za
   pouzdan DE (prag min ćelija) i preskoči ga umesto da daješ smeće.
6. Pytest testovi namere (offline): Wilcoxon DE na sintetičkom AnnData izdvaja gen koji je
   namerno gore-regulisan u jednoj grupi, FDR korekcija monotona, filter "premalo ćelija"
   preskače tip umesto da pukne, enrichment poziv se gradi za neprazan gen-set.

Smoke paritet + seed. Pokreni full Stage 5 za 1-2 tipa ćelije i pokaži volcano + top pathway
tabelu + obrazloženje interne validacije.
```
**Definition of done:** de StageSpec (portovan iz legacy), cell-level Wilcoxon DE + volcano + enrichment, interna validacija (bez izmišljene reference), guard za male grupe, zeleni testovi.

---

## Prompt 6 — Zadatak 6: Size-specific efekti (set-logika + biologija)

**Cilj:** finalni biološki zaključak — šta je jedinstveno za 40nm, jedinstveno za 200nm,
deljeno, i šta se javlja **samo u smeši** (emergentno).

```
Dodaj Stage 6 (size_effects) u STAGE_REGISTRY. Ovo je glavni biološki zaključak projekta.
ULAZ su DE CSV tabele iz Stage 5 (ne h5ad checkpoint), IZLAZ su tabele/figure (ne novi h5ad)
-- status preko reprezentativnog CSV-a, kao faze 4–5. Ovo je čisto agregacija nad Stage 5
izlazom; ne učitavaj 2.7 GB AnnData ako ti ne treba.

1. Iz Stage 5 DE rezultata, po tipu ćelije, definiši set-pripadnost gena/puteva:
   - unique 40nm (sig u 40nm vs control, ne u 200nm)
   - unique 200nm
   - shared (sig u oba pojedinačna)
   - mixture-emergent (sig u smeši, a NIJE sig ni u 40nm ni u 200nm pojedinačno)
   Implementiraj ovu set-logiku eksplicitno i pokrij je testovima namere (ovo je tačno tip
   tihe logičke greške koju lako promašim -- testiraj sve 4 kategorije na sintetičkim setovima).
2. Figure: Venn ili UpSet (40nm/200nm/mixture) po tipu ćelije + tabela emergentnih gena/puteva.
   Full i smoke.
3. Biološka interpretacija (glavni nalaz projekta): da li male vs velike čestice izazivaju
   različit odgovor i da li smeša radi nešto što nijedna veličina sama ne radi -- 1 pasus,
   konkretno, vezano za gene/puteve i tipove ćelija.
4. CSV svih kategorija u results/tables/.

Smoke paritet + seed. Pokreni full Stage 6 i pokaži mi UpSet figuru + listu mixture-emergent
gena sa interpretacijom.
```
**Definition of done:** size_effects StageSpec, 4 kategorije sa testovima namere, Venn/UpSet + tabele, biološki zaključak.

---

## Prompt 7 — Deliverable: 3–5 dodatnih analiza (10 poena)

**Cilj:** predloži i implementiraj 3–5 dodatnih uvida (postavka eksplicitno traži).

```
Postavka traži 3-5 dodatnih analiza/uvida (10 poena). Prvo mi PREDLOŽI 5 kandidata sa
jednom rečenicom obrazloženja i procenom truda, pa kad izaberem implementiraj odabrane kao
zasebne stage-ove/skripte (isti registar/konvencija, smoke paritet, testovi gde ima logike).

Dobre kandidate razmotri:
- cell-cell komunikacija (CellPhoneDB/CellChat-stil ligand-receptor) izložen vs kontrola
- modul/score stresa i inflamacije (oksidativni stres, NF-kB, interferon) po tipu i uzorku
- dose-response: da li je odgovor na smešu aditivan (40nm+200nm) ili sinergičan/antagonistički
- pseudotime/aktivacioni gradijent u monocitima ili T ćelijama
- robustnost: stabilnost klastera/anotacije na variranje LEIDEN_RES i QC pragova
- globalna saglasnost anotacije sa CoDi referencom (concordance po tipu i uzorku, gde se
  metodi sistematski razilaze i da li to nosi biološki signal)

Svaka analiza: figura+tabela u results/ i 2-4 rečenice interpretacije.
```
**Definition of done:** 3–5 odabranih analiza implementirano, svaka sa figurom/tabelom i interpretacijom.

---

## Prompt 8 — Deliverable: repo + environment + README za reprodukciju (20 poena)

**Cilj:** najveći deliverable — sve mora da se reprodukuje jednom komandom, README kompletan.

```
Finalizuj reproduktivnost (20 poena -- najveći deliverable).

1. Prođi requirements.txt da je tačan i potpun za SVE faze (1-6 + bonus); fiksiraj verzije
   koje su zaista korišćene. Proveri da `pip install -r requirements.txt` u čistom venv-u
   povlači sve.
2. Prepiši README.md u kompletan: kratak biološki uvod, setup, kako se reprodukuju rezultati
   (full i smoke), objašnjenje SVAKE faze (1-6) u 2-3 rečenice, bonus analize, i
   "Deliverables checklist" sa statusom. Ažuriraj Status tabelu (trenutna je zastarela --
   tvrdi da su 3-6 "in progress" iako izlazi postoje).
3. Dodaj jedan "reproduce all" ulaz: skripta/Make/`run_pipeline.py --all` koja redom pokrene
   sve faze i regeneriše results/. Mora da prođe i u --smoke-test varijanti za sekunde
   (to je de-facto CI/end-to-end test celog pipeline-a).
4. Pokreni ceo pytest paket i sve smoke faze; sve mora biti zeleno. Izlistaj results/ manifest
   (koje figure/tabele postoje) na kraju README-a.
5. Commit-uj na GitHub sa urednom istorijom.

Pokaži mi finalni README i izlaz "reproduce all --smoke-test".
```
**Definition of done:** čist venv install radi, README kompletan + tačan status, "reproduce all" prolazi (full+smoke), testovi zeleni, push.

---

## Prompt 9 — Deliverable: PowerPoint slajdovi (10 poena)

**Cilj:** prezentacija svih rezultata sa vizuelizacijama, sačuvana u repo-u.

```
Napravi PowerPoint (.pptx) sa svim rezultatima i vizuelizacijama (10 poena), sačuvaj ga u
repo (npr. results/slides/ ili docs/). Koristi figure koje pipeline već generiše u
results/figures/. Struktura:
1. Naslov + pitanje (veličina nanoplastike i imuni odgovor) i dataset (4 uzorka).
2. Metod: pipeline faze 1-6 kao jedan dijagram toka.
3. Po jedan-dva slajda po zadatku: QC, integracija/UMAP, anotacija (+ agreement vs Azimuth),
   composition, DE+pathway, size-specific (UpSet + glavni nalaz).
4. Bonus analize.
5. Zaključak: glavni biološki nalaz + ograničenja.
Tekst na engleskom, sažeto (bullet-i, ne pasusi). Embeduj prave PNG figure, ne placeholder-e.
```
**Definition of done:** `.pptx` u repo-u sa pravim figurama, pokriva svih 6 zadataka + bonus + zaključak.

---

## Prompt 10 — Deliverable: video skripta (20 poena)

**Cilj:** narativna skripta za 5–10 min video prezentaciju (snimanje radiš ti).

```
Napiši skriptu za video prezentaciju od 5-10 minuta (deliverable nosi 20 poena) na osnovu
slajdova i rezultata. Po sekcijama, sa procenom trajanja svake (ukupno 5-10 min):
- hook/pitanje (~30s), dataset i dizajn (~45s), pipeline pregled (~1min),
- ključni rezultati po zadatku sa konkretnim brojevima iz naših tabela/figura,
- glavni biološki nalaz (size-specific + mixture-emergent),
- ograničenja i šta bi sledeće (~30s).
Za svaku sekciju navedi koji slajd/figura je na ekranu. Ton: jasan, tehnički ali pristupačan.
Engleski (ili srpski ako tako odlučim). Daj mi i 2-3 rečenice opisa za YouTube + naslov.
```
**Definition of done:** tajmovana skripta vezana za slajdove + YouTube naslov/opis; spremno za snimanje.

---

### Kako koristiti

Idi redom: **0 → 10**. Posle svakog prompta proveri „Definition of done", pogledaj
figure/tabele i tek onda kreni dalje. Faze 0–2 su brze (postojeći kod), 3–6 su jezgro
analize, 7–10 su deliverable-i. Ako neki korak otvori novi pod-problem, tretiraj ga kao
mini-prompt u istom duhu (jedan inkrement, verifikacija, testovi) — kao što je PSZ dobio
watchdog, throttling i differential audit kao zasebne korake.
