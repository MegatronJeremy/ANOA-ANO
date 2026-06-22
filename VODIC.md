# Vodič kroz projekat — šta je ovo i kako da ga tumačim

Ovaj dokument objašnjava ceo projekat „od nule", jezikom koji ne pretpostavlja da si pratio
kurs. Cilj: da razumeš šta je urađeno, šta svaki rezultat znači, i da možeš da odbraniš projekat.
Pojmovi su objašnjeni u **[Rečniku](#2-rečnik-pojmova)**, a na kraju je **[priprema za
odbranu](#8-pitanja-za-odbranu-i-odgovori)**.

---

## Sadržaj

1. [Velika slika — šta i zašto](#1-velika-slika--šta-i-zašto)
2. [Rečnik pojmova](#2-rečnik-pojmova)
3. [Kako se uklapa u gradivo kursa](#3-kako-se-uklapa-u-gradivo-kursa)
4. [Pipeline korak po korak (6 zadataka)](#4-pipeline-korak-po-korak)
5. [Bonus analize](#5-bonus-analize)
6. [Glavni nalazi — priča projekta](#6-glavni-nalazi--priča-projekta)
7. [Da li imam sve rezultate?](#7-da-li-imam-sve-rezultate)
8. [Pitanja za odbranu i odgovori](#8-pitanja-za-odbranu-i-odgovori)
9. [Kako da pokrenem / pokažem](#9-kako-da-pokrenem--pokažem)

---

## 1. Velika slika — šta i zašto

**Biološko pitanje.** Nanoplastika (sićušne plastične čestice) dospeva u krv i tamo se sreće sa
imunim ćelijama. Pitanje je: **da li veličina čestice menja imuni odgovor?** Konkretno —
reaguju li imune ćelije drugačije na male (40 nm) vs velike (200 nm) čestice, i da li smeša
obe veličine radi nešto što nijedna sama ne radi?

**Kako merimo odgovor.** Koristimo **scRNA-seq** (single-cell RNA sequencing) — tehnologiju
koja za **svaku pojedinačnu ćeliju** izmeri koliko je koji gen aktivan (koliko RNK pravi).
Pošto su geni „uključeni/isključeni" različito u različitim stanjima ćelije, promena u
ekspresiji gena = otisak imunog odgovora.

**Podaci.** 4 uzorka krvnih imunih ćelija (PBMC) jednog donora:

| Uzorak | Izloženost |
|---|---|
| `PSNP_40nm` | 40 nm čestice |
| `PSNP_200nm` | 200 nm čestice |
| `PSNP_mixture` | 40 + 200 nm smeša |
| `control` | bez izloženosti (referenca) |

Svaki uzorak ima ~8.000 ćelija (ukupno ~34.000), i za svaku ćeliju ~20.000 gena. To je velika
tabela brojeva (ćelije × geni) koju pipeline pretvara u biološki zaključak.

> **Zašto „single-cell" a ne prosečno?** Da si izmerio prosek preko svih ćelija (to se zove
> *bulk RNA-seq*), promena u monocitima bi se „izgubila" ako T ćelije rade suprotno. Single-cell
> razdvaja ćelije po tipu, pa vidiš ko tačno reaguje. (Lekcija 4 kursa upravo o tome.)

---

## 2. Rečnik pojmova

Pročitaj ovo prvo — sve ostalo se oslanja na ove pojmove.

- **Gen / ekspresija gena.** Gen je „recept" za protein. *Ekspresija* = koliko se taj recept
  trenutno koristi, mereno količinom RNK. Visoka ekspresija = gen je „upaljen".
- **RNK-seq (RNA-seq).** Tehnologija koja prebroji RNK molekule → koliko je koji gen aktivan.
- **scRNA-seq.** RNA-seq po pojedinačnoj ćeliji (a ne prosek preko miliona ćelija).
- **Ćelija × gen matrica.** Glavna tabela: svaki red = ćelija, svaka kolona = gen, vrednost =
  koliko je tog gena izmereno (broj „count-ova"). To je `.X` u podacima.
- **AnnData / `.h5ad`.** Format fajla koji drži tu matricu (`.X`) + opis ćelija (`.obs`, npr.
  iz kog uzorka je ćelija, koji joj je tip) + opis gena (`.var`). (Lekcija 4/5 kursa.)
- **QC (quality control).** Izbacivanje loših ćelija (npr. mrtve, prazne kapljice) pre analize.
- **Normalizacija.** Izjednačavanje ćelija da bi bile uporedive (neke su uhvaćene „dublje" od
  drugih). Posle se radi `log` da veliki brojevi ne dominiraju.
- **HVG (highly variable genes).** ~2000 gena koji se najviše razlikuju među ćelijama —
  na njima se gradi analiza (ostali su uglavnom „šum/pozadina").
- **Batch effect / integracija.** Tehničke razlike između uzoraka koje nisu biologija. Ako se
  ne isprave, ćelije se grupišu po *uzorku* umesto po *tipu*. **Harmony** je metod koji to
  ispravlja (integracija).
- **PCA.** Sažimanje ~2000 dimenzija (gena) u ~30 najinformativnijih — da algoritmi rade brže
  i čistije. (Linearna redukcija dimenzija; Lekcija 4.)
- **UMAP.** 2D „mapa" svih ćelija gde su slične ćelije blizu. To su one šarene tačkaste slike.
  Služi samo za *vizuelizaciju* (ne za odluke).
- **Klaster (Leiden).** Algoritam grupiše ćelije u grupe sličnih (klastere). Svaki klaster ≈
  jedan tip ćelije.
- **Anotacija / tip ćelije / lineage.** Davanje imena klasterima: T ćelija, B ćelija, NK,
  Monocit... *Lineage* = krupna grupa (mi koristimo 5: T, B, NK, Monocyte, Other).
- **Marker gen.** Gen karakterističan za tip ćelije (npr. `MS4A1` za B ćelije, `CD14` za
  monocite). Ako je upaljen → ta ćelija je verovatno taj tip.
- **DE (differential expression).** Poređenje: koji geni se menjaju kod izloženog uzorka u
  odnosu na kontrolu, **unutar istog tipa ćelije**.
- **log2 fold-change (log2FC).** Koliko se nešto promenilo, u „duplo" jedinicama. +1 = duplo
  više, −1 = duplo manje, 0 = bez promene. +1.35 ≈ 2.5× više.
- **p-value / FDR.** Verovatnoća da je razlika slučajna. Mali p (npr. < 0.05) = verovatno
  stvarno. **FDR** je p ispravljen jer testiramo hiljade gena odjednom.
- **Pathway / enrichment.** Geni ne rade sami — pripadaju „putevima" (npr. *inflamacija*,
  *oksidativni stres*). Enrichment kaže: lista promenjenih gena pripada kom putu → biološko
  značenje umesto gole liste imena. (Baze: GO, KEGG, Reactome.)
- **Pseudoreplikacija (važno za ovaj projekat).** Tretiranje ćelija jedne osobe kao da su
  nezavisni uzorci. Pošto imamo **jednog donora**, statistika je tu *eksploratorna*, ne tvrd
  dokaz — i mi to svuda otvoreno kažemo.

---

## 3. Kako se uklapa u gradivo kursa

Lekcija 4 kursa („Single cell RNA sequencing") na poslednjem slajdu prikazuje **„Single-cell
RNA downstream analysis workflow"** — i to je **tačno** ono što naš pipeline radi, korak po korak:

| Korak iz kursa (Lekcija 4) | Naša faza |
|---|---|
| QC | Stage 1 |
| Normalisation | Stage 1 |
| Feature selection (HVG) | Stage 1 |
| Data correction (batch) | Stage 2 (Harmony) |
| Dimensionality reduction & visualisation (PCA, UMAP) | Stage 2 |
| Clustering | Stage 2 (Leiden) |
| Annotation | Stage 3 |
| Compositional analysis (Enriched vs Control) | Stage 4 |
| Differential expression | Stage 5 |
| (Gene dynamics / trajectory — opciono) | dotaknuto u bonusu |

Dakle, projekat je „udžbenički" scRNA-seq tok primenjen na konkretno pitanje (nanoplastika).
Takođe, *zlatno pravilo bioinformatike* iz Lekcije 1 — **„Never trust your tools or data"** —
je razlog zašto svuda radimo **proveru naspram nezavisne reference** (vidi Stage 3) i imamo
**37 automatskih testova**.

---

## 4. Pipeline korak po korak

Svaka faza: **šta radi → šta je rezultat → gde je fajl → kako da ga tumačiš.** Sve se pokreće
istom komandom (`python run_pipeline.py --stage <ime>` ili meni `.\run.ps1`).

### Stage 1 — QC i preprocessing
- **Šta radi:** izbaci loše ćelije (premalo gena, previše mitohondrijske RNK = umiruća ćelija),
  normalizuje, i izabere ~2000 HVG. Pragovi su *opravdani stvarnim brojevima* (npr. gornja
  granica gena ≈ 99. percentil — iznad su verovatni dubleti).
- **Rezultat / fajl:** `results/figures/01_qc_violin_before/after.png` (raspodele pre/posle),
  `01_qc_scatter_counts_mito.png`.
- **Kako tumačiš:** „posle" violina treba da bude „čistija" od „pre" — odsečeni su repovi
  (ekstremi). To znači da u dalju analizu ulaze samo zdrave ćelije.

### Stage 2 — Integracija i klasterovanje
- **Šta radi:** PCA → Harmony (ispravi batch) → UMAP (mapa) → Leiden (klasteri).
- **Rezultat / fajl:** `02_umap_pre_harmony_by_sample.png` vs `02_umap_by_sample.png`
  (pre/posle), `02_umap_by_cluster.png`, tabela `02_cluster_composition.csv`.
- **Kako tumačiš:** Na slici **„pre"** ćelije su razdvojene po boji (uzorku) — to je tehnička
  greška. Na **„posle"** boje su izmešane — sad se ćelije grupišu po *biologiji*, ne po uzorku.
  To je dokaz da integracija radi.

### Stage 3 — Anotacija tipova ćelija
- **Šta radi:** alat `celltypist` + marker geni daju ime svakom klasteru (T/B/NK/Monocyte/...).
- **Provera (ključno!):** poredimo našu anotaciju sa **dve nezavisne reference** koje su već u
  podacima — *Azimuth* (urađen u R-u) i *CoDi*. Ako se sve tri slažu, verujemo imenima.
- **Rezultat / fajl:** `03_umap_lineage.png` (mapa obojena po tipu), `03_marker_dotplot.png`,
  `results/tables/03_annotation_agreement.csv`.
- **Brojevi (pun dataset):** slaganje naše anotacije sa Azimuth/CoDi: **T ćelije ~100%**,
  **B ćelije 90–95%**, **Monociti 84–93%**. **NK ćelije se slabo slažu** (0.4% / 11%) —
  iskreno navodimo: NK i citotoksične T ćelije su transkripciono jako slične, pa je to poznata
  „teška granica"; različiti alati ih različito zovu. To NE kvari glavni zaključak (koji se
  oslanja na monocite i T ćelije, gde je slaganje visoko).

### Stage 4 — Composition analiza
- **Šta radi:** koliki je *udeo* svakog tipa ćelije u svakom uzorku, i kako se menja u odnosu
  na kontrolu (log2 fold-change).
- **Rezultat / fajl:** `04_composition_stacked.png`, `04_composition_grouped.png`,
  `results/tables/04_composition_proportions.csv`, `..._relative_to_control.csv`.
- **Ključni broj:** udeo **monocita** skoči sa **3.3% (kontrola)** na **8.5% (200 nm)** — više
  nego dupliranje (log2FC ≈ +1.35). Znači velike čestice privlače/aktiviraju monocite.

### Stage 5 — Diferencijalna ekspresija + pathway
- **Šta radi:** za svaki tip ćelije, koji geni se menjaju kod svakog izloženog uzorka vs
  kontrola (Wilcoxon test), pa **enrichment** (kom biološkom putu pripadaju).
- **Rezultat / fajl:** `05_volcano_<tip>_<uzorak>.png` (volcano grafici), po jedna CSV tabela
  za svaku kombinaciju `05_DE_<tip>_<uzorak>_vs_control.csv`, `05_pathway_enrichment.csv`.
- **Kako tumačiš volcano:** svaka tačka = gen. Desno = pojačan, levo = utišan; više gore =
  statistički sigurnije. Crvene tačke = značajni geni.
- **Brojevi:** **monociti reaguju najjače** (864 značajna gena na 200 nm), kod limfocita je
  red veličine **200 nm > 40 nm > smeša**.

### Stage 6 — Size-specific efekti (glavni zaključak)
- **Šta radi:** po tipu ćelije razvrsta gene u 4 grupe: jedinstveni za 40 nm, jedinstveni za
  200 nm, deljeni (oba), i **mixture-emergent** (značajni *samo u smeši*, ni u jednoj veličini
  posebno).
- **Rezultat / fajl:** `06_size_categories.png`, `results/tables/06_size_specific_summary.csv`
  i `06_size_specific_genes.csv`.
- **Brojevi (monociti):** unique_40nm **174**, unique_200nm **340**, shared **524**,
  **mixture_emergent 180**. U svakom tipu je `unique_200nm > unique_40nm` → veće čestice prave
  više jedinstvenih promena; a smeša pravi *emergentan* odgovor najjače u monocitima.

---

## 5. Bonus analize

Pet dodatnih analiza (deliverable traži 3–5), svaka sa figurom i tabelom (`results/.../07_*`):

1. **Module scoring** — koliko je upaljen program *oksidativnog stresa / inflamacije /
   interferona* po uzorku i tipu. **Nalaz:** monociti imaju jasno povišen stres+inflamaciju
   (skorovi ~1.3–1.5 vs ~0.5 ili negativno kod ostalih). Potvrđuje DE priču.
2. **Mixture aditivnost** — da li je smeša = 40 nm + 200 nm? **Nalaz:** 84–93% gena je
   *sub-aditivno* (smeša slabija od zbira) → odgovor na smešu nije prosto sabiranje.
3. **Robusnost klasterovanja** — da li se zaključci menjaju ako promenimo parametar rezolucije?
   **Nalaz:** klasteri stabilni (ARI 0.76–0.96) → rezultati nisu artefakt proizvoljnog izbora.
4. **Ligand-receptor** — gruba mera ćelijske komunikacije izložen vs kontrola.
5. **Dose-response** — magnituda poremećaja (broj značajnih gena) vs veličina čestice.

---

## 6. Glavni nalazi — priča projekta

Tri nezavisne analize pričaju **istu priču**, što joj daje težinu:

1. **Monociti su glavni „reaktor".** Najviše značajnih gena (Stage 5), najveći skok udela
   (Stage 4: 3.3%→8.5% na 200 nm), i najviši skorovi stresa/inflamacije (bonus 1).
2. **Veće čestice (200 nm) izazivaju jači i specifičniji odgovor** od malih (40 nm) — više
   jedinstvenih gena u svakom tipu ćelije (Stage 6, dose-response).
3. **Smeša nije prost zbir.** Kod limfocita je *slabija* od pojedinačnih veličina (sub-aditivna,
   84–93%), ali u **monocitima pravi emergentan odgovor** (180 gena značajno samo u smeši) —
   nešto što nijedna veličina sama ne izaziva.

**Jednom rečenicom:** *veličina nanoplastike menja imuni odgovor, monociti reaguju najjače, a
kombinacija veličina deluje drugačije nego svaka veličina pojedinačno.*

---

## 7. Da li imam sve rezultate?

Da. Sve je generisano na **punim podacima** (33.232 ćelije) i na GitHub-u:

- ✅ **Kod + README + reproducibilnost** (20p): ceo repo, `requirements.txt`, jedna komanda
  `python run_pipeline.py --all` regeneriše sve. 37 automatskih testova prolazi.
- ✅ **Svih 6 zadataka** (figure + tabele u `results/figures/` i `results/tables/`).
- ✅ **5 bonus analiza** (10p) — `results/.../07_*`.
- ✅ **PowerPoint** (10p) — `results/slides/GI_nanoplastic.pptx` (11 slajdova).
- ✅ **Video skripta** (20p) — `VIDEO_SCRIPT.md` (tekst spreman; **snimanje je na tebi**).

**Jedino što ostaje ručno tebi:** (1) snimi video po skripti i okači na YouTube; (2) prijavi/
pošalji projekat na mejl iz postavke (`vladimir.kovacevic@etf.bg.ac.rs`).

Brza provera da je sve na mestu: `.\run.ps1 check` (ispiše da li su podaci, checkpoint-i i
biblioteke spremni).

---

## 8. Pitanja za odbranu i odgovori

Verovatna pitanja i kratki, sigurni odgovori:

- **„Zašto single-cell a ne bulk RNA-seq?"** — Jer prosek preko svih ćelija sakriva ko reaguje;
  monocitni signal bi se izgubio u proseku sa T ćelijama. Single-cell razdvaja po tipu.
- **„Zašto Harmony / čemu integracija?"** — Bez nje se ćelije grupišu po uzorku (tehnička
  greška), a ne po biologiji. Pokazujem to pre/posle UMAP-om.
- **„Kako znaš da su tipovi ćelija tačni?"** — Ne verujem jednom alatu: poredim sa dve nezavisne
  reference (Azimuth, CoDi). Slaganje ~90–100% za T/B/Monocyte. (Zlatno pravilo: ne veruj alatu.)
- **„Zašto NK ćelije slabo se slažu?"** — NK i citotoksične T su transkripciono slične; poznata
  teška granica, alati ih različito zovu. Glavni zaključak se na njih ne oslanja.
- **„Zašto Wilcoxon a ne DESeq2/pseudobulk?"** — Jer imamo **jednog donora, bez bioloških
  replika**. Metode sa replikama (DESeq2, propeller) ovde nisu validne; cell-level Wilcoxon je
  ispravan izbor, a p-vrednosti tretiram kao *eksploratorne* uz caveat o pseudoreplikaciji.
- **„Koliko su rezultati pouzdani?"** — Tri nezavisne analize daju istu priču (monociti, 200 nm,
  emergentna smeša); klasteri robustni (ARI 0.76–0.96); ali bez replika je zaključak
  *eksploratoran*, ne klinički dokaz.
- **„Šta je log2 fold-change +1.35?"** — ~2.5× više nego u kontroli.
- **„Šta je glavni nalaz?"** — Veličina čestice menja imuni odgovor; monociti reaguju najjače;
  smeša nije prost zbir (emergentno u monocitima, sub-aditivno u limfocitima).

---

## 9. Kako da pokrenem / pokažem

```powershell
.\run.ps1 check        # provera okruzenja (sta je spremno)
.\run.ps1              # interaktivni meni (izaberi fazu strelicama)
.\run.ps1 all          # pokreni ceo pipeline redom (pun, ~15-20 min)
.\run.ps1 all -Smoke   # brza verzija na uzorku (~2-3 min) - za demo
.\run.ps1 test         # 37 automatskih testova
.\run.ps1 slides       # ponovo izgradi PowerPoint iz figura
```

Figure su u `results/figures/`, tabele u `results/tables/`, slajdovi u `results/slides/`.
Tehnički detalji svake faze su u `README.md`; razvojni tok je u `PROMPTOVI.md`.
