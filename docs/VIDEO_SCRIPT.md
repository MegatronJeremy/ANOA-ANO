# Video presentation script (5-10 min)

Narration for the results video. Slide numbers match the current 16-slide deck
(`results/slides/GI_nanoplastic.pptx`). Aim for about 6 to 7 minutes when read
at a steady pace. Written to be spoken, not read off the screen.

---

### 0. Intro  ·  ~20s  ·  [Slide 1: title]

"Hi, my name is Vuk Djordjević. I'm a master's student at the School of Electrical
Engineering in Belgrade, and this is my project for the Genome Informatics course,
2026. The project is a single-cell analysis of how human immune cells respond to
nanoplastic particles of different sizes."

### 1. The question  ·  ~35s  ·  [Slide 1: title]

"Nanoplastics are tiny plastic particles, and they've now been found in human blood.
Once they're in the blood, they come into direct contact with our immune cells. The
question I looked at is whether particle size matters. Do small and large nanoplastics
cause different immune responses, and does a mixture of both do something that neither
size does on its own? To answer this I used single-cell RNA sequencing, which lets me
read which genes are active in each individual cell."

### 2. Dataset and design  ·  ~45s  ·  [Slide 2: question and dataset]

"The data is four samples from a single donor. Peripheral blood immune cells were
exposed to 40-nanometre particles, to 200-nanometre particles, to a mixture of both,
and one sample was left unexposed as a control. That's about 34,000 cells in total.
One thing I want to be upfront about: there is only one donor and one sample per
condition, so there are no biological replicates. That limitation shaped how I did the
statistics, and I flag it throughout the project."

### 3. Pipeline overview  ·  ~45s  ·  [Slide 3: pipeline]

"The whole analysis runs through one pipeline of six stages: quality control,
integration, cell-type annotation, composition analysis, differential expression, and
size-specific effects. A single command reproduces all of it from the raw data, and
each stage has automated tests."

### 4. Quality control  ·  ~40s  ·  [Slides 4-5: QC violins, then counts-vs-mito scatter]

"The first stage is quality control, where I remove low-quality cells before trusting
anything. This slide shows three violin plots, one per quality metric: the number of
genes detected per cell, the total reads per cell, and the percent of mitochondrial
reads. The height of each violin is the value; its width is how many cells sit at that
value. Instead of guessing thresholds, I set them from the real distribution: the upper
gene cutoff near the 99th percentile removes likely doublets, meaning two cells read as
one, and the mitochondrial cutoff between the 95th and 99th percentile removes dying
cells. The next slide is a scatter of the same filtering: every dot is a cell, the
x-axis is total reads, the y-axis is percent mitochondrial, and the lines are the
thresholds drawn over the actual data."

### 5. Integration  ·  ~45s  ·  [Slide 6: before and after UMAP]

"The second stage is integration. Both panels are UMAP plots, which are a two-dimensional
map where cells with similar gene expression land near each other. The axes, UMAP-1 and
UMAP-2, have no units; only relative position matters. The colour is which sample a cell
came from. On the left, before correction, the cells separate by sample, and that's a
technical effect, not real biology. On the right, after Harmony correction, the colours
mix together and the cells group by cell type instead. One honest caveat: because each
condition is a single sample, the batch effect and the treatment effect are tied
together, so for the actual gene-level results later I work on the uncorrected data,
where this doesn't bias anything."

### 6. Annotation  ·  ~45s  ·  [Slide 7: lineage UMAP and dotplot]

"The third stage is labelling the cell types. On the left is the same UMAP, now coloured
by cell-type lineage, so each coloured island is T cells, B cells, NK cells, or
monocytes. On the right is a marker dot-plot: the rows are cell types, the columns are
known marker genes, the colour of each dot is how strongly that gene is expressed, and
the size of the dot is the fraction of cells expressing it, so a big dark dot means that
gene is a clean marker for that type. I assigned the types with a tool called celltypist
plus these markers, then cross-checked against two independent methods already in the
data, Azimuth and CoDi. My labels agreed with both at about 93 percent, so I'm confident
in the identities before comparing conditions."

### 7. Composition  ·  ~40s  ·  [Slide 8: composition bars]

"The fourth stage compares how the proportions of each cell type change across
conditions. Both charts have proportion of cells on the y-axis, as a fraction of the
sample. The left is stacked bars, one bar per sample split into coloured cell-type bands;
the right is the same data grouped by lineage so you can compare a type across samples.
The clearest change is monocytes under 200-nanometre exposure: roughly two and a half
times more than the control, a log2 fold-change of about plus one-point-three-five, while
the other cell types barely move. Since there are no replicates, I report these as
descriptive proportions and fold-changes and treat any statistical test as exploratory."

### 8. Differential expression and the main result  ·  ~75s  ·  [Slide 9: dose-response and size categories], [Slide 10: key finding]

"The fifth stage asks a deeper question: within each cell type, which genes actually
change on exposure? This slide has two charts. On the left, the dose-response plot: the
x-axis is particle size, the y-axis is the number of significant genes that change.
Monocytes stack the highest, because they're the immune cells that eat and clear foreign
particles, so it makes sense they react the strongest, and the 200-nanometre bar is
higher than the 40-nanometre one: bigger particle, bigger response. On the right, the
size-categories chart: grouped bars per cell type, y-axis is again the number of changed
genes, split into four groups, genes unique to 40 nanometres, unique to 200 nanometres,
shared by both, and mixture-emergent, meaning significant only in the mixture.

And that fourth group is the main finding. In monocytes, about 180 genes respond only
when both sizes are present, and to neither size on its own. In the lymphocytes the
mixture is actually weaker than either single size. So the combination of the two sizes
does something the individual sizes do not."

### 9. What the mixture does biologically  ·  ~55s  ·  [Slide 11: emergent inflammation]

"So what is that new thing the mixture does? I ran pathway enrichment on those 180
mixture-only monocyte genes, and the result points clearly at inflammation. The top
pathways, with very strong statistical support, are interleukin, cytokine, and TNF
signalling, along with the response to lipopolysaccharide. That last one is the program
a monocyte runs during a bacterial infection, except here it's set off by plastic. Two
of the strongest genes that go up in the mixture are RIPK2 and TRAF1, which are core
innate-immune and TNF-signalling genes. My reading is that 40 and 200 nanometres each
cause small changes that stay below a threshold, but together they push the monocyte
past that threshold into an inflammatory state. Since real-world exposure is always to a
mix of particle sizes, testing single sizes alone could underestimate the risk."

### 10. Additional analyses  ·  ~50s  ·  [Slide 12: list], [Slides 13-15: the five figures]

"For the additional-insights part of the project I implemented five extra analyses,
which this slide lists, and the next slides show each one. The first is module scoring,
shown as a heatmap where the rows are samples, the columns are gene programs like stress
and inflammation, and the colour is the mean score for that program in that sample; the
inflammation program lights up under exposure. The second is a mixture-additivity test,
one scatter per cell type where the x-axis is the response you'd expect if you simply
added 40 and 200 nanometres, and the y-axis is the response actually observed in the
mixture; the diagonal is 'mixture equals sum of parts', and the monocyte points sit
above that line, meaning the mixture does more than the sum. The third is a clustering
robustness check: the x-axis is the clustering resolution, and the green line is the
agreement with the baseline clustering, which stays high, so the cell groups aren't an
artefact of one setting. The fourth is a ligand-receptor analysis, where the x-axis is
the log2 fold-change versus control for signalling pairs, showing how cell-to-cell
communication shifts on exposure. And the fifth, dose-response, I
already showed earlier on the differential-expression slide. All five point the same
way: a size-dependent, non-additive response centred on monocytes."

### 11. Limitations and close  ·  ~30s  ·  [Slide 16: limitations]

"To be clear about the limits: one donor and no replicates mean the differential
expression is done per cell and is exploratory, with a pseudoreplication caveat, and
there's no external gene-level reference, so I validated the results internally and
biologically. Everything reproduces from the repository with a single command and a full
test suite. Thanks for watching."

---

**YouTube title:** Nanoplastic size and mixtures change the immune response: a single-cell study

**YouTube description (2-3 sentences):**
A single-cell RNA-seq analysis of human immune cells exposed to 40 nm, 200 nm, and mixed
polystyrene nanoplastics versus an unexposed control. Monocytes respond the strongest,
larger particles change more cell-type-specific genes, and the 40+200 nm mixture drives
an inflammation-heavy monocyte response that neither size produces alone. Code, figures,
and a fully reproducible pipeline are on GitHub.
