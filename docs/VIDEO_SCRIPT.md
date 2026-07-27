# Video presentation script (5-10 min)

Narration for the results video. Slide numbers match the current 15-slide deck
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

### 4. Quality control  ·  ~35s  ·  [Slides 4-5: QC violins, then counts-vs-mito scatter]

"The first stage is quality control, where I remove low-quality cells before trusting
anything. Instead of picking thresholds by guessing, I looked at the real distribution
of all 34,000 cells and set the cutoffs from that. The upper gene cutoff sits near the
99th percentile, which removes likely doublets, meaning two cells read as one. The
mitochondrial cutoff sits between the 95th and 99th percentile, which removes dying
cells. The scatter plot on the next slide shows those thresholds drawn over the actual
data."

### 5. Integration  ·  ~45s  ·  [Slide 6: before and after UMAP]

"The second stage is integration. On the left is the data before batch correction. You
can see the cells separate by which sample they came from, and that's a technical
effect, not real biology. On the right is after Harmony correction. Now the samples mix
together and the cells group by cell type instead. I should be honest about one thing
here: because each condition is a single sample, the batch effect and the treatment
effect are tied together. So for the actual gene-level results later, I work on the
uncorrected data, where this doesn't bias anything."

### 6. Annotation  ·  ~40s  ·  [Slide 7: lineage UMAP and dotplot]

"The third stage is labelling the cell types. I used a tool called celltypist together
with known marker genes to assign types like T cells, B cells, NK cells, and monocytes.
To check the labels, I compared them against two other methods that were already in the
data: an Azimuth annotation and the CoDi method. My labels agreed with both at about 93
percent, so I'm confident in the cell identities before comparing conditions."

### 7. Composition  ·  ~40s  ·  [Slide 8: composition bars]

"The fourth stage compares how the proportions of each cell type change across
conditions. The clearest change is in the monocytes under 200-nanometre exposure. There
are roughly two and a half times more monocytes than in the control, while the other
cell types barely move. Since there are no replicates, I report these as descriptive
proportions and fold-changes, and I treat any statistical test here as exploratory."

### 8. Differential expression and the main result  ·  ~70s  ·  [Slide 9: dose-response and size categories], [Slide 10: key finding]

"The fifth stage asks a deeper question: within each cell type, which genes actually
change on exposure? Two patterns stand out. First, monocytes respond the most, with
hundreds of genes changing. Monocytes are the immune cells that eat and clear foreign
particles, so it makes sense they react the strongest. Second, the 200-nanometre
particles change more cell-type-specific genes than the 40-nanometre ones. Bigger
particles, bigger response.

The sixth stage is the main finding. I split the changed genes into four groups: genes
that respond only to 40 nanometres, only to 200 nanometres, to both, and genes that
respond only in the mixture. In monocytes, about 180 genes respond only when both sizes
are present, and to neither size on its own. In the lymphocytes the mixture is actually
weaker than either single size. So the combination of the two sizes does something the
individual sizes do not."

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

### 10. Additional analyses  ·  ~50s  ·  [Slide 12: list of all 5], [Slides 13-14: module scores, then mixture additivity]

"For the additional-insights part of the project I implemented five extra analyses,
which this slide lists. The first is module scoring, where I score stress and
inflammation gene programs per sample. The second is a mixture-additivity test, where I
check whether the mixture response equals the sum of the two single sizes; the plot on
the following slides shows it doesn't, the mixture sits above that expectation. The third
is a clustering robustness check, which confirms the clusters hold up across different
resolution settings. The fourth is a ligand-receptor analysis of how cell-to-cell
communication shifts on exposure. And the fifth is a dose-response analysis of disruption
against particle size. All five point the same way: a size-dependent, non-additive
response centred on monocytes."

### 11. Limitations and close  ·  ~30s  ·  [Slide 15: limitations]

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
