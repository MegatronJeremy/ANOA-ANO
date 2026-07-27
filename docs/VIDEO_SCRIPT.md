# Video presentation script (5-10 min)

Narration for the results video (deliverable). Slide numbers match the current
14-slide deck (`results/slides/GI_nanoplastic.pptx`). Target ~7-8 min total.
English (per the project's language preference).

---

### 1. Hook & question  ·  ~40s  ·  [Slide 1: title]

"Tiny plastic particles — nanoplastics — are now found in human blood. Once there, they
meet our immune cells directly. The open question: does *size* matter? Do small versus large
nanoplastics provoke different immune responses, and does a mixture of both do something
neither size does alone? We used single-cell RNA sequencing of human immune cells to find out."

### 2. Dataset & design  ·  ~45s  ·  [Slide 2: question & dataset]

"We have four samples from a single donor: peripheral blood immune cells exposed to 40-nanometre
particles, to 200-nanometre particles, to a 40-plus-200 mixture, and an unexposed control —
about 34,000 cells in total. One important caveat up front: one donor, one sample per condition,
so there are no biological replicates. That shapes every statistical choice we make, and we're
explicit about it throughout."

### 3. Pipeline overview  ·  ~50s  ·  [Slide 3: pipeline]

"Everything runs through one reproducible pipeline of six stages: quality control, batch
integration, cell-type annotation, composition analysis, differential expression, and
size-specific effects. One command reproduces all of it, and every stage has automated tests."

### 4. QC  ·  ~35s  ·  [Slides 4-5: QC violins, then counts-vs-mito scatter]

"First, quality control — filtering out broken cells before we trust anything. Rather than guess
thresholds, we justified them against the real distribution of all 34,000 cells. The upper gene
cutoff sits around the 99th percentile to drop likely doublets — two cells read as one — and the
mitochondrial cutoff between the 95th and 99th percentile to drop dying cells. The scatter on the
next slide shows those thresholds drawn over the actual data."

### 5. Integration  ·  ~40s  ·  [Slide 6: before/after UMAP]

"Next, integration. On the left, before batch correction, cells separate by sample — that's a
technical artifact, not biology. On the right, after Harmony correction, the samples mix and cells
group by biology instead. One honest caveat: because each condition is a single sample, batch and
treatment are confounded — so we run the actual differential expression on the *uncorrected* data,
where this doesn't bias the result."

### 6. Annotation  ·  ~45s  ·  [Slide 7: lineage UMAP + dotplot]

"We then labelled cell types with celltypist and canonical marker genes. Crucially, we
cross-checked our labels against two *independent* references already in the data — an Azimuth
annotation and the CoDi method. Agreement was about 93% with both, which gives real confidence
in the cell identities before we compare conditions."

### 7. Composition  ·  ~40s  ·  [Slide 8: composition bars]

"With cells labelled, we compared cell-type proportions across conditions. The clearest shift is
in monocytes under 200-nanometre exposure — roughly two-and-a-half times more monocytes than
control, while other cell types barely move. Because we have no replicates, we report these as
descriptive proportions and fold-changes; any statistical test here is exploratory only."

### 8. Differential expression & the main result  ·  ~70s  ·  [Slide 9: dose-response + size categories], [Slide 10: key finding]

"Now the core: within each cell type, which genes change on exposure? Two clear patterns emerge.
First, monocytes — the body's particle-eating cells — respond most strongly, with hundreds of
differentially expressed genes. Second, the 200-nanometre particles drive more cell-type-unique
genes than the 40-nanometre ones — bigger particles, bigger and more distinct response.

And the headline finding: the mixture is not simply the sum of its parts. In monocytes, the
mixture produces an *emergent* response — around 180 genes significant only when both sizes are
present, in neither size alone. Yet in lymphocytes the mixture is actually *weaker* than either
single size. So the combination genuinely does something new."

### 9. What the mixture does — biological interpretation  ·  ~55s  ·  [Slide 11: emergent inflammation]

"So what *is* that new thing? We ran pathway enrichment on those 180 mixture-only monocyte genes,
and the answer is strikingly clear: inflammation. The top pathways — at adjusted p-values around
ten-to-the-minus-twenty-one — are interleukin, cytokine and TNF signalling, and 'response to
lipopolysaccharide' — essentially the program a monocyte runs during a bacterial infection, here
triggered by plastic. The top mixture-only genes going up include RIPK2 and TRAF1 — core
innate-immune and TNF-signalling genes. Our interpretation: 40 and 200 nanometre particles each
cause modest, sub-threshold changes, but together they push the monocyte over an inflammatory
activation threshold — an emergent, non-additive effect. And since real exposure is always to
mixtures of sizes, studying single sizes alone may under-estimate the inflammatory risk."

### 10. Bonus analyses  ·  ~40s  ·  [Slides 12-13: module scores, then mixture additivity]

"We added five further analyses: scoring stress and inflammation gene programs, testing whether
the mixture is additive, checking that clustering is robust to parameter choices, a ligand-
receptor communication shift, and a dose-response of disruption versus particle size. They
reinforce the same story — a size-dependent, non-additive response centred on monocytes. The
additivity plot in particular confirms the mixture sits above the 40-plus-200 expectation."

### 11. Limitations & close  ·  ~30s  ·  [Slide 14: limitations]

"To be clear about limits: one donor and no replicates mean our differential expression is
cell-level and exploratory, with a pseudoreplication caveat, and there's no external gene-level
ground truth — so we validate internally and biologically. Everything is reproducible from the
repository with a single command and a full test suite. Thank you."

---

**YouTube title:** Do nanoplastic size and mixtures change the immune response? A single-cell study

**YouTube description (2-3 sentences):**
A single-cell RNA-seq analysis of human immune cells exposed to 40 nm, 200 nm, and mixed
polystyrene nanoplastics versus control. We find monocytes respond most strongly, larger
particles drive more cell-type-unique changes, and the 40+200 nm mixture produces an emergent,
inflammation-dominated monocyte response absent in either size alone. Code, figures, and a fully
reproducible pipeline are on GitHub.
