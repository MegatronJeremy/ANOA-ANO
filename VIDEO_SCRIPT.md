# Video presentation script (5-10 min)

Narration for the results video (deliverable). Each section notes the slide on screen
and a target duration; total ~7 min. English (per the project's language preference).

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

### 4. QC  ·  ~30s  ·  [Slide 4: QC figures]

"First, quality control. Rather than guess thresholds, we justified them against the real
distribution of all 34,000 cells — for example, the upper gene cutoff sits around the 99th
percentile to drop likely doublets, and the mitochondrial cutoff between the 95th and 99th
percentile to drop dying cells."

### 5. Integration  ·  ~40s  ·  [Slide 5: before/after UMAP]

"Next, integration. On the left, before batch correction, cells separate by sample — that's a
technical artifact. On the right, after Harmony correction, the samples mix and cells group by
biology instead. We confirmed this quantitatively: clusters are blends of all four samples."

### 6. Annotation  ·  ~45s  ·  [Slide 6: lineage UMAP + dotplot]

"We then labelled cell types with celltypist and canonical marker genes. Crucially, we
cross-checked our labels against two *independent* references already in the data — an Azimuth
annotation and the CoDi method. Agreement was about 93% with both, which gives real confidence
in the cell identities before we compare conditions."

### 7. Composition  ·  ~40s  ·  [Slide 7: composition bars]

"With cells labelled, we compared cell-type proportions across conditions. The 200-nanometre
exposure shows the largest compositional shift relative to control. Because we have no
replicates, we report these as descriptive proportions and fold-changes — any statistical test
here is exploratory only."

### 8. Differential expression & the main result  ·  ~70s  ·  [Slide 8: dose-response + size categories], [Slide 9: key finding]

"Now the core: within each cell type, which genes change on exposure? Two clear patterns emerge.
First, monocytes — the body's particle-eating cells — respond most strongly, with hundreds of
differentially expressed genes. Second, the 200-nanometre particles drive more cell-type-unique
genes than the 40-nanometre ones.

And the headline finding: the mixture is not simply the sum of its parts. In monocytes, the
mixture produces an *emergent* response — genes that are significant only when both sizes are
present, in neither size alone. Yet in lymphocytes the mixture is actually *weaker* than either
single size. So the combination genuinely does something new."

### 9. Bonus analyses  ·  ~40s  ·  [Slide 10: module scores + additivity]

"We added five further analyses: scoring stress and inflammation gene programs, testing whether
the mixture is additive, checking that clustering is robust to parameter choices, a ligand-
receptor communication shift, and a dose-response of disruption versus particle size. These
reinforce the same story — a size-dependent, non-additive response centred on monocytes."

### 10. Limitations & close  ·  ~30s  ·  [Slide 11: limitations]

"To be clear about limits: one donor and no replicates mean our differential expression is
cell-level and exploratory, with a pseudoreplication caveat, and there's no external gene-level
ground truth — so we validate internally and biologically. Everything is reproducible from the
repository with a single command and a full test suite. Thank you."

---

**YouTube title:** Do nanoplastic size and mixtures change the immune response? A single-cell study

**YouTube description (2-3 sentences):**
A single-cell RNA-seq analysis of human immune cells exposed to 40 nm, 200 nm, and mixed
polystyrene nanoplastics versus control. We find monocytes respond most strongly, larger
particles drive more cell-type-unique changes, and the mixture produces an emergent monocyte
response absent in either size alone. Code, figures, and a fully reproducible pipeline are on
GitHub.
