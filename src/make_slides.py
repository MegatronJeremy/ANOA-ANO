"""
Build the results slide deck (deliverable) from the figures the pipeline saves.

Reproducible: re-run after `--all` to rebuild results/slides/GI_nanoplastic.pptx
from whatever is currently in results/full/figures/. Missing figures are skipped with
a note rather than crashing, so a partial run still yields a deck.

    python -m src.make_slides          # or:  .\run.ps1 slides
"""
from pathlib import Path

from pptx import Presentation
from pptx.util import Inches, Pt

from . import config as cfg

SLIDES_DIR = cfg.PROJECT_ROOT / "results" / "slides"
SLIDES_DIR.mkdir(parents=True, exist_ok=True)

# 16:9 widescreen (the modern default). Layout math below derives from these
# so figure placement stays correct if the dimensions ever change.
SLIDE_W, SLIDE_H = 13.333, 7.5
MARGIN = 0.4   # left/right page margin (inches)

# The deck is a deliverable built from the FULL run's figures. Resolve to the
# full-run location explicitly (independent of any smoke run's cfg state).
FIG = cfg.RESULTS_DIR / "full" / "figures"


def _stretch(ph):
    """Widen a placeholder to the usable width WITHOUT disturbing its vertical
    position. python-pptx's built-in layouts size placeholders for a 4:3 slide
    (~9in wide); on our 16:9 slide that leaves text hugging the left with a dead
    band on the right. Capture the inherited top/height first: setting left/width
    materializes the shape's xfrm, and any position attribute left unset then
    defaults to 0 -- which previously yanked every text box to top=0, height=0."""
    top, height = ph.top, ph.height        # inherited from the layout; read before mutating
    ph.left = Inches(MARGIN)
    ph.width = Inches(SLIDE_W - 2 * MARGIN)
    ph.top = top
    ph.height = height


def _notes(slide, text):
    """Attach speaker notes (visible in PowerPoint's Presenter View, not on the
    slide itself). Used to carry the narration + 'what you're seeing' per slide."""
    if text:
        slide.notes_slide.notes_text_frame.text = text.strip()


def _title_slide(prs, title, subtitle, notes=""):
    s = prs.slides.add_slide(prs.slide_layouts[0])
    s.shapes.title.text = title
    s.placeholders[1].text = subtitle
    for ph in s.placeholders:
        _stretch(ph)
    _notes(s, notes)


def _bullets_slide(prs, title, bullets, notes=""):
    s = prs.slides.add_slide(prs.slide_layouts[1])
    s.shapes.title.text = title
    _stretch(s.shapes.title)
    body = s.placeholders[1]
    _stretch(body)
    tf = body.text_frame
    tf.clear()
    for i, b in enumerate(bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = b
        p.font.size = Pt(18)
    _notes(s, notes)


def _place_fit(slide, path, box_left, box_top, box_w, box_h):
    """Add a picture scaled to FIT inside (box_w x box_h) without distortion or
    overflow -- letterbox by whichever dimension binds, then centre it in the box.
    Sizing by height alone (the old behaviour) let wide figures (violin strips
    are 3.75:1, additivity 5:1) blow past the slide edge and overlap."""
    from PIL import Image
    px_w, px_h = Image.open(str(path)).size
    aspect = px_w / px_h            # width / height of the actual image
    box_aspect = box_w / box_h
    if aspect >= box_aspect:
        w, h = box_w, box_w / aspect   # width-bound (wide image)
    else:
        w, h = box_h * aspect, box_h   # height-bound (tall/square image)
    left = box_left + (box_w - w) / 2  # centre within the box
    top = box_top + (box_h - h) / 2
    slide.shapes.add_picture(str(path), Inches(left), Inches(top),
                             width=Inches(w), height=Inches(h))


def _figure_slide(prs, title, fig_names, caption="", notes=""):
    """One slide, title + up to two figures side by side (skips any missing).
    Positions derive from SLIDE_W / MARGIN so the layout tracks the slide size."""
    paths = [FIG / f for f in fig_names if (FIG / f).exists()]
    s = prs.slides.add_slide(prs.slide_layouts[5])
    s.shapes.title.text = title
    _stretch(s.shapes.title)
    usable_w = SLIDE_W - 2 * MARGIN
    if not paths:
        tb = s.shapes.add_textbox(Inches(MARGIN), Inches(2.5), Inches(usable_w), Inches(1)).text_frame
        tb.text = "(figure not generated yet -- run the pipeline)"
        _notes(s, notes)
        return
    top, box_h = 1.5, 5.0            # vertical band for figures (below title, above caption)
    if len(paths) == 1:
        _place_fit(s, paths[0], MARGIN, top, usable_w, box_h)   # one figure, full usable width
    else:
        gap = 0.3
        half = (usable_w - gap) / 2
        _place_fit(s, paths[0], MARGIN, top, half, box_h)              # left half
        _place_fit(s, paths[1], MARGIN + half + gap, top, half, box_h) # right half
    if caption:
        tb = s.shapes.add_textbox(Inches(MARGIN), Inches(6.9), Inches(usable_w), Inches(0.5)).text_frame
        tb.text = caption
        tb.paragraphs[0].font.size = Pt(12)
    _notes(s, notes)


def build():
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)

    _title_slide(prs, "Single-Cell Immune Response to Nanoplastic Particles",
                 "scRNA-seq of human PBMCs exposed to 40 nm / 200 nm / mixed polystyrene "
                 "nanoparticles vs control  -  Genomics Informatics, ETF",
                 notes=(
                     "Hi, my name is Vuk Djordjevic. I'm a master's student at the School of "
                     "Electrical Engineering in Belgrade, and this is my project for Genome "
                     "Informatics, 2026: a single-cell analysis of how human immune cells respond "
                     "to nanoplastic particles of different sizes. Nanoplastics are tiny plastic "
                     "particles, now found in human blood in direct contact with our immune cells. "
                     "My question is whether size matters: do small and large nanoplastics cause "
                     "different responses, and does a mixture of both do something neither size does "
                     "alone? To answer it I used single-cell RNA sequencing, which reads which genes "
                     "are active in each individual cell."))
    _bullets_slide(prs, "Question & dataset", [
        "Do small vs large nanoplastics provoke different immune responses?",
        "Does the 40+200 nm mixture do something neither size does alone?",
        "4 samples, one donor: PSNP_40nm, PSNP_200nm, PSNP_mixture, control",
        "~34,000 PBMCs, AnnData/.h5ad (Zenodo 10.5281/zenodo.15866724)",
    ], notes=(
        "The data is four samples from a single donor: peripheral blood immune cells exposed to "
        "40-nanometre particles, to 200-nanometre particles, to a mixture of both, and one "
        "unexposed control. About 34,000 cells in total. One thing to be upfront about: one donor, "
        "one sample per condition, so no biological replicates. That shaped how I did the "
        "statistics, and I flag it throughout."))
    _bullets_slide(prs, "Pipeline (6 stages, one reproducible driver)", [
        "1. QC & preprocessing  -  thresholds justified on real percentiles",
        "2. Integration  -  Harmony batch correction, UMAP, Leiden",
        "3. Annotation  -  celltypist, cross-checked vs Azimuth & CoDi (~93%)",
        "4. Composition  -  proportions + log2 fold-change vs control",
        "5. Differential expression  -  per-lineage Wilcoxon + pathway enrichment",
        "6. Size-specific effects  -  unique / shared / mixture-emergent genes",
    ], notes=(
        "The whole analysis runs through one pipeline of six stages: quality control, integration, "
        "cell-type annotation, composition analysis, differential expression, and size-specific "
        "effects. A single command reproduces all of it from the raw data, and each stage has "
        "automated tests."))
    # QC split across two slides -- the violin strip (3.75:1) and the scatter
    # are both too wide to share one slide without shrinking to unreadable.
    _figure_slide(prs, "QC & preprocessing (1/2): per-sample QC metrics", ["01_qc_violin_after.png"],
                  "Violins of the 3 QC metrics after filtering. Thresholds justified on the real distribution "
                  "(max_genes ~p99, max_mito between p95-p99).",
                  notes=(
                      "The first stage is quality control: removing low-quality cells before trusting "
                      "anything. These three plots are violins, basically histograms on their side, fatter "
                      "where more cells sit; along the bottom are my four samples. Left is how many genes each "
                      "cell expressed, middle is total reads per cell, right is the percent of mitochondrial "
                      "reads, a sign of a dying cell. This is after cleaning, and the shapes look healthy and "
                      "similar across all samples. I set the cutoffs from the real distribution, not by "
                      "guessing: drop cells with too many genes, usually two stuck together, and cells with too "
                      "much mitochondria, which are dying."))
    _figure_slide(prs, "QC & preprocessing (2/2): counts vs %mito", ["01_qc_scatter_counts_mito.png"],
                  "total_counts vs %mito (pre-filter); threshold lines drawn. Drops likely doublets and dying cells.",
                  notes=(
                      "This is the same cleaning step, drawn a different way. Every dot here is one cell. "
                      "Left-to-right is how many reads the cell has, and bottom-to-top is how mitochondrial it "
                      "is. The red dashed line is my cutoff at 15 percent. Almost all the cells form a dense "
                      "cloud hugging the bottom, with low mitochondria, which is exactly what healthy cells "
                      "look like. The handful of dots floating up above the red line are dying or broken cells, "
                      "and those get thrown out."))
    _figure_slide(prs, "Integration removes the batch effect", ["02_umap_pre_harmony_by_sample.png", "02_umap_by_sample.png"],
                  "Before (left) vs after Harmony (right), coloured by sample: samples mix after correction.",
                  notes=(
                      "The second stage is integration. Both pictures are a UMAP: each cell has thousands of "
                      "genes, and UMAP squashes that down to one dot on a map, so cells with similar gene "
                      "activity land close together. The axes are just map coordinates; only closeness matters. "
                      "Each dot is a cell, coloured by sample. I want the colours well mixed: if cells clumped "
                      "by colour, the sample would matter more than the biology, a technical artefact. Here "
                      "they're nicely blended, so cells group by cell type, not batch. One honest caveat: with "
                      "a single sample per condition, batch and treatment are tied together, so for the "
                      "gene-level results later I use the uncorrected data, where this can't bias anything."))
    _figure_slide(prs, "Cell-type annotation", ["03_umap_lineage.png", "03_marker_dotplot.png"],
                  "celltypist lineages; agreement with Azimuth 92.7% and CoDi 93.1% (independent references).",
                  notes=(
                      "The third stage is labelling the cell types. On the left, the same cell-map coloured by "
                      "cell type: T cells the big purple mass on the right, monocytes the orange island "
                      "top-left, B cells blue at the bottom, NK cells green. Clean separate islands mean the "
                      "types are real. On the right, a dot-plot that proves the labels: rows are clusters, "
                      "columns are marker genes like CD3D for T cells or CD14 for monocytes. A big dark dot "
                      "means that gene is strongly on in that cluster, and the dark dots line up on the "
                      "diagonal, so each cluster lights up its own marker. I labelled with a tool called "
                      "celltypist and cross-checked against two independent methods; they agreed at about 93 "
                      "percent."))
    _figure_slide(prs, "Composition shifts vs control", ["04_composition_stacked.png", "04_composition_grouped.png"],
                  "Cell-type proportions per sample; PSNP_200nm shows the largest compositional shift.",
                  notes=(
                      "The fourth stage asks whether exposure changed how many of each cell type there are. "
                      "The bars show what fraction of a sample is each cell type, and the four bars per group "
                      "are the four samples. T cells dominate every sample, the tall bars on the right at "
                      "around 80 percent, which is normal for blood. The thing to notice is the Monocyte "
                      "group: the green 200-nanometre bar is noticeably taller, roughly two and a half times "
                      "the control's monocyte fraction, while everything else barely moves. So 200 nanometres "
                      "specifically bumps up the monocytes. With no replicates, I report this as a descriptive "
                      "shift and treat any test as exploratory."))
    _figure_slide(prs, "Differential expression: dose-response", ["07_dose_response.png", "06_size_categories.png"],
                  "Monocytes respond most; 200 nm drives more unique genes than 40 nm; mixture-emergent in monocytes.",
                  notes=(
                      "The fifth stage asks which genes actually change on exposure, within each cell type. In "
                      "both bar charts, bar height is the number of genes that significantly changed. The left "
                      "chart groups by condition, coloured by cell type: the orange Monocyte bars tower over "
                      "everything, because monocytes eat and clear foreign particles, and the 200-nanometre bar "
                      "beats the 40-nanometre one, so bigger particles, bigger response. The right chart is the "
                      "sixth stage, showing where those genes overlap. Four bars per cell type: genes reacting "
                      "only to 40, only to 200, to both, and the key one in red, genes that fire only in the "
                      "mixture and in neither size alone. In monocytes that red bar is about 180 genes, a whole "
                      "response that only appears when both sizes are combined."))
    _bullets_slide(prs, "Key biological finding", [
        "Monocytes show the strongest transcriptional response to nanoplastic.",
        "200 nm particles drive MORE lineage-unique genes than 40 nm.",
        "The 40+200 nm mixture produces an EMERGENT monocyte response -",
        "  genes significant only in the mixture, in neither single size.",
        "In lymphocytes the mixture is weaker than either size (sub-additive).",
    ], notes=(
        "So let me pull the main findings together. First, monocytes show the strongest response "
        "to nanoplastic of any cell type. Second, the 200-nanometre particles drive more "
        "cell-type-unique genes than the 40-nanometre ones. And third, the biggest one: the 40 "
        "plus 200-nanometre mixture produces an emergent monocyte response, about 180 genes that "
        "are significant only in the mixture and in neither single size. Interestingly, in the "
        "lymphocytes the opposite happens, the mixture is actually weaker than either size alone. "
        "Either way, the combination of the two sizes does something the individual sizes do "
        "not."))
    _bullets_slide(prs, "What the mixture does: emergent inflammation", [
        "The 180 mixture-only monocyte genes are dominated by INFLAMMATION.",
        "Pathway enrichment (p ~ 1e-21): interleukin, cytokine & TNF signalling,",
        "  and 'response to lipopolysaccharide' - i.e. a bacterial-attack-like program.",
        "Top mixture-only genes up: RIPK2, TRAF1, PIK3CB (innate/TNF signalling).",
        "Interpretation: two sub-threshold exposures together push monocytes over",
        "  an inflammatory activation threshold - an emergent, non-additive effect.",
        "Real exposure is always to mixtures -> single-size studies may under-estimate risk.",
    ], notes=(
        "So what is that new thing? I ran pathway enrichment on those 180 mixture-only monocyte "
        "genes, and it points clearly at inflammation. The top pathways, with strong statistical "
        "support, are interleukin, cytokine, and TNF signalling, plus the response to "
        "lipopolysaccharide. That last one is the program a monocyte runs during a bacterial "
        "infection, except here it's set off by plastic. Two of the strongest genes going up are "
        "RIPK2 and TRAF1, core innate-immune and TNF genes. My reading: 40 and 200 nanometres each "
        "cause small, sub-threshold changes, but together they push the monocyte over the edge into "
        "an inflammatory state. And since real exposure is always to a mix of sizes, testing single "
        "sizes alone could underestimate the risk."))
    # Bonus split across two slides -- the additivity plot (5:1) is far too wide
    # to sit beside the module-scores heatmap on one slide.
    _bullets_slide(prs, "Additional analyses (5 implemented)", [
        "1. Module scoring - stress & inflammation gene programs per sample (shown next).",
        "2. Mixture additivity - is the mixture = 40nm + 200nm? (shown next).",
        "3. Clustering robustness - clusters stable across Leiden resolutions (ARI).",
        "4. Ligand-receptor - shift in cell-to-cell communication on exposure.",
        "5. Dose-response - disruption magnitude vs particle size.",
        "All five reinforce the same story: size-dependent, non-additive, monocyte-centred.",
    ], notes=(
        "For the additional-insights part of the project I implemented five extra analyses, which "
        "this slide lists, and the next slides show each one. They are module scoring, a mixture "
        "additivity test, a clustering robustness check, a ligand-receptor communication analysis, "
        "and dose-response. I'll walk through them in the next few slides. All five point the same "
        "way: a size-dependent, non-additive response centred on the monocytes."))
    _figure_slide(prs, "Additional analyses: stress/inflammation module scores", ["07_module_scores.png"],
                  "Analysis 1 of 5. Per-sample mean module scores for stress and inflammation gene programs.",
                  notes=(
                      "The first additional analysis is module scoring: a colour grid where each row is a "
                      "sample and each column is a biological program, a named group of genes: oxidative "
                      "stress, NF-kB inflammation, interferon, heat shock. Red means the program is turned up. "
                      "The whole oxidative-stress column goes deep red in all three exposed samples but not the "
                      "control, so plastic reliably stresses the cells, and the inflammation column turns pink "
                      "specifically under 200 nanometres. So the size effect shows up here too, in independent "
                      "gene programs."))
    _figure_slide(prs, "Additional analyses: mixture additivity", ["07_mixture_additivity.png"],
                  "Analysis 2 of 5. Observed mixture response vs the 40+200 nm additive expectation, per lineage.",
                  notes=(
                      "The second analysis is a direct test of whether the mixture is just 40 plus 200 added "
                      "together. One scatter per cell type, each dot a gene: left-to-right is what I'd predict "
                      "if the mixture were the two sizes summed, bottom-to-top is what it actually did, and the "
                      "red diagonal is where the prediction would be perfect. If every gene sat on that line, "
                      "the mixture would be boringly additive. But in the Monocyte panel there are clear blobs "
                      "of genes off the line, genes the mixture switched on or off that you'd never predict "
                      "from the single sizes. That off-the-line behaviour is exactly the emergent, non-additive "
                      "effect."))
    _figure_slide(prs, "Additional analyses: robustness & communication",
                  ["07_clustering_robustness.png", "07_ligand_receptor.png"],
                  "Analyses 3 & 4 of 5. Left: clustering robustness (ARI across Leiden resolutions). "
                  "Right: ligand-receptor communication shift on exposure. (Analysis 5, dose-response, is on the DE slide.)",
                  notes=(
                      "The third and fourth analyses. The left plot is a sanity check on my clustering: the "
                      "green dashed line is how much the grouping still agrees with my default setting, where "
                      "1.0 is identical. It stays high, so my cell groups aren't a fluke of one setting. The "
                      "right plot is cell-to-cell messaging: each row is a signalling pair, and a bar goes "
                      "right if that message got louder after exposure. The standout is IL1B, a classic "
                      "inflammation alarm, shooting far right in all exposed samples and furthest under 200 "
                      "nanometres, so the cells are shouting inflammation at each other. The fifth analysis, "
                      "dose-response, I already showed on the earlier DE slide."))
    _bullets_slide(prs, "Limitations & reproducibility", [
        "One donor, one sample per condition -> NO biological replicates.",
        "DE is cell-level Wilcoxon (not pseudobulk); p-values exploratory, pseudoreplication caveat.",
        "No external gene-level ground truth; DE validated internally + biologically.",
        "Everything reproduces: `python run_pipeline.py --all`; 37 offline intent tests.",
    ], notes=(
        "Finally, to be clear about the limits. One donor and no replicates mean the differential "
        "expression is done per cell and is exploratory, with a pseudoreplication caveat. There's "
        "no external gene-level reference, so I validated the results internally and biologically. "
        "But everything reproduces from the repository with a single command and a full test "
        "suite. Thanks for watching."))

    out = SLIDES_DIR / "GI_nanoplastic.pptx"
    prs.save(out)
    print(f"saved slide deck -> {out}  ({len(list(prs.slides))} slides)")
    return out


if __name__ == "__main__":
    build()
