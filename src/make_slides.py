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
                     "Intro: I'm Vuk Djordjevic, master's student at the School of Electrical "
                     "Engineering, Genome Informatics 2026. This is a single-cell RNA-seq study of "
                     "how human immune cells respond to nanoplastic particles of different sizes. "
                     "The core question: does particle size matter, and does a mixture of sizes do "
                     "something neither size does alone?"))
    _bullets_slide(prs, "Question & dataset", [
        "Do small vs large nanoplastics provoke different immune responses?",
        "Does the 40+200 nm mixture do something neither size does alone?",
        "4 samples, one donor: PSNP_40nm, PSNP_200nm, PSNP_mixture, control",
        "~34,000 PBMCs, AnnData/.h5ad (Zenodo 10.5281/zenodo.15866724)",
    ], notes=(
        "The data is 4 samples from one donor: peripheral blood immune cells exposed to 40 nm "
        "particles, 200 nm particles, a 40+200 nm mixture, and an unexposed control. About "
        "34,000 cells total. Be upfront: one donor, one sample per condition, so no biological "
        "replicates -- that shaped the statistics and I flag it throughout."))
    _bullets_slide(prs, "Pipeline (6 stages, one reproducible driver)", [
        "1. QC & preprocessing  -  thresholds justified on real percentiles",
        "2. Integration  -  Harmony batch correction, UMAP, Leiden",
        "3. Annotation  -  celltypist, cross-checked vs Azimuth & CoDi (~93%)",
        "4. Composition  -  proportions + log2 fold-change vs control",
        "5. Differential expression  -  per-lineage Wilcoxon + pathway enrichment",
        "6. Size-specific effects  -  unique / shared / mixture-emergent genes",
    ], notes=(
        "The analysis is one pipeline of six stages: QC, integration, annotation, composition, "
        "differential expression, and size-specific effects. A single command reproduces all of "
        "it from raw data, and each stage has automated tests."))
    # QC split across two slides -- the violin strip (3.75:1) and the scatter
    # are both too wide to share one slide without shrinking to unreadable.
    _figure_slide(prs, "QC & preprocessing (1/2): per-sample QC metrics", ["01_qc_violin_after.png"],
                  "Violins of the 3 QC metrics after filtering. Thresholds justified on the real distribution "
                  "(max_genes ~p99, max_mito between p95-p99).",
                  notes=(
                      "What you're looking at: three 'violin' plots side by side, one per quality check. A "
                      "violin is just a histogram turned on its side and mirrored -- fatter where more cells "
                      "sit. Along the bottom of each are my four samples (40 nm, 200 nm, mixture, control). "
                      "Left panel = how many genes each cell expressed, middle = total reads per cell, right = "
                      "percent of reads coming from mitochondria (a sign of a dying cell). What to say: this is "
                      "AFTER cleaning. The shapes look healthy and similar across all four samples -- no sample "
                      "is junk. I set the cutoffs from the real data, not by guessing: drop cells with too many "
                      "genes (likely two cells stuck together) or too much mitochondria (dying)."))
    _figure_slide(prs, "QC & preprocessing (2/2): counts vs %mito", ["01_qc_scatter_counts_mito.png"],
                  "total_counts vs %mito (pre-filter); threshold lines drawn. Drops likely doublets and dying cells.",
                  notes=(
                      "What you're looking at: every dot is one cell. Left-to-right is how many reads the cell "
                      "has; bottom-to-top is how mitochondrial it is. The red dashed line is my cutoff at 15 "
                      "percent mito. What to say: almost all cells form a dense cloud hugging the bottom -- "
                      "low mito, which is what healthy cells look like. The handful of dots floating up above "
                      "the red line are dying or broken cells, and those get thrown out. This is the same "
                      "cleaning as the previous slide, just drawn as a cloud so you can see the cutoff."))
    _figure_slide(prs, "Integration removes the batch effect", ["02_umap_pre_harmony_by_sample.png", "02_umap_by_sample.png"],
                  "Before (left) vs after Harmony (right), coloured by sample: samples mix after correction.",
                  notes=(
                      "What you're looking at: both pictures are a UMAP -- think of it as squashing each cell's "
                      "thousands of genes down to a single dot on a 2D map, so that cells with similar genes "
                      "land near each other. The axes are just map coordinates, no units; only closeness "
                      "matters. Each dot is a cell, coloured by which sample it came from. What to say: I want "
                      "the colours to be MIXED -- if cells clumped by colour it would mean the sample it came "
                      "from matters more than its biology, which is a technical artefact. On both sides here "
                      "the colours are already well blended, so batches line up and cells group by TYPE, not "
                      "by sample. (Honest caveat: one sample per condition, so for the gene results later I use "
                      "the uncorrected data to stay safe.)"))
    _figure_slide(prs, "Cell-type annotation", ["03_umap_lineage.png", "03_marker_dotplot.png"],
                  "celltypist lineages; agreement with Azimuth 92.7% and CoDi 93.1% (independent references).",
                  notes=(
                      "What you're looking at: LEFT is the same cell-map, now coloured by what KIND of immune "
                      "cell each dot is. You can literally see the neighbourhoods: T cells are the big purple "
                      "mass on the right, monocytes the orange island top-left, B cells blue at the bottom, NK "
                      "cells green. Clean separate islands = the cell types are real and distinct. RIGHT is a "
                      "dot-plot that PROVES the labels: rows are cell clusters, columns are famous marker genes "
                      "(e.g. CD3D for T cells, CD14 for monocytes). A big dark dot means 'this gene is strongly "
                      "on in this cluster'. The dark dots line up on the diagonal, i.e. each cluster lights up "
                      "its own known marker. What to say: I labelled with a tool called celltypist and it "
                      "agreed with two other independent methods at about 93 percent."))
    _figure_slide(prs, "Composition shifts vs control", ["04_composition_stacked.png", "04_composition_grouped.png"],
                  "Cell-type proportions per sample; PSNP_200nm shows the largest compositional shift.",
                  notes=(
                      "What you're looking at: this asks 'did exposure change HOW MANY of each cell type there "
                      "are?' The bars show what fraction of a sample is each cell type; the four coloured bars "
                      "in each group are the four samples. What to say: T cells dominate every sample (the tall "
                      "bars on the right, ~80 percent) -- that's normal for blood. The eye-catch is the "
                      "Monocyte group: the 200 nm bar (green) is noticeably taller than the others, roughly "
                      "two-and-a-half times the control's monocyte fraction, while everything else barely "
                      "moves. So 200 nm exposure specifically bumps up monocytes. No replicates, so I call this "
                      "a descriptive shift, not a p-value."))
    _figure_slide(prs, "Differential expression: dose-response", ["07_dose_response.png", "06_size_categories.png"],
                  "Monocytes respond most; 200 nm drives more unique genes than 40 nm; mixture-emergent in monocytes.",
                  notes=(
                      "What you're looking at: two bar charts. Both use bar HEIGHT = number of genes that "
                      "significantly changed (taller = bigger reaction). LEFT groups by condition (40 nm, 200 "
                      "nm, mixture), coloured by cell type. The story jumps out: the orange Monocyte bars "
                      "tower over everything -- monocytes are the cells that eat foreign junk, so they react "
                      "hardest -- and 200 nm is taller than 40 nm, meaning bigger particles, bigger response. "
                      "RIGHT zooms into WHERE those genes overlap. For each cell type, four bars: genes that "
                      "react only to 40 nm (blue), only to 200 nm (orange), to both (green), and the key one, "
                      "RED = genes that fire ONLY in the mixture and in neither size alone. What to say: look "
                      "at Monocyte -- that red bar is about 180 genes, a whole response that only appears when "
                      "both sizes are combined."))
    _bullets_slide(prs, "Key biological finding", [
        "Monocytes show the strongest transcriptional response to nanoplastic.",
        "200 nm particles drive MORE lineage-unique genes than 40 nm.",
        "The 40+200 nm mixture produces an EMERGENT monocyte response -",
        "  genes significant only in the mixture, in neither single size.",
        "In lymphocytes the mixture is weaker than either size (sub-additive).",
    ], notes=(
        "This is the headline. Three points: monocytes respond strongest; 200 nm drives more "
        "lineage-unique genes than 40 nm; and the 40+200 nm mixture produces an emergent "
        "monocyte response -- about 180 genes significant only in the mixture, in neither single "
        "size. In lymphocytes the opposite happens: the mixture is weaker than either size alone. "
        "So the combination of sizes does something the individual sizes do not."))
    _bullets_slide(prs, "What the mixture does: emergent inflammation", [
        "The 180 mixture-only monocyte genes are dominated by INFLAMMATION.",
        "Pathway enrichment (p ~ 1e-21): interleukin, cytokine & TNF signalling,",
        "  and 'response to lipopolysaccharide' - i.e. a bacterial-attack-like program.",
        "Top mixture-only genes up: RIPK2, TRAF1, PIK3CB (innate/TNF signalling).",
        "Interpretation: two sub-threshold exposures together push monocytes over",
        "  an inflammatory activation threshold - an emergent, non-additive effect.",
        "Real exposure is always to mixtures -> single-size studies may under-estimate risk.",
    ], notes=(
        "What is that new thing? I ran pathway enrichment on the 180 mixture-only monocyte genes "
        "and it points clearly at inflammation. Top pathways, adjusted p about 1e-21: interleukin, "
        "cytokine and TNF signalling, plus 'response to lipopolysaccharide' -- the program a "
        "monocyte runs during a bacterial infection, here set off by plastic. Strongest genes up "
        "are RIPK2 and TRAF1, core innate-immune and TNF genes. Reading: 40 and 200 nm each cause "
        "sub-threshold changes, but together they push monocytes over an inflammatory threshold. "
        "Real exposure is always to mixtures, so single-size testing may under-estimate the risk."))
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
        "For the additional-insights part I implemented five extra analyses; this slide lists them "
        "and the next slides show each. Module scoring, mixture additivity, clustering robustness, "
        "ligand-receptor communication, and dose-response. All five point the same way: a "
        "size-dependent, non-additive response centred on monocytes."))
    _figure_slide(prs, "Additional analyses: stress/inflammation module scores", ["07_module_scores.png"],
                  "Analysis 1 of 5. Per-sample mean module scores for stress and inflammation gene programs.",
                  notes=(
                      "Analysis 1 of 5. What you're looking at: a colour grid (heatmap). Each ROW is a sample, "
                      "each COLUMN is a biological 'program' -- a named group of genes (oxidative stress, NF-kB "
                      "inflammation, interferon, heat shock). Red = that program is turned UP in that sample, "
                      "blue = turned down, white = neutral. What to say: the whole left column, oxidative "
                      "stress, goes deep red in all three exposed samples but not the control -- so plastic "
                      "reliably stresses the cells. And the NF-kB inflammation column turns pink specifically "
                      "under 200 nm. So the size effect I keep pointing at shows up here too, in independent "
                      "gene programs."))
    _figure_slide(prs, "Additional analyses: mixture additivity", ["07_mixture_additivity.png"],
                  "Analysis 2 of 5. Observed mixture response vs the 40+200 nm additive expectation, per lineage.",
                  notes=(
                      "Analysis 2 of 5 -- the direct test of 'is the mixture just 40 plus 200?'. What you're "
                      "looking at: one scatter per cell type. Each dot is a gene. Left-to-right is what I'd "
                      "PREDICT if the mixture were simply the two sizes added together; bottom-to-top is what "
                      "the mixture ACTUALLY did. The red diagonal line is 'prediction was perfect'. What to "
                      "say: if everything sat on the line, the mixture would be boringly additive. But look at "
                      "the Monocyte panel -- there are clear blobs of genes sitting OFF the line, especially "
                      "the arms pointing up-and-left and the flat streaks near zero. Those are genes the "
                      "mixture switched on (or off) that you'd never predict from the single sizes. That "
                      "off-the-line behaviour IS the emergent, non-additive effect."))
    _figure_slide(prs, "Additional analyses: robustness & communication",
                  ["07_clustering_robustness.png", "07_ligand_receptor.png"],
                  "Analyses 3 & 4 of 5. Left: clustering robustness (ARI across Leiden resolutions). "
                  "Right: ligand-receptor communication shift on exposure. (Analysis 5, dose-response, is on the DE slide.)",
                  notes=(
                      "Analyses 3 and 4, two plots. LEFT is a sanity check on my clustering. Bottom axis = how "
                      "finely I chose to cut the cells into groups ('resolution'). Blue line (rising) = you get "
                      "more clusters if you cut finer -- obvious. The one that matters is the GREEN dashed "
                      "line: how much the grouping still agrees with my default setting, where 1.0 = identical. "
                      "It stays high near my chosen resolution of 1.0, so my cell groups aren't a fluke of one "
                      "knob setting. RIGHT is cell-to-cell messaging. Each row is a signalling pair (a "
                      "messenger and its receptor); bars go RIGHT if that message got louder after exposure, "
                      "LEFT if quieter. What to say: the standout is IL1B-IL1R1, a classic inflammation alarm "
                      "signal, shooting far right in all exposed samples and furthest under 200 nm -- so the "
                      "cells aren't just changing internally, they're shouting inflammation at each other. "
                      "(Analysis 5, dose-response, was the left chart on the earlier DE slide.)"))
    _bullets_slide(prs, "Limitations & reproducibility", [
        "One donor, one sample per condition -> NO biological replicates.",
        "DE is cell-level Wilcoxon (not pseudobulk); p-values exploratory, pseudoreplication caveat.",
        "No external gene-level ground truth; DE validated internally + biologically.",
        "Everything reproduces: `python run_pipeline.py --all`; 37 offline intent tests.",
    ], notes=(
        "The limits, stated plainly: one donor and no replicates mean the differential expression "
        "is per-cell and exploratory, with a pseudoreplication caveat. No external gene-level "
        "reference, so I validated internally and biologically. Everything reproduces from the "
        "repository with a single command and a 37-test suite. Thanks for watching."))

    out = SLIDES_DIR / "GI_nanoplastic.pptx"
    prs.save(out)
    print(f"saved slide deck -> {out}  ({len(list(prs.slides))} slides)")
    return out


if __name__ == "__main__":
    build()
