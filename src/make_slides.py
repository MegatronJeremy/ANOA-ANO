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


def _title_slide(prs, title, subtitle):
    s = prs.slides.add_slide(prs.slide_layouts[0])
    s.shapes.title.text = title
    s.placeholders[1].text = subtitle
    for ph in s.placeholders:
        _stretch(ph)


def _bullets_slide(prs, title, bullets):
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


def _figure_slide(prs, title, fig_names, caption=""):
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


def build():
    prs = Presentation()
    prs.slide_width = Inches(SLIDE_W)
    prs.slide_height = Inches(SLIDE_H)

    _title_slide(prs, "Single-Cell Immune Response to Nanoplastic Particles",
                 "scRNA-seq of human PBMCs exposed to 40 nm / 200 nm / mixed polystyrene "
                 "nanoparticles vs control  -  Genomics Informatics, ETF")
    _bullets_slide(prs, "Question & dataset", [
        "Do small vs large nanoplastics provoke different immune responses?",
        "Does the 40+200 nm mixture do something neither size does alone?",
        "4 samples, one donor: PSNP_40nm, PSNP_200nm, PSNP_mixture, control",
        "~34,000 PBMCs, AnnData/.h5ad (Zenodo 10.5281/zenodo.15866724)",
    ])
    _bullets_slide(prs, "Pipeline (6 stages, one reproducible driver)", [
        "1. QC & preprocessing  -  thresholds justified on real percentiles",
        "2. Integration  -  Harmony batch correction, UMAP, Leiden",
        "3. Annotation  -  celltypist, cross-checked vs Azimuth & CoDi (~93%)",
        "4. Composition  -  proportions + log2 fold-change vs control",
        "5. Differential expression  -  per-lineage Wilcoxon + pathway enrichment",
        "6. Size-specific effects  -  unique / shared / mixture-emergent genes",
    ])
    # QC split across two slides -- the violin strip (3.75:1) and the scatter
    # are both too wide to share one slide without shrinking to unreadable.
    _figure_slide(prs, "QC & preprocessing (1/2): per-sample QC metrics", ["01_qc_violin_after.png"],
                  "Violins of the 3 QC metrics after filtering. Thresholds justified on the real distribution "
                  "(max_genes ~p99, max_mito between p95-p99).")
    _figure_slide(prs, "QC & preprocessing (2/2): counts vs %mito", ["01_qc_scatter_counts_mito.png"],
                  "total_counts vs %mito (pre-filter); threshold lines drawn. Drops likely doublets and dying cells.")
    _figure_slide(prs, "Integration removes the batch effect", ["02_umap_pre_harmony_by_sample.png", "02_umap_by_sample.png"],
                  "Before (left) vs after Harmony (right), coloured by sample: samples mix after correction.")
    _figure_slide(prs, "Cell-type annotation", ["03_umap_lineage.png", "03_marker_dotplot.png"],
                  "celltypist lineages; agreement with Azimuth 92.7% and CoDi 93.1% (independent references).")
    _figure_slide(prs, "Composition shifts vs control", ["04_composition_stacked.png", "04_composition_grouped.png"],
                  "Cell-type proportions per sample; PSNP_200nm shows the largest compositional shift.")
    _figure_slide(prs, "Differential expression: dose-response", ["07_dose_response.png", "06_size_categories.png"],
                  "Monocytes respond most; 200 nm drives more unique genes than 40 nm; mixture-emergent in monocytes.")
    _bullets_slide(prs, "Key biological finding", [
        "Monocytes show the strongest transcriptional response to nanoplastic.",
        "200 nm particles drive MORE lineage-unique genes than 40 nm.",
        "The 40+200 nm mixture produces an EMERGENT monocyte response -",
        "  genes significant only in the mixture, in neither single size.",
        "In lymphocytes the mixture is weaker than either size (sub-additive).",
    ])
    _bullets_slide(prs, "What the mixture does: emergent inflammation", [
        "The 180 mixture-only monocyte genes are dominated by INFLAMMATION.",
        "Pathway enrichment (p ~ 1e-21): interleukin, cytokine & TNF signalling,",
        "  and 'response to lipopolysaccharide' - i.e. a bacterial-attack-like program.",
        "Top mixture-only genes up: RIPK2, TRAF1, PIK3CB (innate/TNF signalling).",
        "Interpretation: two sub-threshold exposures together push monocytes over",
        "  an inflammatory activation threshold - an emergent, non-additive effect.",
        "Real exposure is always to mixtures -> single-size studies may under-estimate risk.",
    ])
    # Bonus split across two slides -- the additivity plot (5:1) is far too wide
    # to sit beside the module-scores heatmap on one slide.
    _figure_slide(prs, "Bonus analyses (1/2): stress/inflammation module scores", ["07_module_scores.png"],
                  "Per-sample mean module scores for stress and inflammation gene programs.")
    _figure_slide(prs, "Bonus analyses (2/2): mixture additivity", ["07_mixture_additivity.png"],
                  "Observed mixture response vs the 40+200 nm additive expectation, per lineage.")
    _bullets_slide(prs, "Limitations & reproducibility", [
        "One donor, one sample per condition -> NO biological replicates.",
        "DE is cell-level Wilcoxon (not pseudobulk); p-values exploratory, pseudoreplication caveat.",
        "No external gene-level ground truth; DE validated internally + biologically.",
        "Everything reproduces: `python run_pipeline.py --all`; 37 offline intent tests.",
    ])

    out = SLIDES_DIR / "GI_nanoplastic.pptx"
    prs.save(out)
    print(f"saved slide deck -> {out}  ({len(list(prs.slides))} slides)")
    return out


if __name__ == "__main__":
    build()
