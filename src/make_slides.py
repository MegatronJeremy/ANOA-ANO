"""
Build the results slide deck (deliverable) from the figures the pipeline saves.

Reproducible: re-run after `--all` to rebuild results/slides/GI_nanoplastic.pptx
from whatever is currently in results/figures/. Missing figures are skipped with
a note rather than crashing, so a partial run still yields a deck.

    python -m src.make_slides          # or:  .\run.ps1 slides
"""
from pathlib import Path

from pptx import Presentation
from pptx.util import Inches, Pt

from . import config as cfg

SLIDES_DIR = cfg.PROJECT_ROOT / "results" / "slides"
SLIDES_DIR.mkdir(parents=True, exist_ok=True)
FIG = cfg.FIG_DIR


def _title_slide(prs, title, subtitle):
    s = prs.slides.add_slide(prs.slide_layouts[0])
    s.shapes.title.text = title
    s.placeholders[1].text = subtitle


def _bullets_slide(prs, title, bullets):
    s = prs.slides.add_slide(prs.slide_layouts[1])
    s.shapes.title.text = title
    tf = s.placeholders[1].text_frame
    tf.clear()
    for i, b in enumerate(bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = b
        p.font.size = Pt(18)


def _figure_slide(prs, title, fig_names, caption=""):
    """One slide, title + up to two figures side by side (skips any missing)."""
    paths = [FIG / f for f in fig_names if (FIG / f).exists()]
    s = prs.slides.add_slide(prs.slide_layouts[5])
    s.shapes.title.text = title
    if not paths:
        tb = s.shapes.add_textbox(Inches(1), Inches(2.5), Inches(8), Inches(1)).text_frame
        tb.text = "(figure not generated yet -- run the pipeline)"
        return
    if len(paths) == 1:
        s.shapes.add_picture(str(paths[0]), Inches(1.2), Inches(1.5), height=Inches(5.2))
    else:
        s.shapes.add_picture(str(paths[0]), Inches(0.3), Inches(1.6), height=Inches(4.6))
        s.shapes.add_picture(str(paths[1]), Inches(5.1), Inches(1.6), height=Inches(4.6))
    if caption:
        tb = s.shapes.add_textbox(Inches(0.5), Inches(6.7), Inches(9), Inches(0.6)).text_frame
        tb.text = caption
        tb.paragraphs[0].font.size = Pt(12)


def build():
    prs = Presentation()
    prs.slide_width = Inches(10)
    prs.slide_height = Inches(7.5)

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
    _figure_slide(prs, "QC & preprocessing", ["01_qc_violin_after.png", "01_qc_scatter_counts_mito.png"],
                  "Thresholds: max_genes~p99, max_mito between p95-p99; justified on the real distribution.")
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
    _figure_slide(prs, "Bonus analyses", ["07_module_scores.png", "07_mixture_additivity.png"],
                  "Stress/inflammation module scores; mixture additivity (observed vs 40+200 expectation).")
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
