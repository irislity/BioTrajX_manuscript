#!/usr/bin/env python3
"""
audit_fontsizes.py

Measures the EFFECTIVE font size of text inside each panel of the final
assembled supplementary/main figures, accounting for the scale-down (or
scale-up) applied when assemble_figure.py rasterises each source panel PDF
into its slot on the page.

Panel content in the assembled PDF is a raster image (baked in by
assemble_figure.py), so font sizes can't be read back from the final PDF
directly. Instead this script:
  1. Reuses assemble_figure.py's own layout code to get each panel's
     rendered (width, height) in the final page, in px.
  2. Reads each panel's *source* PDF (pre-assembly) native size + all text
     spans (with native font size) via pymupdf.
  3. scale = min(rendered_w/native_w, rendered_h/native_h)   [preserveAspectRatio="meet"]
  4. effective_size = native_size * scale

Reports, per panel: min/max effective size and any text below 7pt.

Usage:
  python3 audit_fontsizes.py
"""
import sys
from pathlib import Path

import pymupdf
sys.path.insert(0, str(Path(__file__).parent))
import assemble_figure as af

ASSEMBLY = Path(__file__).parent


def panel_text_sizes(pdf_path: Path):
    doc = pymupdf.open(str(pdf_path))
    page = doc[0]
    spans = []
    d = page.get_text("dict")
    for block in d.get("blocks", []):
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                txt = span["text"].strip()
                if not txt:
                    continue
                spans.append((round(span["size"], 2), "Bold" in span.get("font", "") or "bold" in span.get("font", "").lower(), txt))
    return spans


def audit_layout(name: str, layout, page_w: int, page_h: int, margin: int, gap: int):
    available_w = page_w - 2 * margin
    all_widths = [af.panel_widths(r, available_w, gap) for r in layout]
    all_heights = [af.row_height(r, ws) for r, ws in zip(layout, all_widths)]

    print(f"\n{'='*80}\n{name}  (page {page_w}x{page_h}px)\n{'='*80}")
    for row, widths, rh in zip(layout, all_widths, all_heights):
        for panel, pw in zip(row.panels, widths):
            label = panel.label or "?"
            native_w, native_h = af.image_size(panel.file)  # px at pdf default dpi conversion inside
            # image_size for pdf returns px at dpi=150 scale in af._pdf_page_size_px; we need native PT size instead
            doc = pymupdf.open(str(panel.file))
            rect = doc[0].rect
            native_w_pt, native_h_pt = rect.width, rect.height

            rendered_w_pt = pw * 0.75
            rendered_h_pt = rh * 0.75
            scale = min(rendered_w_pt / native_w_pt, rendered_h_pt / native_h_pt)

            spans = panel_text_sizes(panel.file)
            if not spans:
                print(f"  [{label}] {panel.file.name}: no text found, scale={scale:.3f}")
                continue
            eff = [(round(sz * scale, 2), bold, txt) for sz, bold, txt in spans]
            min_eff = min(e[0] for e in eff)
            max_eff = max(e[0] for e in eff)
            below7 = [e for e in eff if e[0] < 7.0]
            below8 = [e for e in eff if 7.0 <= e[0] < 8.0]
            print(f"  [{label}] {panel.file.name}  scale={scale:.3f}  "
                  f"eff_range=[{min_eff:.2f}, {max_eff:.2f}]pt  n_spans={len(eff)}")
            if below7:
                uniq = sorted(set((s, t[:30]) for s, b, t in below7))
                print(f"      BELOW 7pt ({len(below7)} spans): {uniq[:6]}{' ...' if len(uniq)>6 else ''}")


def audit_glob(name, files_glob, cols, page_w=794, page_h=1123, margin=30, gap=14):
    files = sorted(ASSEMBLY.glob(files_glob), key=af.natural_sort_key)
    if not files:
        print(f"\n{name}: NO FILES matched {files_glob}")
        return
    layout = af.auto_grid(files, cols, available_w=page_w - 2 * margin, gap=gap)
    audit_layout(name, layout, page_w, page_h, margin, gap)


def audit_multiglob(name, globs, cols, page_w=794, page_h=1123, margin=30, gap=14):
    files = []
    for g in globs:
        files.extend(sorted(ASSEMBLY.glob(g), key=af.natural_sort_key))
    if not files:
        print(f"\n{name}: NO FILES matched {globs}")
        return
    layout = af.auto_grid(files, cols, available_w=page_w - 2 * margin, gap=gap)
    audit_layout(name, layout, page_w, page_h, margin, gap)


def audit_json_layout(name, json_path, page_w=794, page_h=1123, margin=30, gap=14, **overrides):
    import json
    page_w = overrides.get("page_w", page_w)
    page_h = overrides.get("page_h", page_h)
    margin = overrides.get("margin", margin)
    gap = overrides.get("gap", gap)
    raw = json.loads((ASSEMBLY / json_path).read_text())
    rows_raw = raw.get("rows", raw) if isinstance(raw, dict) else raw
    layout = af.parse_layout_rows(rows_raw, ASSEMBLY)
    audit_layout(name, layout, page_w, page_h, margin, gap)


if __name__ == "__main__":
    # All figures assemble onto true A4 portrait (assemble_figure.py's default
    # PAGE_W=794, PAGE_H=1123 -- see assemble_supplementary.sh), so every call
    # below uses those defaults; only S1/S11 need explicit page_h overrides
    # because assemble_supplementary.sh hand-tunes their row heights to fit.
    # Mirrors the exact glob/cols/layout args each figure is assembled with
    # in assemble_supplementary.sh -- keep these in sync with that file.
    audit_json_layout("S1", "s1_layout.json")
    audit_glob("S2", "plots/[0-9][0-9]*_a_umap.pdf", cols=4)
    audit_multiglob("S3", ["plots/[0-9][0-9]*_b_gene_trends.pdf",
                            "plots/[0-9][0-9]*_b_gene_trends_AB.pdf"], cols=3)
    audit_multiglob("S4", ["plots/[0-9][0-9]*_c_D_metric.pdf",
                            "plots/[0-9][0-9]*_c_D_metric_AB.pdf"], cols=4)
    audit_multiglob("S5", ["plots/[0-9][0-9]*_d_O_metric.pdf",
                            "plots/[0-9][0-9]*_d_O_metric_AB.pdf"], cols=4)
    audit_multiglob("S6", ["plots/[0-9][0-9]*_e_E_metric.pdf",
                            "plots/[0-9][0-9]*_e_E_metric_AB.pdf"], cols=4)
    audit_multiglob("S7", ["plots/[0-9][0-9]*_f_ground_truth.pdf",
                            "plots/[0-9][0-9]*_f_ground_truth_AB.pdf"], cols=2)
    audit_json_layout("S8", "s8_layout.json")
    audit_json_layout("S9", "s9_layout.json")
    audit_json_layout("S10", "s10_layout.json")
    audit_json_layout("S11", "s11_layout.json")
