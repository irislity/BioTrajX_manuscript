#!/usr/bin/env python3
"""
make_figure_1.py

Assembles Figure 1 for the BioTrajX manuscript.

Panels:
  A — BioTrajX workflow overview     (manuscript/figures/main/panel_A.svg)
  B — D / O / E metric schematic     (man/figures/DOE.jpg)
  C — DOE score heatmap              (manuscript/figures/real/S11/S11_e_doe_heatmap.pdf)
  D — Day-of-infection recovery vs DOE (manuscript/figures/real/S8/S8_d_day_corr_vs_doe.pdf)

Layout:
  Row 1: Panel A  (full width)
  Row 2: Panel B  (full width)
  Row 3: Panel C  (left) | Panel D (right)

Output:
  manuscript/figures/main/pdf/Figure1.pdf
  manuscript/figures/main/pdf/Figure1.png

Usage (from repo root):
  python3 manuscript/scripts/main/make_figure_1.py
"""

import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image as _PIL
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "Arial"
# Make $...$ mathtext (used for the D_comp/E_comp subscripts) render in
# Arial too, instead of matplotlib's default Computer-Modern math font.
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.rm"] = "Arial"
matplotlib.rcParams["mathtext.it"] = "Arial:italic"
matplotlib.rcParams["mathtext.bf"] = "Arial:bold"
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.colors

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO    = Path(__file__).resolve().parents[3]
MAIN    = REPO / "manuscript" / "figures" / "main"
S8      = REPO / "manuscript" / "figures" / "real" / "S8"
S11     = REPO / "manuscript" / "figures" / "real" / "S11"

PANEL_A     = MAIN / "panel_A.pdf"
DOE_JPG     = MAIN / "DOE.jpg"
PANEL_C     = S11  / "S11_e_doe_heatmap.pdf"
PANEL_D     = S8   / "S8_d_day_corr_vs_doe.pdf"

# ── Panel A amber "Composite DOE scores" box, reused for Panel B's box ──
DOE_BOX_BORDER = "#f4c542"
DOE_BOX_FILL   = "#fff2cc"
DOE_BOX_TEXT   = "#000000"

PANEL_B_TITLES  = ["Directionality (D)", "Order Consistency (O)", "Endpoint Validity (E)"]
PANEL_B_CAPTIONS = [["Spearman", "correlation"], ["Monotonic", "Trend"], ["Precision@k"]]
PANEL_B_NOTE = "DOE = mean(D$_{comp}$, O, E$_{comp}$)"

PANEL_C_TITLE    = "BioTrajX provides a quantitative criterion for trajectory root selection"
PANEL_C_SUBTITLE = "Monocle3 roots spanning major CD8+ T-cell states"
PANEL_D_TITLE    = "Higher DOE scores correspond to better recovery of infection timing"
PANEL_D_SUBTITLE = "Independent LCMV time-course validation (GSE131847)"

# Unified panel-letter size (A–D identical) and Panel C/D narrative-title size
# (kept at/under Panel B's metric-title size so no panel dominates the others).
LETTER_FS      = 13.0
B_TITLE_FS     = 15.0
CD_TITLE_FS    = 12.5
CD_SUB_FS      = 10.0

OUT_DIR = MAIN / "pdf"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PNG = str(OUT_DIR / "Figure1.png")
OUT_PDF = str(OUT_DIR / "Figure1.pdf")

DPI = 300


# ---------------------------------------------------------------------------
# Panel regeneration
# ---------------------------------------------------------------------------

def _check_panels() -> None:
    missing = []
    for path, label in [
        (PANEL_A, "A (workflow SVG)"),
        (DOE_JPG, "B (DOE schematic image)"),
        (PANEL_C, "C (DOE heatmap)"),
        (PANEL_D, "D (day-corr-vs-DOE)"),
    ]:
        if not path.exists():
            missing.append(f"  Panel {label}: {path}")
    if missing:
        print("Missing panels:", file=sys.stderr)
        for m in missing:
            print(m, file=sys.stderr)
        sys.exit(1)




# ---------------------------------------------------------------------------
# Conversion helpers
# ---------------------------------------------------------------------------

def _trim_edges(img, edges="top", threshold: int = 245,
                min_px: int = 10, top_pad: int = 0, bottom_pad: int = 4):
    """Crop near-white margins from the specified edges ('top', 'bottom', or 'both')."""
    import numpy as np
    arr = np.array(img.convert("RGB"))
    row_nz = np.sum(np.any(arr < threshold, axis=2), axis=1)
    content_rows = np.where(row_nz > min_px)[0]
    if not content_rows.size:
        return img
    rmin = max(0, content_rows[0] - top_pad)                          if edges in ("top",    "both") else 0
    rmax = min(arr.shape[0], content_rows[-1] + bottom_pad + 1)       if edges in ("bottom", "both") else img.height
    return img.crop((0, rmin, img.width, rmax))


def _to_png(src: Path, out_dir: str, *, trim: str = "",
            dpi: int = DPI, trim_threshold: int = 245,
            top_pad: int = 0, bottom_pad: int = 4) -> Path:
    import pymupdf
    from PIL import Image as _PIL
    import cairosvg
    out = Path(out_dir) / f"{src.stem}.png"
    suffix = src.suffix.lower()
    if suffix in (".png", ".jpg", ".jpeg"):
        img = _PIL.open(str(src))
        if trim:
            img = _trim_edges(img, trim, threshold=trim_threshold,
                              top_pad=top_pad, bottom_pad=bottom_pad)
        img.save(str(out))
        return out
    if suffix == ".svg":
        cairosvg.svg2png(url=str(src), write_to=str(out), dpi=dpi)
    elif suffix == ".pdf":
        doc = pymupdf.open(str(src))
        mat = pymupdf.Matrix(dpi / 72, dpi / 72)
        pix = doc[0].get_pixmap(matrix=mat, alpha=False)
        pix.save(str(out))
    else:
        raise ValueError(f"Unsupported format: {suffix}")
    if trim:
        img = _PIL.open(str(out))
        img = _trim_edges(img, trim, threshold=trim_threshold,
                          top_pad=top_pad, bottom_pad=bottom_pad)
        img.save(str(out))
    return out


def _aspect(path: Path) -> float:
    import pymupdf
    from PIL import Image as _PIL
    suffix = path.suffix.lower()
    if suffix == ".pdf":
        doc = pymupdf.open(str(path))
        r = doc[0].rect
        return r.width / r.height
    else:
        img = _PIL.open(str(path))
        return img.size[0] / img.size[1]


def _measure_width_in(text: str, fontsize: float,
                       fontweight: str = "normal", fontstyle: str = "normal") -> float:
    """Rendered width of `text` in inches, using the actual matplotlib font/renderer."""
    fig_m = plt.figure()
    fig_m.canvas.draw()
    renderer = fig_m.canvas.get_renderer()
    t = fig_m.text(0, 0, text, fontsize=fontsize, fontweight=fontweight, fontstyle=fontstyle)
    fig_m.canvas.draw()
    w = t.get_window_extent(renderer=renderer).width / fig_m.dpi
    plt.close(fig_m)
    return w


def _wrap_to_width(text: str, max_width_in: float, fontsize: float,
                    fontweight: str = "normal", fontstyle: str = "normal",
                    first_line_width_in: float | None = None) -> list[str]:
    """Greedy word-wrap `text` so each line's rendered width fits max_width_in,
    measuring with the actual matplotlib font/renderer used at save time.
    If `first_line_width_in` is given, only the first line uses that budget
    (e.g. to leave room for a letter/marker preceding the paragraph)."""
    words = text.split()
    lines: list[str] = []
    current = ""
    fig_m = plt.figure()
    fig_m.canvas.draw()
    renderer = fig_m.canvas.get_renderer()

    def width_in(s: str) -> float:
        t = fig_m.text(0, 0, s, fontsize=fontsize, fontweight=fontweight, fontstyle=fontstyle)
        fig_m.canvas.draw()
        w = t.get_window_extent(renderer=renderer).width / fig_m.dpi
        t.remove()
        return w

    for word in words:
        trial = (current + " " + word).strip()
        limit = (first_line_width_in if (not lines and first_line_width_in is not None)
                 else max_width_in)
        if current == "" or width_in(trial) <= limit:
            current = trial
        else:
            lines.append(current)
            current = word
    if current:
        lines.append(current)
    plt.close(fig_m)
    return lines


def _prepare_doe_crops():
    """Crop the three D/O/E illustrations out of DOE.jpg, with their baked
    titles removed (redrawn as vector Arial) and the baked "Spearman
    correlation" caption erased (also redrawn as vector Arial). The
    "Pseudotime" label baked onto the O-panel arrow is kept as-is, per
    design: it stays attached to the arrow graphic."""
    arr = np.array(_PIL.open(str(DOE_JPG)).convert("RGB")).copy()
    arr[200:352, 322:535] = 255  # erase baked "Spearman correlation"
    # Erase the baked "Endpoint Validity (E)" title (its descenders dip as
    # low as y~78, into the same band the top inset circle starts in) —
    # everything except the circle's own x-span, which is left untouched.
    arr[0:96, 1221:1420] = 255
    arr[0:96, 1545:1832] = 255

    def crop_padded(bbox, pad=6):
        x0, x1, y0, y1 = bbox
        x0 = max(0, x0 - pad); y0 = max(0, y0 - pad)
        x1 = min(arr.shape[1], x1 + pad); y1 = min(arr.shape[0], y1 + pad)
        return arr[y0:y1, x0:x1]

    D_img = crop_padded((36, 389, 106, 430))
    O_img = crop_padded((695, 1095, 149, 402))
    # y0=79 (not 96) — the top inset circle's arc starts at y~81, above where
    # the old crop cut it off.
    E_img = crop_padded((1276, 1743, 79, 463))
    return [D_img, O_img, E_img]


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def assemble(panel_a: Path, panel_c: Path, panel_d: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        png = {
            "A": _to_png(panel_a, tmp, trim="both",
                         trim_threshold=253, top_pad=0, bottom_pad=6),
            "C": _to_png(panel_c, tmp),
            "D": _to_png(panel_d, tmp),
        }
        doe_crops = _prepare_doe_crops()   # [D_img, O_img, E_img] numpy arrays

        # ── Figure sizing ──────────────────────────────────────────────
        FIG_W    = 8.3
        ML, MR   = 0.015, 0.015
        MT, MB   = 0.010, 0.010
        ROW_GAP  = 0.22   # inches: gap between A↔B rows
        ROW_GAP2 = 0.20   # inches: gap between B and C/D rows
        COL_GAP  = 0.08   # inches: gap between C and D

        # Panel B (D/O/E schematic) layout
        B_COL_GAP      = 0.10   # inches: gap between the D/O/E sub-columns
        B_TITLE_LH     = B_TITLE_FS * 1.25 / 72
        B_TITLE_ROW_H  = 0.08 + B_TITLE_LH + 0.06
        B_IMG_ROW_H    = 1.35   # inches: D/O illustration height
        B_GAP_IMG_CAP  = 0.06
        B_CAP_FS       = 13.0
        B_CAP_LH       = B_CAP_FS * 1.25 / 72
        B_CAPTION_ROW_H = 2 * B_CAP_LH + 0.04
        # Shared (image + caption) content height at the deepest caption (2
        # lines, D/O). A column with a shorter caption (E's 1-line
        # "Precision@k") gets a taller illustration instead — same total
        # height, more of it spent on the image — so E renders larger while
        # every column still starts at the same top row.
        B_MAX_CONTENT_H = B_IMG_ROW_H + B_GAP_IMG_CAP + B_CAPTION_ROW_H
        B_ILLUS_BAND_H  = B_TITLE_ROW_H + 0.05 + B_MAX_CONTENT_H
        B_GAP_TO_BOX   = 0.10
        B_BOX_H        = 0.36
        h_b = B_ILLUS_BAND_H + B_GAP_TO_BOX + B_BOX_H

        # Panel C/D headers: unified letter (LETTER_FS, matches A/B) followed
        # by a narrative title (CD_TITLE_FS, kept <= Panel B's metric titles)
        CD_TITLE_LH   = CD_TITLE_FS * 1.25 / 72     # line height, inches
        CD_SUB_LH     = CD_SUB_FS * 1.3 / 72        # line height, inches
        CD_TOP_PAD       = 0.03   # inches: above first title line
        CD_TITLE_SUB_GAP = 0.05   # inches: between title block and subtitle
        CD_SUB_IMG_GAP   = 0.08   # inches: between subtitle and panel image
        CD_LETTER_TITLE_GAP = 0.05   # inches: between the "C"/"D" letter and its title

        cw = (1.0 - ML - MR) * FIG_W   # content width in inches

        r_a = _aspect(png["A"])
        r_c = _aspect(png["C"])
        r_d = _aspect(png["D"])

        h_a    = cw / r_a
        h_cd_img = (cw - COL_GAP) / (r_c + r_d)
        fc     = h_cd_img * r_c          # width of C in inches
        fd     = h_cd_img * r_d          # width of D in inches

        # Wrap the C/D headers to their own column widths, using the actual
        # font/renderer so line counts (and thus header height) are exact.
        # The first line of the title budgets extra room for the letter.
        letter_w_c = _measure_width_in("C", LETTER_FS, fontweight="bold")
        letter_w_d = _measure_width_in("D", LETTER_FS, fontweight="bold")
        c_title_lines = _wrap_to_width(PANEL_C_TITLE, fc, CD_TITLE_FS, fontweight="bold",
                                        first_line_width_in=fc - letter_w_c - CD_LETTER_TITLE_GAP)
        d_title_lines = _wrap_to_width(PANEL_D_TITLE, fd, CD_TITLE_FS, fontweight="bold",
                                        first_line_width_in=fd - letter_w_d - CD_LETTER_TITLE_GAP)
        c_sub_lines   = _wrap_to_width(PANEL_C_SUBTITLE, fc, CD_SUB_FS, fontstyle="italic")
        d_sub_lines   = _wrap_to_width(PANEL_D_SUBTITLE, fd, CD_SUB_FS, fontstyle="italic")
        n_title_lines = max(len(c_title_lines), len(d_title_lines))
        n_sub_lines   = max(len(c_sub_lines), len(d_sub_lines))

        cd_header_h = (CD_TOP_PAD + n_title_lines * CD_TITLE_LH + CD_TITLE_SUB_GAP
                       + n_sub_lines * CD_SUB_LH + CD_SUB_IMG_GAP)
        h_cd = cd_header_h + h_cd_img   # total row height: header + image

        content_h = h_a + ROW_GAP + h_b + ROW_GAP2 + h_cd
        FIG_H     = content_h / (1.0 - MT - MB)

        def fh(x): return x / FIG_H
        def fw(x): return x / FIG_W

        # Y positions (bottom of each panel, from figure bottom)
        y_a  = 1.0 - MT - fh(h_a)
        y_b  = y_a  - fh(ROW_GAP)  - fh(h_b)
        y_cd = y_b  - fh(ROW_GAP2) - fh(h_cd)

        fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI)

        # ── Panel A ──────────────────────────────────────────────────
        ax_a = fig.add_axes([ML, y_a, 1 - ML - MR, fh(h_a)])
        ax_a.imshow(mpimg.imread(str(png["A"])))
        ax_a.axis("off")
        fig.text(ML + 0.006, y_a + fh(h_a) + 0.006, "A",
                 fontsize=LETTER_FS, fontweight="bold", ha="left", va="bottom", color="black")

        # ── Panel B: D / O / E schematic, rebuilt from cleaned DOE.jpg crops ──
        img_cd_frac = h_cd_img / h_cd   # fraction of C/D height occupied by the plot image
        sub_col_w = (cw - 2 * B_COL_GAP) / 3
        y_illus = y_b + fh(B_GAP_TO_BOX + B_BOX_H)

        fig.text(ML + 0.006, y_illus + fh(B_ILLUS_BAND_H) + 0.006, "B",
                 fontsize=LETTER_FS, fontweight="bold", ha="left", va="bottom", color="black")

        # All three captions share the same (image + caption) height budget, so
        # a shorter caption (E's single-line "Precision@k") converts directly
        # into a taller illustration — this is what makes E render larger
        # while keeping every column's total content height identical.
        cap_max_lines = max(len(c) for c in PANEL_B_CAPTIONS)
        for i, (title, im, cap) in enumerate(zip(PANEL_B_TITLES, doe_crops, PANEL_B_CAPTIONS)):
            x0 = ML + i * (sub_col_w + B_COL_GAP)
            ax = fig.add_axes([fw(x0), y_illus, fw(sub_col_w), fh(B_ILLUS_BAND_H)])
            ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
            ax.text(0.5, 1.0 - (0.08 / B_ILLUS_BAND_H), title, ha="center", va="top",
                    fontsize=B_TITLE_FS, fontweight="bold", color="black")
            r = im.shape[1] / im.shape[0]
            this_img_h = B_IMG_ROW_H + (cap_max_lines - len(cap)) * B_CAP_LH
            img_h_frac = this_img_h / B_ILLUS_BAND_H
            img_w_frac = min((this_img_h * r) / sub_col_w, 1.0)
            img_top_frac = 1 - (B_TITLE_ROW_H + 0.05) / B_ILLUS_BAND_H
            img_x0 = (1 - img_w_frac) / 2
            img_ax = ax.inset_axes([img_x0, img_top_frac - img_h_frac, img_w_frac, img_h_frac])
            img_ax.imshow(im)
            img_ax.axis("off")
            cap_top = img_top_frac - img_h_frac - (B_GAP_IMG_CAP / B_ILLUS_BAND_H)
            line_h_frac = B_CAP_LH / B_ILLUS_BAND_H
            y = cap_top
            for line in cap:
                ax.text(0.5, y, line, ha="center", va="top",
                        fontsize=B_CAP_FS, color="black")
                y -= line_h_frac

        box_ax = fig.add_axes([ML, y_b, fw(cw), fh(B_BOX_H)])
        box_ax.set_xlim(0, 1); box_ax.set_ylim(0, 1); box_ax.axis("off")
        note_w = 0.40
        note_x0 = (1 - note_w) / 2
        box_ax.add_patch(mpatches.FancyBboxPatch(
            (note_x0, 0.06), note_w, 0.88,
            boxstyle="round,pad=0.02,rounding_size=0.05",
            # Fill at 0.4 opacity (border fully opaque) — matches the exact
            # opacity/border split on Panel A's "Composite DOE scores" box.
            facecolor=matplotlib.colors.to_rgba(DOE_BOX_FILL, alpha=0.4),
            edgecolor=DOE_BOX_BORDER, linewidth=1.4,
            transform=box_ax.transAxes, clip_on=False,
        ))
        box_ax.text(0.5, 0.5, PANEL_B_NOTE, ha="center", va="center",
                    fontsize=12.5, color=DOE_BOX_TEXT, fontweight="bold",
                    transform=box_ax.transAxes)

        # ── Panels C and D ───────────────────────────────────────────
        def _place_header(ax, letter, letter_w, title_lines, sub_lines, h_total, col_w):
            y = 1.0 - CD_TOP_PAD / h_total
            ax.text(0.0, y, letter, ha="left", va="top",
                    fontsize=LETTER_FS, fontweight="bold", color="black")
            if title_lines:
                x_title = (letter_w + CD_LETTER_TITLE_GAP) / col_w
                ax.text(x_title, y, title_lines[0], ha="left", va="top",
                        fontsize=CD_TITLE_FS, fontweight="bold", color="black")
                y -= CD_TITLE_LH / h_total
                for line in title_lines[1:]:
                    ax.text(0.0, y, line, ha="left", va="top",
                            fontsize=CD_TITLE_FS, fontweight="bold", color="black")
                    y -= CD_TITLE_LH / h_total
            y -= CD_TITLE_SUB_GAP / h_total
            for line in sub_lines:
                ax.text(0.0, y, line, ha="left", va="top",
                        fontsize=CD_SUB_FS, fontstyle="italic", color="#555555")
                y -= CD_SUB_LH / h_total

        cd_specs = [
            ("C", fc, ML,                    c_title_lines, c_sub_lines, letter_w_c),
            ("D", fd, ML + fw(fc + COL_GAP), d_title_lines, d_sub_lines, letter_w_d),
        ]
        for label, col_w_in, x0, title_lines, sub_lines, letter_w in cd_specs:
            ax = fig.add_axes([x0, y_cd, fw(col_w_in), fh(h_cd)])
            ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
            img_ax = ax.inset_axes([0, 0, 1, img_cd_frac])
            img_ax.imshow(mpimg.imread(str(png[label])))
            img_ax.axis("off")
            for spine in img_ax.spines.values():
                spine.set_visible(True)
                spine.set_edgecolor("#CCCCCC")
                spine.set_linewidth(0.8)
            _place_header(ax, label, letter_w, title_lines, sub_lines, h_cd, col_w_in)

        for out in (OUT_PNG, OUT_PDF):
            fig.savefig(out, dpi=DPI, bbox_inches="tight",
                        facecolor="white", pad_inches=0.03)
            print(f"Saved: {out}")

        plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    _check_panels()
    print("Assembling Figure 1 ...")
    assemble(PANEL_A, PANEL_C, PANEL_D)
    print("Done.")
