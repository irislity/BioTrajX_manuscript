#!/usr/bin/env python3
"""
General-purpose manuscript figure assembler.

Arranges any collection of PNG/PDF images into a multi-panel figure PDF
(A4 portrait by default) with auto-lettered panels and titles derived
from filenames.  Layout can be uniform auto-grid, or fully custom via a
Python block at the top of this file, or an external JSON file.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
QUICK START
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  # Auto-grid: every PNG in ./plots/, 3 columns, A4 portrait
  python assemble_figure.py --input-dir plots --cols 3 --out figure.pdf

  # Specific file list
  python assemble_figure.py plots/fig_a.png plots/fig_b.png plots/fig_c.png

  # Pattern glob
  python assemble_figure.py --glob "plots/*_heatmap.png" --cols 2

  # Custom layout (JSON)
  python assemble_figure.py --layout layout.json

  # A4 landscape
  python assemble_figure.py --glob "plots/*.png" --landscape

  # Custom page size (in SVG px; 1 inch = 96 px)
  python assemble_figure.py --glob "plots/*.png" --width 1587 --height 2245   # A3 portrait

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
MANUAL LAYOUT FORMAT
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Set MANUAL_LAYOUT below to a list of rows.  Each row is a list of panel dicts:

  {
    "file":   "path/to/image.png",   # relative to this script or absolute
    "weight": 1.0,                   # relative width within the row (default 1.0)
    "label":  "A",                   # override auto-letter  (optional)
    "title":  "My Panel Title",      # override auto-title   (optional)
  }

Optionally wrap a row in {"panels": [...], "row_height": 200} to fix its height.

EXAMPLE — row 1 has one wide panel (weight 2) + two normal panels:
  MANUAL_LAYOUT = [
      [
          {"file": "plots/wide_panel.png",  "weight": 2.0},
          {"file": "plots/panel_b.png"},
          {"file": "plots/panel_c.png"},
      ],
      [
          {"file": "plots/panel_d.png"},
          {"file": "plots/panel_e.png"},
          {"file": "plots/panel_f.png"},
          {"file": "plots/panel_g.png"},
      ],
  ]

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
JSON LAYOUT FORMAT  (--layout file.json)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Same structure as MANUAL_LAYOUT, serialised as JSON.
Top level can be a list of rows, or {"rows": [...], "figure_title": "Fig 1"}.

Dependencies: pip install svgutils cairosvg Pillow lxml
"""

# ============================================================
#  USER-EDITABLE DEFAULTS
# ============================================================

# A4 portrait at 96 SVG px/inch.  A4 landscape = swap these two.
PAGE_W: int = 794
PAGE_H: int = 1123

MARGIN: int  = 30    # outer margin on all four sides (SVG px)
GAP: int     = 14    # gap between panels (SVG px)
LABEL_H: int = 36    # tall enough for two-line titles (scenario + branch)    # height reserved above each panel for letter + title
FONT_LABEL: int = 15
FONT_TITLE: int = 11

# Optional: hardcode a layout here (bypasses auto-grid and --glob).
# Set to None to use CLI flags instead.
MANUAL_LAYOUT = None

# ============================================================
#  IMPLEMENTATION
# ============================================================

import argparse
import base64
import glob as glob_module
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import cairosvg
import svgutils.transform as sg
from lxml import etree

SVG_NS   = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
PANEL_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class PanelSpec:
    file: Path
    weight: float = 1.0
    label: Optional[str] = None
    title: Optional[str] = None


@dataclass
class RowSpec:
    panels: list[PanelSpec] = field(default_factory=list)
    row_height: Optional[int] = None
    fixed_panel_w: Optional[float] = None  # if set, every panel in this row uses this width


# ---------------------------------------------------------------------------
# Filename → title
# ---------------------------------------------------------------------------

# Strip tokens from filenames that rarely make useful panel titles.
_STRIP_WORDS = frozenset({"png", "pdf", "svg", "fig", "figure", "plot", "panel"})

# Panel-type suffixes added by visualize_fixtures.R (e.g. "_a_umap", "_c_D_metric").
# Optional trailing branch tag "_AB" or "_AC" is captured in group 2.
_PANEL_SUFFIX_RE = re.compile(
    r"_[a-z]_(umap|gene_trends|[A-Za-z]+_metric|ground_truth|marker_scores|heatmap)(?:_(AB|AC))?$",
    re.IGNORECASE,
)

# Map branch codes to arrow notation
_BRANCH_LABEL = {"AB": "A\u2192B", "AC": "A\u2192C"}


def auto_title(path: Path) -> str:
    """Derive a human-readable title from a filename.

    Logic:
      1. Take the stem (no extension).
      2. Strip a leading numeric prefix (e.g. '01_', '02-').
      3. Strip trailing panel-type suffixes (e.g. '_a_umap', '_c_D_metric').
         If a branch tag (_AB / _AC) is present it is captured and appended
         as a second title line: "Scenario Name\\nA→B".
      4. Split on underscores / hyphens / camelCase boundaries.
      5. Drop noise tokens (see _STRIP_WORDS).
      6. Title-case the remaining words.
    """
    stem = path.stem
    # Remove leading numeric prefix: "01_" or "01-" or "01 "
    stem = re.sub(r"^\d+[\s_\-]+", "", stem)
    # Remove trailing panel-type suffix, capturing any branch tag
    m = _PANEL_SUFFIX_RE.search(stem)
    branch_line = _BRANCH_LABEL.get((m.group(2) or "").upper()) if m else None
    stem = _PANEL_SUFFIX_RE.sub("", stem)
    # Split on underscores, hyphens, dots, and camelCase boundaries
    tokens = re.sub(r"([a-z])([A-Z])", r"\1 \2", stem)
    tokens = re.split(r"[\s_\-\.]+", tokens)
    words = [t for t in tokens if t.lower() not in _STRIP_WORDS and t]
    title = " ".join(w.capitalize() for w in words) if words else stem
    # Restore known all-caps acronyms that title-case mangles
    title = re.sub(r"\bZinb\b", "ZINB", title, flags=re.IGNORECASE)
    # Long titles (3+ words) wrap onto two lines, split at the midpoint, so
    # they don't overflow a panel's width in dense multi-panel grids.
    title_words = title.split(" ")
    if len(title_words) >= 3:
        mid = (len(title_words) + 1) // 2
        title = " ".join(title_words[:mid]) + "\n" + " ".join(title_words[mid:])
    if branch_line:
        title = f"{title}\n{branch_line}"
    return title


# ---------------------------------------------------------------------------
# Image helpers
# ---------------------------------------------------------------------------

def _pdf_to_png_bytes(path: Path, dpi: int) -> bytes:
    """Rasterise the first page of a PDF to PNG bytes using pymupdf."""
    import pymupdf  # pip install pymupdf
    doc  = pymupdf.open(str(path))
    page = doc[0]
    mat  = pymupdf.Matrix(dpi / 72, dpi / 72)   # 72 pt/in → dpi px/in
    pix  = page.get_pixmap(matrix=mat, alpha=False)
    return pix.tobytes("png")


def _pdf_page_size_px(path: Path, dpi: int = 150) -> tuple[int, int]:
    """Return (width, height) in pixels for the first page of a PDF at dpi."""
    import pymupdf
    doc  = pymupdf.open(str(path))
    rect = doc[0].rect          # points (1 pt = 1/72 in)
    w    = int(rect.width  * dpi / 72)
    h    = int(rect.height * dpi / 72)
    return w, h


def image_size(path: Path) -> tuple[int, int]:
    """Return (width, height) pixels; handles raster images and PDFs."""
    suffix = path.suffix.lower()
    if suffix == ".pdf":
        return _pdf_page_size_px(path)
    try:
        from PIL import Image
        with Image.open(path) as img:
            return img.size
    except Exception:
        return 800, 600


def aspect_ratio(path: Path) -> float:
    w, h = image_size(path)
    return h / w


def load_png_bytes(path: Path, raster_dpi: int = 150) -> bytes:
    """Return PNG bytes for any supported input (raster image or PDF)."""
    suffix = path.suffix.lower()
    if suffix == ".pdf":
        return _pdf_to_png_bytes(path, raster_dpi)
    if suffix in (".png", ".jpg", ".jpeg", ".tiff", ".bmp"):
        return path.read_bytes()
    raise ValueError(f"Unsupported image format: {suffix!r}  ({path.name})")


def b64_uri(raw: bytes, mime: str = "image/png") -> str:
    return f"data:{mime};base64," + base64.b64encode(raw).decode()


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

SUPPORTED_EXTS = {".png", ".jpg", ".jpeg", ".tiff", ".bmp", ".pdf"}


def natural_sort_key(p: Path) -> tuple:
    parts = re.split(r"(\d+)", p.name)
    return tuple(int(x) if x.isdigit() else x.lower() for x in parts)


def collect_files(paths: list[str], pattern: Optional[str], input_dir: Optional[Path]) -> list[Path]:
    """Gather files from explicit paths, a glob pattern, or a directory."""
    found: list[Path] = []

    if paths:
        for p in paths:
            fp = Path(p)
            if not fp.exists():
                sys.exit(f"ERROR: file not found: {fp}")
            found.append(fp)
        return found

    if pattern:
        matched = sorted(glob_module.glob(pattern, recursive=True), key=lambda s: natural_sort_key(Path(s)))
        if not matched:
            sys.exit(f"ERROR: no files matched glob pattern: {pattern!r}")
        return [Path(m) for m in matched]

    if input_dir:
        all_files = sorted(input_dir.iterdir(), key=natural_sort_key)
        found = [f for f in all_files if f.is_file() and f.suffix.lower() in SUPPORTED_EXTS]
        if not found:
            sys.exit(f"ERROR: no supported image files in {input_dir}")
        return found

    sys.exit("ERROR: supply positional file paths, --glob, or --input-dir.")


def resolve_file(spec: str, base_dir: Path) -> Path:
    """Resolve a bare filename, relative path, or glob (must match exactly 1 file)."""
    p = Path(spec)
    if p.is_absolute():
        if p.exists():
            return p
        sys.exit(f"ERROR: file not found: {p}")
    # Try relative to base_dir, then CWD
    for root in (base_dir, Path.cwd()):
        candidate = root / p
        if candidate.exists():
            return candidate
        # Treat as glob
        matches = sorted(root.glob(spec), key=natural_sort_key)
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            sys.exit(f"ERROR: glob {spec!r} matched {len(matches)} files; be more specific.")
    sys.exit(f"ERROR: cannot resolve {spec!r} relative to {base_dir} or {Path.cwd()}")


# ---------------------------------------------------------------------------
# Layout resolution
# ---------------------------------------------------------------------------

def auto_grid(files: list[Path], ncols: int,
              available_w: float = 0, gap: float = 0) -> list[RowSpec]:
    # Compute a canonical panel width so every row (including the short last one)
    # uses the same column width and the grid looks uniform.
    canonical_pw: Optional[float] = None
    if available_w > 0 and ncols > 0:
        canonical_pw = (available_w - gap * (ncols - 1)) / ncols

    rows = []
    for i in range(0, len(files), ncols):
        chunk = files[i:i + ncols]
        rows.append(RowSpec(
            panels=[PanelSpec(file=f) for f in chunk],
            fixed_panel_w=canonical_pw,
        ))
    return rows


def parse_layout_rows(raw_rows: list, base_dir: Path) -> list[RowSpec]:
    rows = []
    for raw in raw_rows:
        if isinstance(raw, dict):
            row_height = raw.get("row_height")
            panels_raw = raw.get("panels", [])
        else:
            row_height = None
            panels_raw = raw

        fixed_panel_w = float(raw["fixed_panel_w"]) if isinstance(raw, dict) and raw.get("fixed_panel_w") is not None else None
        panels = []
        for pdict in panels_raw:
            fp = resolve_file(pdict["file"], base_dir)
            panels.append(PanelSpec(
                file=fp,
                weight=float(pdict.get("weight", 1.0)),
                label=pdict.get("label"),
                title=pdict.get("title"),
            ))
        rows.append(RowSpec(panels=panels, row_height=row_height, fixed_panel_w=fixed_panel_w))
    return rows


def resolve_layout(
    args,
    files: list[Path],
    base_dir: Path,
    available_w: float = 0,
) -> tuple[list[RowSpec], Optional[str]]:
    """Return (layout_rows, optional_figure_title)."""

    if MANUAL_LAYOUT is not None:
        return parse_layout_rows(MANUAL_LAYOUT, base_dir), None

    if args.layout:
        raw = json.loads(Path(args.layout).read_text())
        if isinstance(raw, dict):
            rows_raw = raw.get("rows", [])
            fig_title = raw.get("figure_title")
        else:
            rows_raw = raw
            fig_title = None
        return parse_layout_rows(rows_raw, base_dir), fig_title

    return auto_grid(files, args.cols, available_w=available_w, gap=args.gap), None


# ---------------------------------------------------------------------------
# Row geometry
# ---------------------------------------------------------------------------

def panel_widths(row: RowSpec, available_w: float, gap: float) -> list[float]:
    n = len(row.panels)
    if n == 0:
        return []
    # fixed_panel_w: used by auto_grid to keep all panels the same width even in
    # a short last row, so the grid looks uniform.
    if row.fixed_panel_w is not None:
        return [row.fixed_panel_w] * n
    usable = available_w - gap * (n - 1)
    total_weight = sum(p.weight for p in row.panels)
    return [usable * p.weight / total_weight for p in row.panels]


def row_height(row: RowSpec, widths: list[float]) -> float:
    if row.row_height is not None:
        return float(row.row_height)
    if not widths:
        return 100.0
    return max(w * aspect_ratio(p.file) for p, w in zip(row.panels, widths))


# ---------------------------------------------------------------------------
# SVG builder
# ---------------------------------------------------------------------------

def build_svg(
    layout: list[RowSpec],
    *,
    page_w: int,
    page_h: int,
    margin: int,
    gap: int,
    label_h: int,
    font_label: int,
    font_title: int,
    figure_title: Optional[str] = None,
    show_labels: bool = True,
    show_titles: bool = True,
) -> sg.SVGFigure:

    available_w = page_w - 2 * margin

    # Pre-compute geometry
    all_widths  = [panel_widths(r, available_w, gap) for r in layout]
    all_heights = [row_height(r, ws) for r, ws in zip(layout, all_widths)]

    # Extra space at top for an optional figure title
    title_h = (font_label + 12) if figure_title else 0

    content_h = (
        title_h
        + sum(h + (label_h if show_labels or show_titles else 0) for h in all_heights)
        + gap * max(len(layout) - 1, 0)
    )

    if content_h + 2 * margin > page_h:
        print(
            f"WARNING: content ({content_h + 2*margin:.0f} px) taller than page "
            f"({page_h} px). Try --landscape, more --cols, or a smaller margin/gap.",
            file=sys.stderr,
        )

    # Centre grid vertically
    y0 = (page_h - content_h) // 2

    fig = sg.SVGFigure()
    fig.set_size((f"{page_w}px", f"{page_h}px"))
    root = fig.root

    # White background
    etree.SubElement(root, f"{{{SVG_NS}}}rect", attrib={
        "x": "0", "y": "0",
        "width": str(page_w), "height": str(page_h),
        "fill": "white",
    })

    # Optional figure title
    if figure_title:
        ft = etree.SubElement(root, f"{{{SVG_NS}}}text", attrib={
            "x": str(margin),
            "y": str(y0 + font_label + 2),
            "font-size": str(font_label + 2),
            "font-weight": "bold",
            "font-family": "Times New Roman, serif",
            "text-anchor": "start",
            "fill": "#000000",
        })
        ft.text = figure_title

    panel_counter = 0
    effective_label_h = label_h if (show_labels or show_titles) else 0
    cursor_y = y0 + title_h

    for row, widths, rh in zip(layout, all_widths, all_heights):
        cursor_x = margin

        for panel, pw in zip(row.panels, widths):
            pw_i = int(round(pw))
            ph_i = int(round(rh))
            x    = int(round(cursor_x))
            y    = int(round(cursor_y))

            if show_labels:
                letter = (
                    panel.label
                    or (PANEL_LETTERS[panel_counter] if panel_counter < 26 else str(panel_counter + 1))
                )
                lbl = etree.SubElement(root, f"{{{SVG_NS}}}text", attrib={
                    "x": str(x + 3),
                    "y": str(y + font_label),
                    "font-size": str(font_label),
                    "font-weight": "bold",
                    "font-family": "Arial, Helvetica, sans-serif",
                    "fill": "#000000",
                })
                lbl.text = letter

            if show_titles:
                title = panel.title or auto_title(panel.file)
                lines = title.split("\n")
                line_h = font_title + 2   # px between baselines
                # Position first baseline so the block is vertically centred
                # within the label area (font_label px tall)
                y0 = y + font_label - (len(lines) - 1) * line_h // 2
                ttl = etree.SubElement(root, f"{{{SVG_NS}}}text", attrib={
                    "x": str(x + pw_i // 2),
                    "y": str(y0),
                    "font-size": str(font_title),
                    "font-family": "Arial, Helvetica, sans-serif",
                    "text-anchor": "middle",
                    "fill": "#444444",
                })
                if len(lines) == 1:
                    ttl.text = lines[0]
                else:
                    for i, line in enumerate(lines):
                        ts = etree.SubElement(ttl, f"{{{SVG_NS}}}tspan", attrib={
                            "x": str(x + pw_i // 2),
                            "dy": "0" if i == 0 else str(line_h),
                        })
                        ts.text = line

            # Embed raster image
            raw = load_png_bytes(panel.file)
            href = b64_uri(raw)
            etree.SubElement(root, f"{{{SVG_NS}}}image", attrib={
                "x": str(x),
                "y": str(y + effective_label_h),
                "width":  str(pw_i),
                "height": str(ph_i),
                "href": href,
                f"{{{XLINK_NS}}}href": href,
                "preserveAspectRatio": "xMidYMid meet",
            })

            cursor_x    += pw + gap
            panel_counter += 1

        cursor_y += rh + effective_label_h + gap

    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Assemble images into a manuscript-ready figure PDF (A4 portrait by default).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input
    in_grp = p.add_argument_group("Input (pick one)")
    in_grp.add_argument("files", nargs="*", metavar="FILE",
                        help="Explicit image files to include, in order")
    in_grp.add_argument("--glob", metavar="PATTERN",
                        help='Glob pattern, e.g. "plots/*_heatmap.png"')
    in_grp.add_argument("--input-dir", metavar="DIR",
                        help="Directory — all supported images inside are used")
    in_grp.add_argument("--layout", metavar="JSON",
                        help="JSON file defining rows / panel sizes / titles")

    # Grid (auto-mode only)
    p.add_argument("--cols", type=int, default=4,
                   help="Columns for auto-grid (ignored when --layout or MANUAL_LAYOUT is set)")

    # Page
    pg = p.add_argument_group("Page")
    pg.add_argument("--landscape", action="store_true",
                    help="A4 landscape (1123 × 794 px)")
    pg.add_argument("--width",  type=int, default=None,
                    help="Custom page width in SVG px (96 px = 1 inch)")
    pg.add_argument("--height", type=int, default=None,
                    help="Custom page height in SVG px")
    pg.add_argument("--margin", type=int, default=MARGIN)
    pg.add_argument("--gap",    type=int, default=GAP,
                    help="Gap between panels (SVG px)")

    # Labels / titles
    an = p.add_argument_group("Annotations")
    an.add_argument("--figure-title", metavar="TEXT",
                    help="Optional figure title printed above the grid")
    an.add_argument("--no-labels", action="store_true",
                    help="Suppress bold panel letters (A, B, …)")
    an.add_argument("--no-titles", action="store_true",
                    help="Suppress per-panel filename-derived titles")
    an.add_argument("--label-h",    type=int, default=LABEL_H,
                    help="Height reserved for labels/titles above each panel")
    an.add_argument("--font-label", type=int, default=FONT_LABEL)
    an.add_argument("--font-title", type=int, default=FONT_TITLE)

    # Output
    p.add_argument("--out", default="figure.pdf",
                   help="Output PDF path")
    p.add_argument("--out-dir", default=None,
                   help="Directory that will contain pdf/ and svg/ subdirs "
                        "(default: parent of --out)")

    return p.parse_args()


def main() -> None:
    args = parse_args()

    # Page dimensions
    page_w = 1123 if args.landscape else PAGE_W
    page_h = 794  if args.landscape else PAGE_H
    if args.width:  page_w = args.width
    if args.height: page_h = args.height

    # Resolve input files (only needed for auto-grid / validation)
    base_dir = Path(args.input_dir).resolve() if args.input_dir else Path.cwd()
    files = collect_files(
        args.files,
        args.glob,
        Path(args.input_dir).resolve() if args.input_dir else None,
    ) if (args.files or args.glob or args.input_dir) else []

    # Resolve layout
    layout, json_title = resolve_layout(
        args, files, base_dir,
        available_w=page_w - 2 * args.margin,
    )

    if not layout or not any(r.panels for r in layout):
        sys.exit("ERROR: no panels to assemble.")

    figure_title = args.figure_title or json_title

    # Print summary
    n_panels = sum(len(r.panels) for r in layout)
    orientation = "landscape" if page_w > page_h else "portrait"
    print(f"Page   : {page_w} × {page_h} px  ({orientation})")
    print(f"Panels : {n_panels} across {len(layout)} row(s)")
    available_w = page_w - 2 * args.margin
    for i, row in enumerate(layout):
        ws = panel_widths(row, available_w, args.gap)
        rh = row_height(row, ws)
        print(f"  Row {i+1}: {len(row.panels)} panel(s), height ≈ {rh:.0f} px  "
              f"[widths: {', '.join(f'{w:.0f}' for w in ws)}]")

    # Build and save
    fig = build_svg(
        layout,
        page_w=page_w,
        page_h=page_h,
        margin=args.margin,
        gap=args.gap,
        label_h=args.label_h,
        font_label=args.font_label,
        font_title=args.font_title,
        figure_title=figure_title,
        show_labels=not args.no_labels,
        show_titles=not args.no_titles,
    )

    out_pdf  = Path(args.out).resolve()
    base_out = Path(args.out_dir).resolve() if args.out_dir else out_pdf.parent
    pdf_dir  = base_out / "pdf"
    svg_dir  = base_out / "svg"
    pdf_dir.mkdir(parents=True, exist_ok=True)
    svg_dir.mkdir(parents=True, exist_ok=True)
    out_pdf_final = pdf_dir / out_pdf.name
    out_svg = svg_dir / out_pdf.with_suffix(".svg").name

    fig.save(str(out_svg))
    print(f"\nSVG → {out_svg}")

    # scale=1: SVG CSS px map to PDF pts at 0.75 ratio → exact physical A4
    cairosvg.svg2pdf(url=str(out_svg), write_to=str(out_pdf_final), scale=1)
    print(f"PDF → {out_pdf_final}")
    print("Done.")


if __name__ == "__main__":
    main()
