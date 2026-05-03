"""Build a vertical montage by picking one row per experiment.

By default each row is rasterized with :class:`FigureBuilder` (titles, scale bar,
overlays, intensity column) — same as ``Panel_View``. Optionally paste raw
``1_Cyan.tif`` … strips without typography.

Row order follows each experiment's session pickle gallery (same as SUMMARY_MONTAGE).
"""

import glob
import json
import os
from typing import Dict, Iterator, List, Optional, Sequence, Tuple, Union

from PIL import Image

from ._core import FigureBuilder, load_session_record

SelectionPair = Tuple[str, int]

TIFF_BY_KEY: Dict[str, str] = {
    "cyan": "1_Cyan.tif",
    "green": "2_Green.tif",
    "mag": "3_Magenta.tif",
    "merge": "4_Merge.tif",
    "zoom": "5_Zoom.tif",
    "bf": "6_BF.tif",
}


def _session_record_pickle_path(exp_dir: str) -> str:
    pkls = glob.glob(os.path.join(exp_dir, "*_session_record.pkl"))
    if not pkls:
        raise FileNotFoundError(
            f"No *_session_record.pkl in experiment folder: {exp_dir}"
        )
    if len(pkls) > 1:
        pkls.sort()
    return pkls[0]


def list_experiment_dirs(date_dir: str) -> List[str]:
    """Return sorted absolute paths to immediate subdirs that contain a session pickle."""
    if not os.path.isdir(date_dir):
        raise FileNotFoundError(f"Date folder not found: {date_dir}")
    names = sorted(
        entry for entry in os.listdir(date_dir)
        if os.path.isdir(os.path.join(date_dir, entry))
        and glob.glob(os.path.join(date_dir, entry, "*_session_record.pkl"))
    )
    return [os.path.join(date_dir, name) for name in names]


def ordered_image_basenames(exp_dir: str) -> List[str]:
    """Basenames in gallery / montage row order (1-based row index applies here)."""
    pkl = _session_record_pickle_path(exp_dir)
    _, gallery, _ = load_session_record(pkl)
    return [item["name"] for item in gallery]


def panel_png_path(exp_dir: str, basename: str) -> str:
    path = os.path.join(exp_dir, basename, "Panel_View.png")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing Panel_View.png: {path}")
    return path


def _iter_tiff_keys(item: dict) -> Iterator[str]:
    """Mirror FigureBuilder._iter_panel_items but only image keys (no intensity)."""
    channels = item.get("included_channels")
    if channels is None:
        channels = {
            "cyan": item.get("cyan") is not None,
            "green": item.get("green") is not None,
            "mag": item.get("mag") is not None,
        }
    include_c = bool(channels.get("cyan")) and item.get("cyan") is not None
    include_g = bool(channels.get("green")) and item.get("green") is not None
    include_m = bool(channels.get("mag")) and item.get("mag") is not None

    if include_c:
        yield "cyan"
    if include_g:
        yield "green"
    if include_m:
        yield "mag"

    if item.get("merge") is not None:
        yield "merge"
    if item.get("zoom") is not None:
        yield "zoom"
    if item.get("bf") is not None:
        yield "bf"


def _build_row_from_channel_tifs(
    image_dir: str,
    item: dict,
    channel_spacing_px: int,
) -> Tuple[Optional[Image.Image], List[str]]:
    """One horizontal strip from channel TIFFs, common row height, adjustable gaps."""
    spacing = max(int(channel_spacing_px), 0)
    keys = list(_iter_tiff_keys(item))
    if not keys:
        return None, []

    raw: List[Image.Image] = []
    paths: List[str] = []
    for key in keys:
        fname = TIFF_BY_KEY[key]
        path = os.path.join(image_dir, fname)
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"Expected channel TIFF for cross-montage: {path}"
            )
        raw.append(Image.open(path).convert("RGBA"))
        paths.append(path)

    try:
        target_h = max(im.height for im in raw)
        pieces: List[Image.Image] = []
        for im in raw:
            if im.height != target_h:
                new_w = round(im.width * target_h / im.height)
                resized = im.resize((new_w, target_h), Image.LANCZOS)
                im.close()
                pieces.append(resized)
            else:
                pieces.append(im)

        inner_w = sum(p.width for p in pieces) + spacing * (len(pieces) - 1)
        row = Image.new("RGBA", (inner_w, target_h), (0, 0, 0, 0))
        x = 0
        for idx, p in enumerate(pieces):
            row.paste(p, (x, 0))
            x += p.width + (spacing if idx < len(pieces) - 1 else 0)

        for p in pieces:
            p.close()

        return row, paths
    except Exception:
        for im in raw:
            im.close()
        raise


def _experiment_fb_gallery(exp_abs: str, cache: dict) -> Tuple[FigureBuilder, List[dict]]:
    if exp_abs not in cache:
        cfg, gallery, _ = load_session_record(_session_record_pickle_path(exp_abs))
        cache[exp_abs] = (FigureBuilder(cfg), gallery)
    return cache[exp_abs]


def build_custom_montage(
    selections: Sequence[SelectionPair],
    out_base: str,
    gap_px: int = 20,
    gap_horizontal_px: int = 0,
    channel_spacing_px: int = 12,
    labeled_strip: bool = True,
    fallback_panel_view: bool = True,
    zoom_roi_stroke_pt: Optional[float] = None,
    dpi: Union[int, Tuple[int, int]] = 300,
    save_spec: bool = True,
    log=print,
) -> str:
    """Stack one strip per selection (top to bottom).

    With *labeled_strip* True (default), each row is re-rendered with
    ``FigureBuilder`` — same titles, scale bar, overlays and intensity subplot
    as ``Panel_View``. *channel_spacing_px* maps to subplot column spacing versus
    the export *dpi*.

    With *labeled_strip* False, rows paste ``1_Cyan.tif`` … (no typography).
    Missing pieces fall back to ``Panel_View.png`` when permitted.

    *zoom_roi_stroke_pt* — when positive, override the dashed zoom ROI stroke
    width in matplotlib points; otherwise use the exact ``Panel_View`` default.
    The ROI square area remains each session's ``cfg.zoom_size``.
    """
    if not selections:
        raise ValueError("No selections: provide at least one (experiment_dir, row).")

    stroke_opt = None
    if zoom_roi_stroke_pt is not None:
        z = float(zoom_roi_stroke_pt)
        if z > 0:
            stroke_opt = z

    gap_v = max(int(gap_px), 0)
    margin_h = max(int(gap_horizontal_px), 0)
    ch_space = max(int(channel_spacing_px), 0)
    if isinstance(dpi, tuple):
        dpi_t = dpi
    else:
        dpi_t = (int(dpi), int(dpi))

    sess_cache: dict = {}
    rows_meta: List[dict] = []
    row_images: List[Image.Image] = []

    for exp_dir, row_1 in selections:
        exp_abs = os.path.abspath(exp_dir)
        if row_1 < 1:
            raise ValueError(
                f"Row index must be >= 1, got {row_1} for {exp_abs}"
            )

        fb, gallery = _experiment_fb_gallery(exp_abs, sess_cache)
        if row_1 > len(gallery):
            raise ValueError(
                f"Experiment {os.path.basename(exp_abs)!r}: row {row_1} "
                f"out of range (1-{len(gallery)} processed images)."
            )
        item = gallery[row_1 - 1]
        basename = item["name"]
        image_dir = os.path.join(exp_abs, basename)

        png_path = os.path.join(exp_abs, basename, "Panel_View.png")
        tif_paths: List[str] = []
        row_strip: Optional[Image.Image] = None
        src_mode = "panel_view"

        export_dpi = float(fb.cfg.dpi)
        axis_sp_in = (
            ch_space / export_dpi
            if export_dpi > 0 else fb.cfg.spacing_inch)

        if labeled_strip:
            try:
                row_strip = fb.render_panel_row_pil(
                    item,
                    axis_spacing_inch=axis_sp_in,
                    output_dpi=None,
                    zoom_roi_stroke_pt=stroke_opt,
                    log=log,
                )
                if row_strip is not None:
                    src_mode = "labeled_matplotlib"
            except Exception as err:
                if not fallback_panel_view:
                    raise
                log(f"   >> Cross-exp: matplotlib row render failed ({err}); fallback.")

        if row_strip is None and not labeled_strip:
            try:
                row_strip, tif_paths = _build_row_from_channel_tifs(
                    image_dir, item, ch_space)
                if row_strip is not None:
                    src_mode = "channel_tifs"
            except FileNotFoundError as err:
                if not fallback_panel_view:
                    raise
                log(f"   >> Cross-exp: raw TIFF row failed ({err}); trying Panel_View.")

        if row_strip is None:
            if not os.path.isfile(png_path):
                raise FileNotFoundError(
                    f"Cross-exp: could not assemble row (render/TIFF fallback); "
                    f"missing {png_path}"
                )
            row_strip = Image.open(png_path).convert("RGBA")
            src_mode = "panel_view"

        rows_meta.append({
            "experiment_folder": os.path.basename(exp_abs),
            "experiment_path": exp_abs,
            "row_1based": row_1,
            "image_basename": basename,
            "assemble_mode": src_mode,
            "panel_png": png_path,
            "channel_tifs": tif_paths,
        })
        row_images.append(row_strip)

    max_w = max(r.width for r in row_images)
    centered_rows: List[Image.Image] = []
    for rimg in row_images:
        if rimg.width < max_w:
            bg = Image.new("RGBA", (max_w, rimg.height), (0, 0, 0, 0))
            bg.paste(rimg, ((max_w - rimg.width) // 2, 0))
            rimg.close()
            centered_rows.append(bg)
        else:
            centered_rows.append(rimg)

    content_w = max_w + 2 * margin_h
    total_h_comb = (
        sum(r.height for r in centered_rows)
        + gap_v * (len(centered_rows) - 1))

    canvas = Image.new("RGBA", (content_w, total_h_comb), (0, 0, 0, 0))
    y = 0
    for idx, rim in enumerate(centered_rows):
        canvas.paste(rim, (margin_h, y))
        y += rim.height + (gap_v if idx < len(centered_rows) - 1 else 0)

    png_path_out = f"{out_base}.png"
    pdf_path_out = f"{out_base}.pdf"
    canvas.save(png_path_out, dpi=dpi_t)
    canvas.save(pdf_path_out, dpi=dpi_t)
    log(f">> Custom cross-experiment montage: {png_path_out} / {pdf_path_out}")

    if save_spec:
        spec_path = f"{out_base}_spec.json"
        payload = {
            "out_base": out_base,
            "gap_vertical_px": gap_v,
            "gap_horizontal_px": margin_h,
            "channel_spacing_px": ch_space,
            "labeled_strip": labeled_strip,
            "fallback_panel_view": fallback_panel_view,
            "zoom_roi_stroke_pt": stroke_opt,
            "dpi": list(dpi_t),
            "panels": rows_meta,
        }
        with open(spec_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        log(f">> Spec saved: {spec_path}")

    for rim in centered_rows:
        rim.close()
    canvas.close()

    return png_path_out
