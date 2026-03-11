"""
_core.py
========
Processing logic extracted from coloc_analyzer.py.

This module contains **zero** Qt or Napari imports.  Every class and
function here is pure Python + NumPy + SciPy + matplotlib so that
it can be tested in isolation and driven from any front-end.
"""

import os
import glob
import datetime
import re
import warnings
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import tifffile
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.font_manager import FontProperties
from skimage.transform import resize
from skimage.measure import profile_line
from scipy.stats import pearsonr

# imagecodecs is required by tifffile to decompress LZW / Deflate (ZIP)
# compressed TIFFs, which is the default lossless export format in Leica LAS X.
try:
    import imagecodecs  # noqa: F401
except ImportError:
    raise ImportError(
        "'imagecodecs' is not installed.\n"
        "  tifffile requires this package to read LZW / Deflate-compressed\n"
        "  TIFFs (the lossless export format used by Leica LAS X).\n"
        "  Install it with:\n"
        "      pip install imagecodecs"
    )

# ================================================================
#  DEFAULT BASE DIRECTORY
# ================================================================
BASE_DIR = os.path.join(os.path.expanduser("~"), "Desktop", "Confocal_Workflow")


# ================================================================
#  SESSION CONFIGURATION
# ================================================================
class SessionConfig:
    """All run-time parameters.  Accepts keyword arguments so it can be
    driven from either the terminal wizard or the Napari dock widget."""

    def __init__(self, *,
                 input_dir=None,
                 exp_name=None,
                 crop_w=360,
                 crop_h=360,
                 do_bf=False,
                 do_zoom=False,
                 zoom_mag=3,
                 lbl_c="Protein-Cyan",
                 lbl_g="Protein-GFP",
                 lbl_m="Protein-mCherry",
                 lbl_merge="Merged",
                 lbl_zoom=None,
                 lbl_bf="BF",
                 do_quant=True,
                 do_intensity_profile=False,
                 panel_w_mm=20.0,
                 spacing_pt=3.0,
                 font_size=5.0,
                 sb_len_mm=2.646,
                 sb_um=10.0,
                 pixel_size_um=None):

        # --- Direct parameters ---
        self.input_dir = input_dir or os.path.join(BASE_DIR, "Input_Raw")
        today = datetime.datetime.now().strftime("%Y-%m-%d")
        self.exp_name = exp_name or f"Experiment_{today}"

        base_out = os.path.join(BASE_DIR, "Output_Coloc")
        self.output_dir = os.path.join(base_out, self.exp_name)
        os.makedirs(self.output_dir, exist_ok=True)

        self.crop_w = int(crop_w)
        self.crop_h = int(crop_h)

        self.do_bf = do_bf
        self.do_zoom = do_zoom
        self.zoom_mag = int(zoom_mag)
        self.zoom_size = self.crop_w // self.zoom_mag if self.do_zoom else 0

        self.lbl_c = lbl_c
        self.lbl_g = lbl_g
        self.lbl_m = lbl_m
        self.lbl_merge = lbl_merge
        self.lbl_zoom = lbl_zoom or (f"{self.zoom_mag}X enlarged" if self.do_zoom else "Enlarged")
        self.lbl_bf = lbl_bf

        self.do_quant = do_quant
        self.do_intensity_profile = do_intensity_profile

        self.panel_w_mm = float(panel_w_mm)
        self.spacing_pt = float(spacing_pt)
        self.font_size = float(font_size)
        self.sb_len_mm = float(sb_len_mm)
        self.sb_um = float(sb_um)
        self.pixel_size_um = float(pixel_size_um) if pixel_size_um else None

        # --- Derived values (identical to coloc_analyzer.py lines 158-183) ---
        self.aspect_ratio = self.crop_h / self.crop_w
        self.panel_h_mm = self.panel_w_mm * self.aspect_ratio

        self.panel_w_inch = self.panel_w_mm / 25.4
        self.panel_h_inch = self.panel_h_mm / 25.4
        self.spacing_inch = self.spacing_pt / 72.0

        self.dpi = self.crop_w / self.panel_w_inch

        # Scale bar: if pixel size is known, compute from µm directly;
        # otherwise fall back to the manual mm-based approach.
        if self.pixel_size_um and self.pixel_size_um > 0:
            self.sb_px_w = self.sb_um / self.pixel_size_um
            # Derive sb_len_mm so downstream print-size calculations stay valid
            self.sb_len_mm = (self.sb_px_w / self.crop_w) * self.panel_w_mm
        else:
            self.sb_px_w = (self.sb_len_mm / self.panel_w_mm) * self.crop_w

        self.sb_px_h = (0.5 / self.panel_w_mm) * self.crop_w
        self.sb_margin = (1.0 / self.panel_w_mm) * self.crop_w

        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['ps.fonttype'] = 42
        plt.rcParams['font.family'] = 'Arial'

        self.font_prop = FontProperties()
        self.font_prop.set_family('Arial')
        self.font_prop.set_size(self.font_size)
        self.font_prop.set_weight('bold')

        self.MIX_CYAN = np.array([0.0, 1.0, 1.0])
        self.MIX_GREEN = np.array([0.0, 1.0, 0.0])
        self.MIX_MAGENTA = np.array([1.0, 0.0, 1.0])

        self.zoom_upscaled_size = self.crop_w if self.do_zoom else 0
        self.zoom_panel_w_inch = self.panel_h_inch if self.do_zoom else 0


# ================================================================
#  XML METADATA DETECTOR
# ================================================================
class XMLMetadataDetector:
    """Parses Leica LAS X exported XML metadata to map channel indices to roles."""

    def __init__(self):
        self.lut_map = {
            'cyan':    'cyan',
            'blue':    'cyan',
            'green':   'green',
            'magenta': 'mag',
            'red':     'mag',
            'gray':    'bf',
            'grey':    'bf',
        }

    def assign_roles(self, files, log=print):
        """Map a list of TIF files to their respective roles using Leica XML.

        Parameters
        ----------
        files : list[str]
            Paths to ``_chXX.tif`` files belonging to one experiment set.
        log : callable
            Function used to emit status messages (default: ``print``).
        """
        result = {'cyan': None, 'green': None, 'mag': None, 'bf': None}

        for tif_path in files:
            base_dir = os.path.dirname(tif_path)
            filename = os.path.basename(tif_path)

            match = re.search(r'(.*)_ch(\d+)\.tif$', filename, re.IGNORECASE)
            if not match:
                continue

            base_name = match.group(1)
            ch_idx = int(match.group(2))

            xml_candidates = [
                os.path.join(base_dir, f"{base_name}_Properties.xml"),
                os.path.join(base_dir, f"{base_name}.xml"),
                os.path.join(base_dir, "MetaData", f"{base_name}_Properties.xml"),
                os.path.join(base_dir, "MetaData", f"{base_name}.xml"),
            ]
            xml_path = next((p for p in xml_candidates if os.path.exists(p)), None)

            if not xml_path:
                log(f"   >> WARNING: No XML metadata found for {filename}.")
                continue

            try:
                tree = ET.parse(xml_path)
                root = tree.getroot()

                channels_node = root.find('.//Channels')
                if channels_node is None:
                    log(f"   >> ERROR: <Channels> tag missing in {os.path.basename(xml_path)}")
                    continue

                channels = channels_node.findall('./ChannelDescription')
                if ch_idx >= len(channels):
                    log(f"   >> ERROR: Index {ch_idx} exceeds XML channel count in {os.path.basename(xml_path)}.")
                    continue

                lut_name = channels[ch_idx].attrib.get('LUTName', '').lower()
                role = self.lut_map.get(lut_name)

                if role in result:
                    result[role] = tif_path
                    log(f"   >> METADATA: Mapped {filename} (LUT: {lut_name.capitalize()}) -> {role.upper()}")
                else:
                    log(f"   >> UNKNOWN LUT: {filename} has unmapped LUT '{lut_name}'. Ignored.")

            except Exception as e:
                log(f"   >> ERROR parsing XML {xml_path}: {e}")

        return result

    def extract_pixel_size(self, files, log=print):
        """Extract pixel size in µm/pixel from Leica XML metadata.

        Searches for the X-dimension descriptor in the XML accompanying
        the supplied TIF files.  Supports both ``_Properties.xml``
        (human-readable, ``Unit="µm"``, ``Voxel`` attribute) and the
        compact ``.xml`` format (``Unit="m"``, compute from Length /
        NumberOfElements).

        Returns
        -------
        float or None
            Pixel size in µm, or ``None`` if metadata is unavailable.
        """
        for tif_path in files:
            base_dir = os.path.dirname(tif_path)
            filename = os.path.basename(tif_path)

            match = re.search(r'(.*)_ch(\d+)\.tif$', filename, re.IGNORECASE)
            if not match:
                continue

            base_name = match.group(1)
            xml_candidates = [
                os.path.join(base_dir, f"{base_name}_Properties.xml"),
                os.path.join(base_dir, f"{base_name}.xml"),
                os.path.join(base_dir, "MetaData", f"{base_name}_Properties.xml"),
                os.path.join(base_dir, "MetaData", f"{base_name}.xml"),
            ]
            xml_path = next((p for p in xml_candidates if os.path.exists(p)), None)
            if not xml_path:
                continue

            try:
                tree = ET.parse(xml_path)
                root = tree.getroot()

                for dim in root.iter('DimensionDescription'):
                    dim_id = dim.attrib.get('DimID', '')
                    # Match X dimension: DimID="X" or DimID="1"
                    if dim_id not in ('X', '1'):
                        continue

                    unit = dim.attrib.get('Unit', '')
                    n_elements = int(dim.attrib.get('NumberOfElements', 0))
                    if n_elements == 0:
                        continue

                    # _Properties.xml: has Voxel attribute directly in µm
                    voxel = dim.attrib.get('Voxel')
                    if voxel is not None:
                        px_um = float(voxel)
                        log(f"   >> METADATA: Pixel size = {px_um:.4f} µm/px "
                            f"(from Voxel attribute)")
                        return px_um

                    # Compact .xml: Length in metres, compute manually
                    length = float(dim.attrib.get('Length', 0))
                    if length <= 0:
                        continue

                    if unit == 'm':
                        px_um = (length / n_elements) * 1e6
                    elif unit in ('µm', 'um', 'μm'):
                        px_um = length / n_elements
                    else:
                        log(f"   >> WARNING: Unknown dimension unit '{unit}' "
                            f"in {os.path.basename(xml_path)}")
                        continue

                    log(f"   >> METADATA: Pixel size = {px_um:.4f} µm/px "
                        f"(from Length/N, unit={unit})")
                    return px_um

            except Exception as e:
                log(f"   >> ERROR extracting pixel size from "
                    f"{os.path.basename(xml_path)}: {e}")

        log("   >> WARNING: Could not detect pixel size from XML metadata. "
            "Using manual scale bar (mm) fallback.")
        return None


# ================================================================
#  COLOCALIZATION QUANTIFIER
# ================================================================
def _threshold_max_entropy(image):
    """Kapur Maximum-Entropy thresholding (Kapur et al., 1985)."""
    hist, bin_edges = np.histogram(image.ravel(), bins=256, density=False)
    hist = hist.astype(np.float64)
    hist = hist / hist.sum()

    cum_sum = np.cumsum(hist)
    eps = np.finfo(np.float64).eps

    best_t = 0
    best_ent = -np.inf

    for t in range(1, 256):
        p_bg = cum_sum[t - 1]
        if p_bg < eps or p_bg > 1.0 - eps:
            continue
        p_fg = 1.0 - p_bg

        h_bg = hist[:t]
        h_bg = h_bg[h_bg > eps]
        ent_bg = -np.sum((h_bg / p_bg) * np.log(h_bg / p_bg))

        h_fg = hist[t:]
        h_fg = h_fg[h_fg > eps]
        ent_fg = -np.sum((h_fg / p_fg) * np.log(h_fg / p_fg))

        total_ent = ent_bg + ent_fg
        if total_ent > best_ent:
            best_ent = total_ent
            best_t = t

    return (bin_edges[best_t] + bin_edges[best_t + 1]) / 2.0


class Quantifier:
    """Compute Pearson's R and Manders' coefficients for a pair of
    normalised single-channel arrays."""

    @staticmethod
    def analyze(name, green_crop, mag_crop, log=print):
        flat_g = green_crop.flatten()
        flat_m = mag_crop.flatten()
        if len(flat_g) < 2:
            return None, None, None

        if np.std(flat_g) == 0 or np.std(flat_m) == 0:
            log("   >> WARNING: One or both channels are constant. Pearson R set to NaN.")
            pearson_r = float('nan')
        else:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pearson_r, _ = pearsonr(flat_g, flat_m)

        try:
            thresh_g = _threshold_max_entropy(green_crop)
            thresh_m = _threshold_max_entropy(mag_crop)
        except ValueError:
            thresh_g = 0.01
            thresh_m = 0.01

        mask_g = green_crop > thresh_g
        mask_m = mag_crop > thresh_m

        total_g = np.sum(green_crop)
        total_m = np.sum(mag_crop)

        m1 = np.sum(green_crop[mask_m]) / total_g if total_g != 0 else 0
        m2 = np.sum(mag_crop[mask_g]) / total_m if total_m != 0 else 0

        stats = {
            'Image_ID':             name,
            'Pearson_R':            round(pearson_r, 3),
            'Manders_M1 (G_in_M)': round(m1, 3),
            'Manders_M2 (M_in_G)': round(m2, 3),
            'Thresh_Green':         round(thresh_g, 4),
            'Thresh_Mag':           round(thresh_m, 4),
        }
        return stats, mask_g, mask_m


# ================================================================
#  UTILITY FUNCTIONS
# ================================================================
def scan_directory(input_dir):
    """Return sorted list of unique experiment base paths found under *input_dir*."""
    search_pattern = os.path.join(input_dir, "**", "*.tif")
    all_tifs = glob.glob(search_pattern, recursive=True)
    unique_bases = set()

    for f in all_tifs:
        if os.path.basename(f).startswith("._"):
            continue
        if "_ch" in f:
            base = f.split("_ch")[0]
            unique_bases.add(base)

    return sorted(unique_bases)


def load_flat(path):
    """Load a TIFF and flatten to a 2-D single-channel array."""
    if not path or not os.path.exists(path):
        return None
    raw = tifffile.imread(path)
    raw = np.squeeze(raw)
    while raw.ndim > 3:
        raw = np.max(raw, axis=0)
    if raw.ndim == 2:
        return raw
    if raw.ndim == 3:
        if raw.shape[-1] == 3:
            return np.max(raw, axis=-1)
        return np.max(raw, axis=0)
    return raw


def auto_contrast_limits(img):
    """Return (vmin, vmax) percentile limits for *img*."""
    if img is None or img.size == 0:
        return (0, 1)
    vmin, vmax = np.percentile(img, 0.1), np.percentile(img, 99.98)
    if vmax <= vmin:
        vmax = vmin + 1
    return (vmin, vmax)


def auto_contrast(img, limits=None):
    """Normalise *img* to [0, 1] using 0.1–99.98 percentile stretch.

    If *limits* is provided as a ``(vmin, vmax)`` tuple, those values
    are used instead of computing percentiles from *img*.
    """
    if img is None or img.size == 0:
        return None
    if limits is not None:
        vmin, vmax = limits
    else:
        vmin, vmax = auto_contrast_limits(img)
    return np.clip((img.astype(float) - vmin) / (vmax - vmin), 0, 1)


def make_rgb(norm, color_mix):
    """Convert a single-channel normalised image to a pseudocolor RGB."""
    if norm is None:
        return None
    return norm[..., np.newaxis] * color_mix


def get_box_center(shape_data):
    """Return (cy, cx) for a rectangle defined by its corner vertices."""
    y_vals = [p[0] for p in shape_data]
    x_vals = [p[1] for p in shape_data]
    return (int((min(y_vals) + max(y_vals)) / 2),
            int((min(x_vals) + max(x_vals)) / 2))


def safe_crop(img, y1, y2, x1, x2):
    """Crop *img* with bounds clamping.  Returns ``None`` if *img* is ``None``."""
    if img is None:
        return None
    h, w = img.shape[:2]
    y1, x1 = max(0, y1), max(0, x1)
    y2, x2 = min(h, y2), min(w, x2)
    if y1 >= y2 or x1 >= x2:
        return None
    return img[y1:y2, x1:x2]


def upscale_channel(crop, target_h, target_w):
    """Bicubic-upscale a 2-D array to *(target_h, target_w)*."""
    if crop is None:
        return None
    return resize(
        crop.astype(np.float64), (target_h, target_w),
        order=3, anti_aliasing=False, preserve_range=True)


# ================================================================
#  FIGURE BUILDER
# ================================================================
class FigureBuilder:
    """Groups all matplotlib output methods.  Depends only on a
    :class:`SessionConfig` instance — no Qt or Napari references."""

    def __init__(self, cfg):
        self.cfg = cfg

    # ----- layout helpers -----
    def _col_widths_inch(self, col_labels):
        widths = []
        for lbl in col_labels:
            if lbl == self.cfg.lbl_zoom:
                widths.append(self.cfg.zoom_panel_w_inch)
            else:
                widths.append(self.cfg.panel_w_inch)
        return widths

    def _make_axes_grid(self, fig, n_rows, n_cols, col_w, total_w, total_h):
        row_h = self.cfg.panel_h_inch
        axes = []
        for r in range(n_rows):
            row_top = r * (row_h + self.cfg.spacing_inch)
            bottom = 1.0 - (row_top + row_h) / total_h
            h_frac = row_h / total_h

            row_axes = []
            x = 0.0
            for ci in range(n_cols):
                left = x / total_w
                w_frac = col_w[ci] / total_w
                ax = fig.add_axes([left, bottom, w_frac, h_frac])
                row_axes.append(ax)
                x += col_w[ci] + self.cfg.spacing_inch
            axes.append(row_axes)
        return axes

    # ----- intensity plot -----
    def _render_intensity_axes(self, ax, intensity_data, mode='DUAL',
                               linewidth=0.8, tick_scale=1.0):
        ax.set_facecolor('none')
        ax.patch.set_alpha(0)

        dist = intensity_data['Distance_px'].values
        if mode == "TRIPLE" and self.cfg.lbl_c in intensity_data:
            ax.plot(dist, intensity_data[self.cfg.lbl_c].values,
                    color='cyan', linewidth=linewidth)
        if self.cfg.lbl_g in intensity_data:
            ax.plot(dist, intensity_data[self.cfg.lbl_g].values,
                    color='lime', linewidth=linewidth)
        if self.cfg.lbl_m in intensity_data:
            ax.plot(dist, intensity_data[self.cfg.lbl_m].values,
                    color='magenta', linewidth=linewidth)

        ax.set_xlim(dist[0], dist[-1])
        ax.set_ylim(0, 1.05)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_edgecolor('black')
        ax.spines['left'].set_linewidth(0.5 * tick_scale)
        ax.spines['bottom'].set_edgecolor('black')
        ax.spines['bottom'].set_linewidth(0.5 * tick_scale)

        ts = self.cfg.font_size * 0.7 * tick_scale
        ax.tick_params(axis='y', colors='black', labelsize=ts,
                       length=2 * tick_scale, pad=1)
        ax.tick_params(axis='x', colors='black', labelsize=ts,
                       length=2 * tick_scale, pad=1, labelbottom=False)
        ax.yaxis.set_label_position('left')
        ax.yaxis.tick_left()

        pos = ax.get_position()
        left_margin = pos.width * 0.10
        ax.set_position([
            pos.x0 + left_margin,
            pos.y0,
            pos.width - left_margin,
            pos.height,
        ])

    def add_scale_bar_patch(self, ax):
        x_start = self.cfg.crop_w - self.cfg.sb_px_w - self.cfg.sb_margin
        y_start = self.cfg.crop_h - self.cfg.sb_margin - self.cfg.sb_px_h
        ax.add_patch(patches.Rectangle(
            (x_start, y_start), self.cfg.sb_px_w, self.cfg.sb_px_h,
            linewidth=0, edgecolor='none', facecolor='white'))

    # ----- individual panel -----
    def save_individual_and_local_panel(self, name, c, g, m, merge, bf,
                                        zoom, zoom_coords, mask_g, mask_m,
                                        mode="DUAL",
                                        intensity_data=None,
                                        line_coords_crop=None):
        save_path = os.path.join(self.cfg.output_dir, name)
        os.makedirs(save_path, exist_ok=True)

        def to_uint8(arr):
            return (arr * 255).astype(np.uint8)

        if mode == "TRIPLE":
            tifffile.imwrite(os.path.join(save_path, "1_Cyan.tif"), to_uint8(c))
        tifffile.imwrite(os.path.join(save_path, "2_Green.tif"), to_uint8(g))
        tifffile.imwrite(os.path.join(save_path, "3_Magenta.tif"), to_uint8(m))
        tifffile.imwrite(os.path.join(save_path, "4_Merge.tif"), to_uint8(merge))
        if zoom is not None:
            tifffile.imwrite(os.path.join(save_path, "5_Zoom.tif"), to_uint8(zoom))
        if self.cfg.do_bf:
            tifffile.imwrite(os.path.join(save_path, "6_BF.tif"), to_uint8(bf))

        if mask_g is not None and mask_m is not None:
            qc_path = os.path.join(save_path, "QC_Masks")
            os.makedirs(qc_path, exist_ok=True)
            tifffile.imwrite(
                os.path.join(qc_path, "Mask_Green_Otsu.tif"), to_uint8(mask_g))
            tifffile.imwrite(
                os.path.join(qc_path, "Mask_Mag_Otsu.tif"), to_uint8(mask_m))

        plot_list = []
        if mode == "TRIPLE":
            plot_list.append((c, self.cfg.lbl_c))
        plot_list.append((g, self.cfg.lbl_g))
        plot_list.append((m, self.cfg.lbl_m))
        plot_list.append((merge, self.cfg.lbl_merge))
        if zoom is not None:
            plot_list.append((zoom, self.cfg.lbl_zoom))
        if self.cfg.do_bf:
            plot_list.append((bf, self.cfg.lbl_bf))

        has_profile = (intensity_data is not None and not intensity_data.empty)
        if has_profile:
            plot_list.append((None, '__intensity__'))

        n_cols = len(plot_list)
        col_labels = [lbl for _, lbl in plot_list]
        col_w = self._col_widths_inch(col_labels)

        total_w = sum(col_w) + (n_cols - 1) * self.cfg.spacing_inch
        total_h = self.cfg.panel_h_inch

        fig = plt.figure(figsize=(total_w, total_h))
        axes_grid = self._make_axes_grid(fig, 1, n_cols, col_w, total_w, total_h)
        axes = axes_grid[0]

        for i, ax in enumerate(axes):
            img_data, label = plot_list[i]

            if label == '__intensity__':
                self._render_intensity_axes(ax, intensity_data, mode=mode,
                                            linewidth=0.8, tick_scale=1.0)
                continue

            ax.imshow(img_data, aspect='auto', interpolation='none')
            ax.axis('off')
            ax.text(0.05, 0.95, label, transform=ax.transAxes, color='white',
                    fontproperties=self.cfg.font_prop, ha='left', va='top')
            if i == 0:
                self.add_scale_bar_patch(ax)

            if label == self.cfg.lbl_merge:
                if zoom_coords is not None:
                    ry, rx = zoom_coords
                    ax.add_patch(patches.Rectangle(
                        (rx, ry), self.cfg.zoom_size, self.cfg.zoom_size,
                        linewidth=1, edgecolor='white',
                        facecolor='none', linestyle='--'))
                if line_coords_crop is not None:
                    (ly0, lx0), (ly1, lx1) = line_coords_crop
                    ax.plot([lx0, lx1], [ly0, ly1],
                            color='white', linewidth=0.8, linestyle='--')

        plt.savefig(os.path.join(save_path, "Panel_View.pdf"), dpi=self.cfg.dpi)
        plt.savefig(os.path.join(save_path, "Panel_View.png"), dpi=self.cfg.dpi)
        plt.close(fig)

    # ----- finalisation -----
    def finalize_analysis(self, gallery, stats_log, log=print):
        if gallery:
            self.build_global_montage(gallery, log=log)

        if self.cfg.do_quant and stats_log:
            df = pd.DataFrame(stats_log)
            csv_path = os.path.join(self.cfg.output_dir,
                                    f"{self.cfg.exp_name}_QUANTIFICATION.csv")
            df.to_csv(csv_path, index=False)
            log(f">> STATS SAVED: {csv_path}")

    # ----- relabeling helpers -----
    def regenerate_all_panels(self, gallery, log=print):
        for item in gallery:
            name = item['name']
            mode = item['mode']
            self.save_individual_and_local_panel(
                name,
                item['cyan'],
                item['green'],
                item['mag'],
                item['merge'],
                item.get('bf'),
                item.get('zoom'),
                item.get('zoom_coords'),
                None, None,
                mode=mode,
                intensity_data=item.get('intensity_data'),
                line_coords_crop=item.get('line_coords'),
            )
        log(f"   >> {len(gallery)} panel view(s) regenerated.")

    def regenerate_intensity_csvs(self, gallery, log=print):
        for item in gallery:
            raw = item.get('raw_profiles')
            if raw is None:
                continue
            new_dict = {'Distance_px': raw['distance']}
            if 'cyan' in raw:
                new_dict[self.cfg.lbl_c] = raw['cyan']
            new_dict[self.cfg.lbl_g] = raw['green']
            new_dict[self.cfg.lbl_m] = raw['mag']
            new_df = pd.DataFrame(new_dict)
            item['intensity_data'] = new_df
            save_dir = os.path.join(self.cfg.output_dir, item['name'])
            new_df.to_csv(
                os.path.join(save_dir, "Intensity_Profile.csv"), index=False)
        if any(item.get('raw_profiles') for item in gallery):
            log("   >> Intensity CSVs updated with new labels.")

    # ----- summary montage -----
    def build_global_montage(self, gallery, log=print):
        log(">>> GENERATING SUMMARY MONTAGES (PDF + PNG)...")

        n_rows = len(gallery)
        has_cyan = any(x['mode'] == "TRIPLE" for x in gallery)
        has_zoom = any(x['zoom'] is not None for x in gallery)
        has_intensity = any(x.get('intensity_data') is not None for x in gallery)

        col_keys = []
        col_labels = []
        if has_cyan:
            col_keys.append('cyan');  col_labels.append(self.cfg.lbl_c)
        col_keys.append('green');     col_labels.append(self.cfg.lbl_g)
        col_keys.append('mag');       col_labels.append(self.cfg.lbl_m)
        col_keys.append('merge');     col_labels.append(self.cfg.lbl_merge)
        if has_zoom:
            col_keys.append('zoom');  col_labels.append(self.cfg.lbl_zoom)
        if self.cfg.do_bf:
            col_keys.append('bf');    col_labels.append(self.cfg.lbl_bf)
        if has_intensity:
            col_keys.append('__intensity__'); col_labels.append('__intensity__')

        n_cols = len(col_keys)
        col_w = self._col_widths_inch(col_labels)
        total_w = sum(col_w) + (n_cols - 1) * self.cfg.spacing_inch
        total_h = (n_rows * self.cfg.panel_h_inch +
                   (n_rows - 1) * self.cfg.spacing_inch)

        fig = plt.figure(figsize=(total_w, total_h))
        axes = self._make_axes_grid(fig, n_rows, n_cols, col_w, total_w, total_h)

        h, w = self.cfg.crop_h, self.cfg.crop_w
        black_main = np.zeros((h, w, 3))
        zup = self.cfg.zoom_upscaled_size if self.cfg.zoom_upscaled_size > 0 else h
        black_zoom = np.zeros((zup, zup, 3))

        for r, item in enumerate(gallery):
            for ci, key in enumerate(col_keys):
                ax = axes[r][ci]

                if key == '__intensity__':
                    idata = item.get('intensity_data')
                    if idata is not None and not idata.empty:
                        self._render_intensity_axes(
                            ax, idata, mode=item.get('mode', 'DUAL'),
                            linewidth=0.5, tick_scale=0.7)
                    else:
                        ax.set_facecolor('none')
                        ax.patch.set_alpha(0)
                        ax.axis('off')
                    continue

                img_data = item.get(key)
                if img_data is None:
                    img_data = black_zoom if key == 'zoom' else black_main

                ax.imshow(img_data, aspect='auto', interpolation='none')
                ax.axis('off')

                if r == 0:
                    ax.text(0.05, 0.95, col_labels[ci], transform=ax.transAxes,
                            color='white', fontproperties=self.cfg.font_prop,
                            ha='left', va='top')

                if r == 0 and ci == 0:
                    self.add_scale_bar_patch(ax)

                if key == 'merge':
                    if item.get('zoom_coords') is not None:
                        ry, rx = item['zoom_coords']
                        ax.add_patch(patches.Rectangle(
                            (rx, ry), self.cfg.zoom_size, self.cfg.zoom_size,
                            linewidth=0.5, edgecolor='white',
                            facecolor='none', linestyle='--'))
                    if item.get('line_coords') is not None:
                        (ly0, lx0), (ly1, lx1) = item['line_coords']
                        ax.plot([lx0, lx1], [ly0, ly1],
                                color='white', linewidth=0.5, linestyle='--')

        base_name = os.path.join(self.cfg.output_dir,
                                 f"{self.cfg.exp_name}_SUMMARY_MONTAGE")
        plt.savefig(f"{base_name}.pdf", dpi=self.cfg.dpi)
        plt.savefig(f"{base_name}.png", dpi=self.cfg.dpi)
        log(f">> MONTAGE SAVED: {base_name}.pdf")
        plt.close(fig)
