"""
coloc_analyzer.py
=================
Interactive fluorescence colocalization analysis tool built on Napari.

Workflow
--------
1. Run this script from a terminal::

       python coloc_analyzer.py

2. Answer the setup wizard prompts (input folder, crop size, labels, etc.).
3. In Napari, drag the **yellow** box over your region of interest.
   If zoom is enabled, drag the **cyan** box over the sub-region to enlarge.
4. Press **[c]** to crop, save, and advance to the next image set.
   Press **[s]** to skip the current image set.
5. After all images are processed, a summary montage (PDF + PNG) and a
   quantification CSV are saved automatically.

Input
-----
Multi-channel confocal images stored as separate single-channel TIFF files
named with a ``_ch<N>.tif`` suffix, e.g.::

    MyExperiment_s001_ch00.tif   (Cyan / first fluorophore)
    MyExperiment_s001_ch01.tif   (Green / second fluorophore)
    MyExperiment_s001_ch02.tif   (Magenta / third fluorophore)
    MyExperiment_s001_ch03.tif   (Brightfield, optional)

Files belonging to the same field of view must share a common prefix before
``_ch``.  Subdirectories are scanned recursively.

Output (per image)
------------------
``<output_dir>/<experiment_name>/<image_basename>/``

* ``1_Cyan.tif``, ``2_Green.tif``, ``3_Magenta.tif``, ``4_Merge.tif``
* ``5_Zoom.tif``  (if zoom is enabled)
* ``6_BF.tif``    (if brightfield is enabled)
* ``Panel_View.pdf`` / ``Panel_View.png``  -- publication-ready single-row panel
* ``QC_Masks/``   -- Otsu threshold masks (if quantification is enabled)

Output (global)
---------------
``<output_dir>/<experiment_name>/``

* ``<experiment_name>_SUMMARY_MONTAGE.pdf`` / ``.png``
* ``<experiment_name>_QUANTIFICATION.csv``  -- Pearson R, Manders M1/M2

Dependencies
------------
See ``requirements.txt``.

License
-------
MIT License -- see ``LICENSE``.
"""

import napari
import os
import glob
import numpy as np
import pandas as pd
import tifffile
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.font_manager import FontProperties
from skimage.filters import threshold_otsu
from skimage.transform import resize
from scipy.stats import pearsonr
import datetime

# ================================================================
#  DEFAULT BASE DIRECTORY
#  Change this to match your system, or enter a custom path in
#  the wizard prompt at runtime.
# ================================================================
BASE_DIR = os.path.join(os.path.expanduser("~"), "Desktop", "Confocal_Workflow")


# ----------------------------------------------------------------
#  Utility
# ----------------------------------------------------------------
def get_user_input(prompt, default_val):
    """Prompt the user for input, returning *default_val* if nothing entered."""
    val = input(f"   {prompt} [{default_val}]: ")
    return val.strip() if val.strip() != "" else default_val


# ================================================================
#  SESSION CONFIGURATION  (interactive wizard)
# ================================================================
class SessionConfig:
    """Collect all run-time parameters through an interactive command-line wizard."""

    def __init__(self):
        print("--- EXPERIMENT SETUP WIZARD ---")

        # 1. INPUT DIRECTORY ------------------------------------------------
        default_in = os.path.join(BASE_DIR, "Input_Raw")
        print("\n--- DATA SOURCE ---")
        self.input_dir = get_user_input("Input Folder Path", default_in)
        self.input_dir = self.input_dir.replace("'", "").replace('"', "").strip()

        if not os.path.exists(self.input_dir):
            print(f"\nCRITICAL ERROR: Folder not found: {self.input_dir}")
            exit()

        # 2. OUTPUT DIRECTORY -----------------------------------------------
        print("\n--- DATA DESTINATION ---")
        today = datetime.datetime.now().strftime("%Y-%m-%d")
        default_exp = f"Experiment_{today}"
        self.exp_name = get_user_input("Experiment Name", default_exp)

        base_out = os.path.join(BASE_DIR, "Output_Coloc")
        self.output_dir = os.path.join(base_out, self.exp_name)
        os.makedirs(self.output_dir, exist_ok=True)

        # 3. CROP GEOMETRY --------------------------------------------------
        print("\n--- CROP GEOMETRY ---")
        self.crop_w = int(get_user_input("Crop Width (px)", 360))
        self.crop_h = int(get_user_input("Crop Height (px)", 360))

        # 4. OPTIONAL FEATURES ----------------------------------------------
        print("\n--- EXTRA FEATURES ---")
        self.do_bf = get_user_input(
            "Include Brightfield Panel? [y/n]", "n").lower() == 'y'
        self.do_zoom = get_user_input(
            "Include Zoom/Enlargement Panel? [y/n]", "n").lower() == 'y'
        self.zoom_mag  = 1
        self.zoom_size = 0   # square side length in pixels

        if self.do_zoom:
            mag = int(get_user_input("   > Magnification Factor (e.g. 3)", "3"))
            self.zoom_mag  = mag
            self.zoom_size = self.crop_w // mag   # square ROI side

        # 5. CHANNEL LABELS -------------------------------------------------
        print("\n--- CHANNEL LABELS ---")
        self.lbl_c = get_user_input("Cyan Label",    "Protein-Cyan")
        self.lbl_g = get_user_input("Green Label",   "Protein-GFP")
        self.lbl_m = get_user_input("Magenta Label", "Protein-mCherry")

        default_zoom_lbl = f"{self.zoom_mag}X enlarged" if self.do_zoom else "Enlarged"
        self.lbl_merge = get_user_input("Merge Label", "Merged")
        if self.do_zoom:
            self.lbl_zoom = get_user_input("Enlargement Label", default_zoom_lbl)
        else:
            self.lbl_zoom = default_zoom_lbl
        self.lbl_bf = get_user_input("Brightfield Label", "BF") if self.do_bf else "BF"

        # 6. ILLUSTRATOR / PRINT SETTINGS -----------------------------------
        print("\n--- ILLUSTRATOR EXPORT SETTINGS ---")
        # panel_w_mm: printed width of one panel column in millimetres
        # sb_len_mm : physical length represented by the scale bar
        self.panel_w_mm  = float(get_user_input("Panel Width (mm)",        20.0))
        self.spacing_pt  = float(get_user_input("Spacing (pt)",             3.0))
        self.font_size   = float(get_user_input("Font Size (pt)",           5.0))
        self.sb_len_mm   = float(get_user_input("Scale Bar Length (mm)",    2.646))

        # 7. QUANTIFICATION -------------------------------------------------
        print("\n--- ANALYSIS MODULE ---")
        q_resp = get_user_input(
            "Enable Quantification (Pearson/Manders)? [y/n]", "y").lower()
        self.do_quant = (q_resp == 'y')

        # --- DERIVED VALUES ------------------------------------------------
        self.aspect_ratio   = self.crop_h / self.crop_w
        self.panel_h_mm     = self.panel_w_mm * self.aspect_ratio

        self.panel_w_inch   = self.panel_w_mm / 25.4
        self.panel_h_inch   = self.panel_h_mm / 25.4
        self.spacing_inch   = self.spacing_pt / 72.0

        self.dpi            = self.crop_w / self.panel_w_inch

        # Scale bar in pixel units
        self.sb_px_w   = (self.sb_len_mm / self.panel_w_mm) * self.crop_w
        self.sb_px_h   = (0.5            / self.panel_w_mm) * self.crop_w
        self.sb_margin = (1.0            / self.panel_w_mm) * self.crop_w

        # Matplotlib font / PDF settings for vector output
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['ps.fonttype']  = 42
        plt.rcParams['font.family']  = 'Arial'

        self.font_prop = FontProperties()
        self.font_prop.set_family('Arial')
        self.font_prop.set_size(self.font_size)
        self.font_prop.set_weight('bold')

        # Pseudocolour mixing constants (RGB)
        self.MIX_CYAN    = np.array([0.0, 1.0, 1.0])
        self.MIX_GREEN   = np.array([0.0, 1.0, 0.0])
        self.MIX_MAGENTA = np.array([1.0, 0.0, 1.0])

        # Zoom upscale: zoom_size * zoom_mag == crop_w  (square)
        self.zoom_upscaled_size = self.crop_w if self.do_zoom else 0
        # Physical width of the zoom panel = panel height  (square in print)
        self.zoom_panel_w_inch  = self.panel_h_inch if self.do_zoom else 0

        # --- SUMMARY -------------------------------------------------------
        print("\n>>> CONFIGURATION COMPLETE")
        if self.do_zoom:
            print(f"    Zoom ROI    : {self.zoom_size}x{self.zoom_size}px (square) -> "
                  f"{self.zoom_upscaled_size}x{self.zoom_upscaled_size}px "
                  f"after {self.zoom_mag}x bicubic upscale")
            print(f"    Zoom panel  : {self.panel_h_mm:.1f}x{self.panel_h_mm:.1f} mm "
                  f"(square, fits row height)")

        print("\n" + "="*50)
        print("  HOW TO USE IN NAPARI")
        print("="*50)
        print("  1. Drag the YELLOW box to your region of interest")
        print("  2. (If zoom enabled) Drag the CYAN box to the")
        print("     sub-region you want enlarged")
        print("  3. Press  [c]  ->  Crop, save & advance to next image")
        print("  4. Press  [s]  ->  Skip this image & advance")
        print("  5. After the last image, montage + CSV are auto-saved")
        print("="*50 + "\n")


# ================================================================
#  RGB CHANNEL DETECTOR
# ================================================================
class RGBDetector:
    """Classify a single-channel TIFF into a fluorophore role based on its
    RGB colour statistics (works with false-colour TIFF exports from most
    confocal acquisition software).

    Returns one of: ``'cyan'``, ``'green'``, ``'mag'``, ``'bf'``,
    ``'empty'``, ``'unknown'``, or ``'error'``.
    """

    def detect(self, path):
        try:
            img = tifffile.imread(path)
            if img.ndim != 3:
                return "unknown"
            # Normalise to channel-last (HxWx3)
            if img.shape[0] == 3:
                img = np.moveaxis(img, 0, -1)

            p99_r = np.percentile(img[..., 0], 99.5)
            p99_g = np.percentile(img[..., 1], 99.5)
            p99_b = np.percentile(img[..., 2], 99.5)

            max_sig = max(p99_r, p99_g, p99_b)
            if max_sig < 1.0:
                return "empty"

            total = p99_r + p99_g + p99_b
            if total == 0:
                return "empty"

            r_score = p99_r / total
            g_score = p99_g / total
            b_score = p99_b / total

            # Brightfield: roughly equal R/G/B
            if (abs(r_score - g_score) < 0.2 and
                    abs(g_score - b_score) < 0.2):
                return "bf"
            if g_score > 0.45:
                return "green"
            if r_score > 0.3 and b_score > 0.3:
                return "mag"
            if r_score > 0.7:
                return "mag"
            if g_score > 0.3 and b_score > 0.3:
                return "cyan"
            return "unknown"

        except Exception:
            return "error"


# ================================================================
#  COLOCALIZATION QUANTIFIER
# ================================================================
class Quantifier:
    """Compute Pearson's R and Manders' coefficients (M1 / M2) for a
    pair of normalised single-channel arrays.

    Both input arrays should be normalised to [0, 1] before calling
    :meth:`analyze`.
    """

    @staticmethod
    def analyze(name, green_crop, mag_crop):
        """
        Parameters
        ----------
        name : str
            Image identifier (used as ``Image_ID`` in the output dict).
        green_crop, mag_crop : ndarray
            2-D float arrays normalised to [0, 1].

        Returns
        -------
        stats : dict or None
        mask_g : ndarray[bool] or None   -- Otsu mask for *green*
        mask_m : ndarray[bool] or None   -- Otsu mask for *magenta*
        """
        flat_g = green_crop.flatten()
        flat_m = mag_crop.flatten()
        if len(flat_g) < 2:
            return None, None, None

        pearson_r, _ = pearsonr(flat_g, flat_m)

        try:
            thresh_g = threshold_otsu(green_crop)
            thresh_m = threshold_otsu(mag_crop)
        except ValueError:
            thresh_g = 0.01
            thresh_m = 0.01

        mask_g = green_crop > thresh_g
        mask_m = mag_crop  > thresh_m

        total_g = np.sum(green_crop)
        total_m = np.sum(mag_crop)

        m1 = np.sum(green_crop[mask_m]) / total_g if total_g != 0 else 0
        m2 = np.sum(mag_crop[mask_g])   / total_m if total_m != 0 else 0

        stats = {
            'Image_ID':              name,
            'Pearson_R':             round(pearson_r, 3),
            'Manders_M1 (G_in_M)':  round(m1, 3),
            'Manders_M2 (M_in_G)':  round(m2, 3),
            'Thresh_Green':          round(thresh_g, 4),
            'Thresh_Mag':            round(thresh_m, 4),
        }
        return stats, mask_g, mask_m


# ================================================================
#  MAIN MANAGER
# ================================================================
class ColocManager:
    """Top-level controller that drives the Napari viewer, handles key
    bindings, orchestrates cropping / quantification, and writes all output.
    """

    def __init__(self):
        self.cfg       = SessionConfig()
        self.detector  = RGBDetector()
        self.gallery   = []
        self.stats_log = []

        self.experiments = self.scan_directory()
        if not self.experiments:
            print("CRITICAL ERROR: No TIF files found in the input folder.")
            return

        self.index             = 0
        self.viewer            = napari.Viewer()
        self.shapes_layer      = None
        self.zoom_shapes_layer = None

        print(f"Found {len(self.experiments)} experiment sets.")
        print("   Keys:  [c] = Crop & Next  |  [s] = Skip\n")

        self.viewer.bind_key('c', self.process_and_next)
        self.viewer.bind_key('s', self.skip_and_next)
        self.load_current_set()
        napari.run()

    # ------------------------------------------------------------------
    #  Directory scanning
    # ------------------------------------------------------------------
    def scan_directory(self):
        """Return a sorted list of unique base paths (everything before '_ch')."""
        search_pattern = os.path.join(self.cfg.input_dir, "**", "*.tif")
        all_tifs       = glob.glob(search_pattern, recursive=True)
        unique_bases   = set()

        for f in all_tifs:
            if os.path.basename(f).startswith("._"):
                continue   # skip macOS resource forks
            if "_ch" in f:
                base = f.split("_ch")[0]
                unique_bases.add(base)

        return sorted(unique_bases)

    # ------------------------------------------------------------------
    #  Image loading helpers
    # ------------------------------------------------------------------
    def load_flat(self, path):
        """Read a TIFF and collapse to a 2-D (HxW) float array via max projection."""
        if not path or not os.path.exists(path):
            return None
        raw = tifffile.imread(path)
        if raw.ndim == 2:
            return raw
        if raw.ndim == 3:
            if raw.shape[-1] == 3:
                return np.max(raw, axis=-1)
            return np.max(raw, axis=0)
        return raw

    def auto_contrast(self, img):
        """Stretch contrast to the 0.1-99.9 percentile range -> [0, 1]."""
        if img is None or img.size == 0:
            return None
        vmin, vmax = np.percentile(img, 0.1), np.percentile(img, 99.9)
        if vmax <= vmin:
            vmax = vmin + 1
        return np.clip((img.astype(float) - vmin) / (vmax - vmin), 0, 1)

    def make_rgb(self, norm, color_mix):
        """Multiply a greyscale [0, 1] image by an RGB colour triplet."""
        if norm is None:
            return None
        return norm[..., np.newaxis] * color_mix

    # ------------------------------------------------------------------
    #  Napari shapes helpers
    # ------------------------------------------------------------------
    def add_default_boxes(self):
        """Place the yellow crop box (and optional cyan zoom box) at image centre."""
        h, w = 1024, 1024
        for layer in self.viewer.layers:
            if isinstance(layer, napari.layers.Image) and layer.data.ndim >= 2:
                h, w = layer.data.shape[-2], layer.data.shape[-1]
                break

        cy, cx = h // 2, w // 2
        rh, rw = self.cfg.crop_h // 2, self.cfg.crop_w // 2

        self.shapes_layer.data = [
            np.array([[cy-rh, cx-rw], [cy-rh, cx+rw],
                      [cy+rh, cx+rw], [cy+rh, cx-rw]])
        ]
        self.shapes_layer.mode = 'select'

        if self.cfg.do_zoom and self.zoom_shapes_layer:
            zh = self.cfg.zoom_size // 2   # square half-side
            self.zoom_shapes_layer.data = [
                np.array([[cy-zh, cx-zh], [cy-zh, cx+zh],
                          [cy+zh, cx+zh], [cy+zh, cx-zh]])
            ]
            self.zoom_shapes_layer.mode = 'select'

        self.viewer.layers.selection.active = self.shapes_layer

    def get_box_center(self, shape_data):
        """Return (y_centre, x_centre) of a quadrilateral shape."""
        y_vals = [p[0] for p in shape_data]
        x_vals = [p[1] for p in shape_data]
        return (int((min(y_vals) + max(y_vals)) / 2),
                int((min(x_vals) + max(x_vals)) / 2))

    # ------------------------------------------------------------------
    #  Image loading into Napari
    # ------------------------------------------------------------------
    def load_current_set(self):
        """Load the next experiment set into the Napari viewer."""
        if self.index >= len(self.experiments):
            print(">> ALL DATASETS PROCESSED.")
            self.finalize_analysis()
            return

        self.current_base_path = self.experiments[self.index]
        self.current_name      = os.path.basename(self.current_base_path)

        files     = sorted(glob.glob(f"{self.current_base_path}_ch*.tif"))
        self.map  = {'cyan': None, 'green': None, 'mag': None, 'bf': None}

        print(f"\n[{self.index+1}/{len(self.experiments)}] {self.current_name}")

        for f in files:
            role = self.detector.detect(f)
            if role in self.map and self.map[role] is None:
                self.map[role] = f

        self.mode = "TRIPLE" if self.map['cyan'] else "DUAL"

        # Clear previous layers
        self.viewer.layers.select_all()
        self.viewer.layers.remove_selected()

        # Add image layers
        if self.map['cyan']:
            self.viewer.add_image(self.load_flat(self.map['cyan']),
                name='Cyan_Ch',  colormap='cyan',     blending='additive')
        if self.map['green']:
            self.viewer.add_image(self.load_flat(self.map['green']),
                name='Green_Ch', colormap='green',    blending='additive')
        if self.map['mag']:
            self.viewer.add_image(self.load_flat(self.map['mag']),
                name='Mag_Ch',   colormap='magenta',  blending='additive')
        if self.cfg.do_bf and self.map['bf']:
            self.viewer.add_image(self.load_flat(self.map['bf']),
                name='BF_Ch',    colormap='gray',
                blending='translucent', opacity=0.4)

        # Add shapes layers
        self.shapes_layer = self.viewer.add_shapes(
            name='1_MAIN_CROP',
            edge_color='yellow', face_color='transparent', edge_width=3)
        if self.cfg.do_zoom:
            self.zoom_shapes_layer = self.viewer.add_shapes(
                name='2_ZOOM_ROI',
                edge_color='cyan', face_color='transparent', edge_width=2)

        self.add_default_boxes()

    # ------------------------------------------------------------------
    #  Key-binding callbacks
    # ------------------------------------------------------------------
    def skip_and_next(self, viewer):
        """[s] Skip the current image and move to the next."""
        print(f"   >> Skipping {self.current_name}...")
        self.index += 1
        self.load_current_set()

    def process_and_next(self, viewer):
        """[c] Crop the current image, run analysis, save outputs, and advance."""
        if self.shapes_layer is None or len(self.shapes_layer.data) == 0:
            print("   >> No crop box found!")
            return

        yc, xc = self.get_box_center(self.shapes_layer.data[-1])
        rh, rw  = self.cfg.crop_h // 2, self.cfg.crop_w // 2
        y1, y2  = yc - rh, yc + rh
        x1, x2  = xc - rw, xc + rw

        try:
            def get_layer_data(name):
                return viewer.layers[name].data if name in viewer.layers else None

            img_c  = get_layer_data('Cyan_Ch')
            img_g  = get_layer_data('Green_Ch')
            img_m  = get_layer_data('Mag_Ch')
            img_bf = get_layer_data('BF_Ch')

            raw_crop_c  = self.safe_crop(img_c,  y1, y2, x1, x2)
            raw_crop_g  = self.safe_crop(img_g,  y1, y2, x1, x2)
            raw_crop_m  = self.safe_crop(img_m,  y1, y2, x1, x2)
            raw_crop_bf = self.safe_crop(img_bf, y1, y2, x1, x2)

            h, w = self.cfg.crop_h, self.cfg.crop_w

            def ensure_shape(arr):
                if arr is None or arr.shape != (h, w):
                    return np.zeros((h, w), dtype=np.float32)
                return arr.astype(np.float32)

            # --- QUANTIFICATION -------------------------------------------
            mask_g = mask_m = None
            if self.cfg.do_quant:
                def norm_for_stats(arr):
                    mmax = np.max(arr)
                    return arr / mmax if mmax != 0 else arr

                stats_g = norm_for_stats(ensure_shape(raw_crop_g))
                stats_m = norm_for_stats(ensure_shape(raw_crop_m))

                stats, mask_g, mask_m = Quantifier.analyze(
                    self.current_name, stats_g, stats_m)
                self.stats_log.append(stats)
                print(f"   >> STATS: R={stats['Pearson_R']}, "
                      f"M1={stats['Manders_M1 (G_in_M)']}")
            else:
                print("   >> Visuals only")

            # --- DISPLAY NORMALISATION ------------------------------------
            norm_c  = ensure_shape(self.auto_contrast(raw_crop_c))
            norm_g  = ensure_shape(self.auto_contrast(raw_crop_g))
            norm_m  = ensure_shape(self.auto_contrast(raw_crop_m))
            norm_bf = ensure_shape(self.auto_contrast(raw_crop_bf))

            c_rgb     = self.make_rgb(norm_c, self.cfg.MIX_CYAN)
            g_rgb     = self.make_rgb(norm_g, self.cfg.MIX_GREEN)
            m_rgb     = self.make_rgb(norm_m, self.cfg.MIX_MAGENTA)
            merge_rgb = np.clip(c_rgb + g_rgb + m_rgb, 0, 1)
            bf_rgb    = np.stack((norm_bf,)*3, axis=-1)

            # --- ZOOM PANEL -----------------------------------------------
            zoom_rgb            = None
            zoom_coords_in_main = None

            if (self.cfg.do_zoom and self.zoom_shapes_layer is not None
                    and len(self.zoom_shapes_layer.data) > 0):

                zyc, zxc = self.get_box_center(self.zoom_shapes_layer.data[-1])
                zs       = self.cfg.zoom_size    # square side
                zh       = zs // 2
                zy1, zy2 = zyc - zh, zyc + zh
                zx1, zx2 = zxc - zh, zxc + zh
                target_px = self.cfg.zoom_upscaled_size  # = crop_w

                def zoom_crop_upscale_contrast(img_data):
                    crop = self.safe_crop(img_data, zy1, zy2, zx1, zx2)
                    if crop is None or crop.shape != (zs, zs):
                        return None
                    upscaled = self.upscale_channel(crop, target_px, target_px)
                    return self.auto_contrast(upscaled)

                zn_c = zoom_crop_upscale_contrast(img_c)
                zn_g = zoom_crop_upscale_contrast(img_g)
                zn_m = zoom_crop_upscale_contrast(img_m)

                z_rgb = np.zeros((target_px, target_px, 3))
                if zn_c is not None: z_rgb += self.make_rgb(zn_c, self.cfg.MIX_CYAN)
                if zn_g is not None: z_rgb += self.make_rgb(zn_g, self.cfg.MIX_GREEN)
                if zn_m is not None: z_rgb += self.make_rgb(zn_m, self.cfg.MIX_MAGENTA)

                zoom_rgb            = np.clip(z_rgb, 0, 1)
                zoom_coords_in_main = (zy1 - y1, zx1 - x1)

                print(f"   >> ZOOM: {zs}x{zs} -> {target_px}x{target_px} "
                      f"({self.cfg.zoom_mag}x bicubic upscale)")

            # --- SAVE -----------------------------------------------------
            self.save_individual_and_local_panel(
                self.current_name,
                c_rgb, g_rgb, m_rgb, merge_rgb, bf_rgb,
                zoom_rgb, zoom_coords_in_main, mask_g, mask_m)

            self.gallery.append({
                'name':        self.current_name,
                'cyan':        c_rgb,
                'green':       g_rgb,
                'mag':         m_rgb,
                'merge':       merge_rgb,
                'bf':          bf_rgb if self.cfg.do_bf else None,
                'zoom':        zoom_rgb,
                'zoom_coords': zoom_coords_in_main,
                'mode':        self.mode,
            })

            self.index += 1
            self.load_current_set()

        except Exception as e:
            print(f"ERROR: {e}")
            import traceback; traceback.print_exc()

    # ------------------------------------------------------------------
    #  Cropping / upscaling utilities
    # ------------------------------------------------------------------
    def safe_crop(self, img, y1, y2, x1, x2):
        """Clip crop coordinates to image bounds and return the sub-array."""
        if img is None:
            return None
        h, w = img.shape[:2]
        y1, x1 = max(0, y1), max(0, x1)
        y2, x2 = min(h, y2), min(w, x2)
        if y1 >= y2 or x1 >= x2:
            return None
        return img[y1:y2, x1:x2]

    def upscale_channel(self, crop, target_h, target_w):
        """Bicubic upscale a 2-D array to (*target_h*, *target_w*)."""
        if crop is None:
            return None
        return resize(
            crop.astype(np.float64), (target_h, target_w),
            order=3, anti_aliasing=False, preserve_range=True)

    # ------------------------------------------------------------------
    #  Figure layout helpers
    # ------------------------------------------------------------------
    def _col_widths_inch(self, col_labels):
        """Return a list of column widths (inches).
        The zoom column is square -> ``zoom_panel_w_inch``.
        All other columns -> ``panel_w_inch``."""
        return [self.cfg.zoom_panel_w_inch
                if lbl == self.cfg.lbl_zoom else self.cfg.panel_w_inch
                for lbl in col_labels]

    def _make_axes_grid(self, fig, n_rows, n_cols, col_w, total_w, total_h):
        """Create a 2-D list of axes using ``fig.add_axes`` for
        pixel-perfect placement with per-column widths."""
        row_h = self.cfg.panel_h_inch
        axes  = []
        for r in range(n_rows):
            row_top = r * (row_h + self.cfg.spacing_inch)
            bottom  = 1.0 - (row_top + row_h) / total_h
            h_frac  = row_h / total_h

            row_axes = []
            x = 0.0
            for ci in range(n_cols):
                left   = x / total_w
                w_frac = col_w[ci] / total_w
                ax     = fig.add_axes([left, bottom, w_frac, h_frac])
                row_axes.append(ax)
                x += col_w[ci] + self.cfg.spacing_inch
            axes.append(row_axes)
        return axes

    def add_scale_bar_patch(self, ax):
        """Draw a white filled rectangle as a scale bar on *ax*."""
        x_start = self.cfg.crop_w - self.cfg.sb_px_w - self.cfg.sb_margin
        y_start = self.cfg.crop_h - self.cfg.sb_margin - self.cfg.sb_px_h
        ax.add_patch(patches.Rectangle(
            (x_start, y_start), self.cfg.sb_px_w, self.cfg.sb_px_h,
            linewidth=0, edgecolor='none', facecolor='white'))

    # ------------------------------------------------------------------
    #  Output saving
    # ------------------------------------------------------------------
    def save_individual_and_local_panel(self, name, c, g, m, merge, bf,
                                        zoom, zoom_coords, mask_g, mask_m):
        """Save individual channel TIFFs, QC masks, and a panel PDF/PNG
        for a single cropped field of view."""
        save_path = os.path.join(self.cfg.output_dir, name)
        os.makedirs(save_path, exist_ok=True)

        def to_uint8(arr):
            return (arr * 255).astype(np.uint8)

        # Individual TIFFs
        if self.mode == "TRIPLE":
            tifffile.imwrite(os.path.join(save_path, "1_Cyan.tif"),    to_uint8(c))
        tifffile.imwrite(os.path.join(save_path, "2_Green.tif"),       to_uint8(g))
        tifffile.imwrite(os.path.join(save_path, "3_Magenta.tif"),     to_uint8(m))
        tifffile.imwrite(os.path.join(save_path, "4_Merge.tif"),       to_uint8(merge))
        if zoom is not None:
            tifffile.imwrite(os.path.join(save_path, "5_Zoom.tif"),    to_uint8(zoom))
        if self.cfg.do_bf:
            tifffile.imwrite(os.path.join(save_path, "6_BF.tif"),      to_uint8(bf))

        # QC Otsu masks
        if mask_g is not None and mask_m is not None:
            qc_path = os.path.join(save_path, "QC_Masks")
            os.makedirs(qc_path, exist_ok=True)
            tifffile.imwrite(
                os.path.join(qc_path, "Mask_Green_Otsu.tif"), to_uint8(mask_g))
            tifffile.imwrite(
                os.path.join(qc_path, "Mask_Mag_Otsu.tif"),   to_uint8(mask_m))

        # Build panel plot list (BF always last)
        plot_list = []
        if self.mode == "TRIPLE":
            plot_list.append((c,     self.cfg.lbl_c))
        plot_list.append((g,         self.cfg.lbl_g))
        plot_list.append((m,         self.cfg.lbl_m))
        plot_list.append((merge,     self.cfg.lbl_merge))
        if zoom is not None:
            plot_list.append((zoom,  self.cfg.lbl_zoom))
        if self.cfg.do_bf:
            plot_list.append((bf,    self.cfg.lbl_bf))

        n_cols     = len(plot_list)
        col_labels = [lbl for _, lbl in plot_list]
        col_w      = self._col_widths_inch(col_labels)

        total_w = sum(col_w) + (n_cols - 1) * self.cfg.spacing_inch
        total_h = self.cfg.panel_h_inch

        fig        = plt.figure(figsize=(total_w, total_h))
        axes_grid  = self._make_axes_grid(fig, 1, n_cols, col_w, total_w, total_h)
        axes       = axes_grid[0]

        for i, ax in enumerate(axes):
            img_data, label = plot_list[i]
            ax.imshow(img_data, aspect='auto')
            ax.axis('off')
            ax.text(0.05, 0.95, label, transform=ax.transAxes, color='white',
                    fontproperties=self.cfg.font_prop, ha='left', va='top')
            if i == 0:
                self.add_scale_bar_patch(ax)

            # Dashed outline on the merge panel showing the zoom ROI
            if label == self.cfg.lbl_merge and zoom_coords is not None:
                ry, rx = zoom_coords
                ax.add_patch(patches.Rectangle(
                    (rx, ry), self.cfg.zoom_size, self.cfg.zoom_size,
                    linewidth=1, edgecolor='white',
                    facecolor='none', linestyle='--'))

        plt.savefig(os.path.join(save_path, "Panel_View.pdf"), dpi=self.cfg.dpi)
        plt.savefig(os.path.join(save_path, "Panel_View.png"), dpi=self.cfg.dpi)
        plt.close(fig)

    # ------------------------------------------------------------------
    #  Finalisation
    # ------------------------------------------------------------------
    def finalize_analysis(self):
        """Called after the last image is processed: build montage and save CSV."""
        if self.gallery:
            self.build_global_montage()

        if self.cfg.do_quant and self.stats_log:
            df       = pd.DataFrame(self.stats_log)
            csv_path = os.path.join(self.cfg.output_dir,
                                    f"{self.cfg.exp_name}_QUANTIFICATION.csv")
            df.to_csv(csv_path, index=False)
            print(f">> STATS SAVED: {csv_path}")

    # ------------------------------------------------------------------
    #  Summary Montage
    # ------------------------------------------------------------------
    def build_global_montage(self):
        """Assemble all cropped panels into a single summary montage (PDF + PNG)."""
        print("\n>>> GENERATING SUMMARY MONTAGES (PDF + PNG)...")

        n_rows   = len(self.gallery)
        has_cyan = any(x['mode'] == "TRIPLE"    for x in self.gallery)
        has_zoom = any(x['zoom'] is not None     for x in self.gallery)

        col_keys   = []
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

        n_cols  = len(col_keys)
        col_w   = self._col_widths_inch(col_labels)
        total_w = sum(col_w) + (n_cols - 1) * self.cfg.spacing_inch
        total_h = (n_rows * self.cfg.panel_h_inch +
                   (n_rows - 1) * self.cfg.spacing_inch)

        fig  = plt.figure(figsize=(total_w, total_h))
        axes = self._make_axes_grid(fig, n_rows, n_cols, col_w, total_w, total_h)

        h, w        = self.cfg.crop_h, self.cfg.crop_w
        black_main  = np.zeros((h, w, 3))
        zup         = self.cfg.zoom_upscaled_size if self.cfg.zoom_upscaled_size > 0 else h
        black_zoom  = np.zeros((zup, zup, 3))

        for r, item in enumerate(self.gallery):
            for ci, key in enumerate(col_keys):
                ax       = axes[r][ci]
                img_data = item.get(key)
                if img_data is None:
                    img_data = black_zoom if key == 'zoom' else black_main

                ax.imshow(img_data, aspect='auto')
                ax.axis('off')

                # Column headers on the first row only
                if r == 0:
                    ax.text(0.05, 0.95, col_labels[ci], transform=ax.transAxes,
                            color='white', fontproperties=self.cfg.font_prop,
                            ha='left', va='top')

                # Scale bar on the top-left panel only
                if r == 0 and ci == 0:
                    self.add_scale_bar_patch(ax)

                # Dashed zoom outline on merge column
                if key == 'merge' and item.get('zoom_coords') is not None:
                    ry, rx = item['zoom_coords']
                    ax.add_patch(patches.Rectangle(
                        (rx, ry), self.cfg.zoom_size, self.cfg.zoom_size,
                        linewidth=0.5, edgecolor='white',
                        facecolor='none', linestyle='--'))

        base_name = os.path.join(self.cfg.output_dir,
                                 f"{self.cfg.exp_name}_SUMMARY_MONTAGE")
        plt.savefig(f"{base_name}.pdf", dpi=self.cfg.dpi)
        plt.savefig(f"{base_name}.png", dpi=self.cfg.dpi)
        print(f">> MONTAGE SAVED: {base_name}.pdf")
        plt.close(fig)


# ================================================================
#  ENTRY POINT
# ================================================================
if __name__ == "__main__":
    ColocManager()
