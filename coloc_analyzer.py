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
named with a ``_ch<N>.tif`` suffix, alongside their exported XML metadata.

Output (per image)
------------------
``<output_dir>/<experiment_name>/<image_basename>/``
"""

import napari
import os
import glob
import numpy as np
import re
import xml.etree.ElementTree as ET

# imagecodecs is required by tifffile to decompress LZW / Deflate (ZIP)
# compressed TIFFs, which is the default lossless export format in Leica LAS X.
try:
    import imagecodecs  # noqa: F401
except ImportError:
    print(
        "\nCRITICAL ERROR: 'imagecodecs' is not installed.\n"
        "  tifffile requires this package to read LZW / Deflate-compressed\n"
        "  TIFFs (the lossless export format used by Leica LAS X).\n"
        "  Install it with:\n"
        "      pip install imagecodecs\n"
    )
    exit(1)
import pickle
import math
from PIL import Image
import pandas as pd
import tifffile
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.font_manager import FontProperties
from skimage.transform import resize
from skimage.measure import profile_line
from scipy.stats import pearsonr
import datetime
import warnings

# ================================================================
#  DEFAULT BASE DIRECTORY
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
        today = datetime.datetime.now().strftime("%Y%m%d")
        self.date_folder = get_user_input("Date Folder (YYYYMMDD)", today)

        default_exp = f"Experiment_{self.date_folder}"
        self.exp_name = get_user_input("Experiment Name", default_exp)

        base_out = os.path.join(BASE_DIR, "Output_Coloc")
        self.date_dir = os.path.join(base_out, self.date_folder)
        self.output_dir = os.path.join(self.date_dir, self.exp_name)
        os.makedirs(self.output_dir, exist_ok=True)

        # 3. CROP GEOMETRY --------------------------------------------------
        print("\n--- CROP GEOMETRY ---")
        self.crop_w = int(get_user_input("Crop Width (px)", 360))
        self.crop_h = int(get_user_input("Crop Height (px)", 360))

        # 4. OPTIONAL FEATURES ----------------------------------------------
        print("\n--- EXTRA FEATURES ---")
        self.do_bf = get_user_input(
            "Include Brightfield Panel?[y/n]", "n").lower() == 'y'
        self.do_zoom = get_user_input(
            "Include Zoom/Enlargement Panel? [y/n]", "n").lower() == 'y'
        self.zoom_mag  = 1
        self.zoom_size = 0   

        if self.do_zoom:
            mag = int(get_user_input("   > Magnification Factor (e.g. 3)", "3"))
            self.zoom_mag  = mag
            self.zoom_size = self.crop_w // mag   

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
        self.panel_w_mm  = float(get_user_input("Panel Width (mm)",        20.0))
        self.spacing_pt  = float(get_user_input("Spacing (pt)",             3.0))
        self.font_size   = float(get_user_input("Font Size (pt)",           5.0))
        self.sb_len_mm   = float(get_user_input("Scale Bar Length (mm)",    2.646))

        # 7. QUANTIFICATION -------------------------------------------------
        print("\n--- ANALYSIS MODULE ---")
        q_resp = get_user_input(
            "Enable Quantification (Pearson/Manders)? [y/n]", "y").lower()
        self.do_quant = (q_resp == 'y')

        # 8. INTENSITY LINE PROFILE -----------------------------------------
        ip_resp = get_user_input(
            "Enable Intensity Line Profile? [y/n]", "n").lower()
        self.do_intensity_profile = (ip_resp == 'y')

        # --- DERIVED VALUES ------------------------------------------------
        self.aspect_ratio   = self.crop_h / self.crop_w
        self.panel_h_mm     = self.panel_w_mm * self.aspect_ratio

        self.panel_w_inch   = self.panel_w_mm / 25.4
        self.panel_h_inch   = self.panel_h_mm / 25.4
        self.spacing_inch   = self.spacing_pt / 72.0

        self.dpi            = self.crop_w / self.panel_w_inch

        self.sb_px_w   = (self.sb_len_mm / self.panel_w_mm) * self.crop_w
        self.sb_px_h   = (0.5            / self.panel_w_mm) * self.crop_w
        self.sb_margin = (1.0            / self.panel_w_mm) * self.crop_w

        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['ps.fonttype']  = 42
        plt.rcParams['font.family']  = 'Arial'

        self.font_prop = FontProperties()
        self.font_prop.set_family('Arial')
        self.font_prop.set_size(self.font_size)
        self.font_prop.set_weight('bold')

        self.MIX_CYAN    = np.array([0.0, 1.0, 1.0])
        self.MIX_GREEN   = np.array([0.0, 1.0, 0.0])
        self.MIX_MAGENTA = np.array([1.0, 0.0, 1.0])

        self.zoom_upscaled_size = self.crop_w if self.do_zoom else 0
        self.zoom_panel_w_inch  = self.panel_h_inch if self.do_zoom else 0

        print("\n>>> CONFIGURATION COMPLETE")
        if self.do_zoom:
            print(f"    Zoom ROI    : {self.zoom_size}×{self.zoom_size}px (square) → "
                  f"{self.zoom_upscaled_size}×{self.zoom_upscaled_size}px "
                  f"after {self.zoom_mag}× bicubic upscale")
            print(f"    Zoom panel  : {self.panel_h_mm:.1f}×{self.panel_h_mm:.1f} mm "
                  f"(square, fits row height)")

        print("\n" + "="*50)
        print("  HOW TO USE IN NAPARI")
        print("="*50)
        print("  1. Drag the YELLOW box to your region of interest")
        print("  2. (If zoom enabled) Drag the CYAN box to the")
        print("     sub-region you want enlarged")
        if self.do_intensity_profile:
            print("  3. Drag the WHITE line for intensity profiling")
            print("  4. Press  [c]  →  Crop, save & advance to next image")
            print("  5. Press  [s]  →  Skip this image & advance")
            print("  6. After the last image, montage + CSV are auto-saved")
        else:
            print("  3. Press  [c]  →  Crop, save & advance to next image")
            print("  4. Press  [s]  →  Skip this image & advance")
            print("  5. After the last image, montage + CSV are auto-saved")
        print("="*50 + "\n")


# ================================================================
#  XML METADATA DETECTOR
# ================================================================
class XMLMetadataDetector:
    """Parses Leica LAS X exported XML metadata to map channel indices to roles.
    Requires _chXX.tif nomenclature and a corresponding .xml file.
    """
    def __init__(self):
        # Strict mapping of Leica LUT naming conventions to script roles
        self.lut_map = {
            'cyan':    'cyan',
            'blue':    'cyan',   # Sometimes CFP/Cyan is exported as Blue LUT
            'green':   'green',
            'magenta': 'mag',
            'red':     'mag',    # mCherry/TRITC often exported as Red LUT
            'gray':    'bf',
            'grey':    'bf'
        }

    def assign_roles(self, files):
        """Map a list of TIF files to their respective roles using Leica XML."""
        result = {'cyan': None, 'green': None, 'mag': None, 'bf': None}
        
        for tif_path in files:
            base_dir = os.path.dirname(tif_path)
            filename = os.path.basename(tif_path)
            
            # Extract base name and channel index
            match = re.search(r'(.*)_ch(\d+)\.tif$', filename, re.IGNORECASE)
            if not match:
                continue
                
            base_name = match.group(1)
            ch_idx    = int(match.group(2))
            
            # Locate the corresponding XML file
            xml_candidates =[
                os.path.join(base_dir, f"{base_name}_Properties.xml"),
                os.path.join(base_dir, f"{base_name}.xml"),
                os.path.join(base_dir, "MetaData", f"{base_name}_Properties.xml"),
                os.path.join(base_dir, "MetaData", f"{base_name}.xml")
            ]
            
            xml_path = next((p for p in xml_candidates if os.path.exists(p)), None)
            
            if not xml_path:
                print(f"   >> WARNING: No XML metadata found for {filename}.")
                continue
                
            # Parse XML and map LUT
            try:
                tree = ET.parse(xml_path)
                root = tree.getroot()
                
                channels_node = root.find('.//Channels')
                if channels_node is None:
                    print(f"   >> ERROR: <Channels> tag missing in {os.path.basename(xml_path)}")
                    continue
                    
                channels = channels_node.findall('./ChannelDescription')
                
                if ch_idx >= len(channels):
                    print(f"   >> ERROR: Index {ch_idx} exceeds XML channel count in {os.path.basename(xml_path)}.")
                    continue
                    
                lut_name = channels[ch_idx].attrib.get('LUTName', '').lower()
                role = self.lut_map.get(lut_name)
                
                if role in result:
                    result[role] = tif_path
                    print(f"   >> METADATA: Mapped {filename} (LUT: {lut_name.capitalize()}) -> {role.upper()}")
                else:
                    print(f"   >> UNKNOWN LUT: {filename} has unmapped LUT '{lut_name}'. Ignored.")
                    
            except Exception as e:
                print(f"   >> ERROR parsing XML {xml_path}: {e}")

        return result


# ================================================================
#  COLOCALIZATION QUANTIFIER
# ================================================================
def _threshold_max_entropy(image):
    """Kapur Maximum-Entropy thresholding (Kapur et al., 1985).

    Selects the threshold that maximises the sum of foreground and
    background entropy, making it more sensitive to sparse puncta
    than Otsu's bimodal assumption.
    """
    # Build normalised 256-bin histogram
    hist, bin_edges = np.histogram(image.ravel(), bins=256, density=False)
    hist = hist.astype(np.float64)
    hist = hist / hist.sum()          # probability distribution

    # Cumulative sums
    cum_sum = np.cumsum(hist)
    # Avoid log(0)
    eps = np.finfo(np.float64).eps

    best_t   = 0
    best_ent = -np.inf

    for t in range(1, 256):
        # Background probability
        p_bg = cum_sum[t - 1]
        if p_bg < eps or p_bg > 1.0 - eps:
            continue
        p_fg = 1.0 - p_bg

        # Background entropy
        h_bg = hist[:t]
        h_bg = h_bg[h_bg > eps]
        ent_bg = -np.sum((h_bg / p_bg) * np.log(h_bg / p_bg))

        # Foreground entropy
        h_fg = hist[t:]
        h_fg = h_fg[h_fg > eps]
        ent_fg = -np.sum((h_fg / p_fg) * np.log(h_fg / p_fg))

        total_ent = ent_bg + ent_fg
        if total_ent > best_ent:
            best_ent = total_ent
            best_t   = t

    # Map bin index back to intensity value
    return (bin_edges[best_t] + bin_edges[best_t + 1]) / 2.0


class Quantifier:
    """Compute Pearson's R and Manders' coefficients (M1 / M2) for a
    pair of normalised single-channel arrays."""

    @staticmethod
    def analyze(name, green_crop, mag_crop):
        flat_g = green_crop.flatten()
        flat_m = mag_crop.flatten()
        if len(flat_g) < 2:
            return None, None, None

        if np.std(flat_g) == 0 or np.std(flat_m) == 0:
            print("   >> WARNING: One or both channels are constant. Pearson R set to NaN.")
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
        # Initialize the new metadata detector
        self.detector  = XMLMetadataDetector()
        self.gallery   =[]
        self.stats_log =[]

        self.experiments = self.scan_directory()
        if not self.experiments:
            print("CRITICAL ERROR: No TIF files found in the input folder.")
            return

        self.index             = 0
        self.viewer            = napari.Viewer()
        self.shapes_layer      = None
        self.zoom_shapes_layer = None
        self.line_shapes_layer = None

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
        search_pattern = os.path.join(self.cfg.input_dir, "**", "*.tif")
        all_tifs       = glob.glob(search_pattern, recursive=True)
        unique_bases   = set()

        for f in all_tifs:
            if os.path.basename(f).startswith("._"):
                continue   
            if "_ch" in f:
                base = f.split("_ch")[0]
                unique_bases.add(base)

        return sorted(unique_bases)

    # ------------------------------------------------------------------
    #  Image loading helpers
    # ------------------------------------------------------------------
    def load_flat(self, path):
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

    def auto_contrast_limits(self, img):
        """Return (vmin, vmax) percentile limits for *img*."""
        if img is None or img.size == 0:
            return (0, 1)
        vmin, vmax = np.percentile(img, 0.1), np.percentile(img, 99.98)
        if vmax <= vmin:
            vmax = vmin + 1
        return (vmin, vmax)

    def auto_contrast(self, img, limits=None):
        if img is None or img.size == 0:
            return None
        if limits is not None:
            vmin, vmax = limits
        else:
            vmin, vmax = self.auto_contrast_limits(img)
        return np.clip((img.astype(float) - vmin) / (vmax - vmin), 0, 1)

    def make_rgb(self, norm, color_mix):
        if norm is None:
            return None
        return norm[..., np.newaxis] * color_mix

    # ------------------------------------------------------------------
    #  Napari shapes helpers
    # ------------------------------------------------------------------
    def add_default_boxes(self):
        h, w = 1024, 1024
        for layer in self.viewer.layers:
            if isinstance(layer, napari.layers.Image) and layer.data.ndim >= 2:
                h, w = layer.data.shape[-2], layer.data.shape[-1]
                break

        cy, cx = h // 2, w // 2
        rh, rw = self.cfg.crop_h // 2, self.cfg.crop_w // 2

        self.shapes_layer.data = [
            np.array([[cy-rh, cx-rw], [cy-rh, cx+rw],[cy+rh, cx+rw], [cy+rh, cx-rw]])
        ]
        self.shapes_layer.mode = 'select'

        if self.cfg.do_zoom and self.zoom_shapes_layer:
            zh = self.cfg.zoom_size // 2  
            self.zoom_shapes_layer.data = [
                np.array([[cy-zh, cx-zh], [cy-zh, cx+zh],
                          [cy+zh, cx+zh], [cy+zh, cx-zh]])
            ]
            self.zoom_shapes_layer.mode = 'select'

        if self.cfg.do_intensity_profile and self.line_shapes_layer is not None:
            self.line_shapes_layer.data = [
                np.array([[cy, cx - rw], [cy, cx + rw]])
            ]
            self.line_shapes_layer.mode = 'select'

        self.viewer.layers.selection.active = self.shapes_layer

    def get_box_center(self, shape_data):
        y_vals = [p[0] for p in shape_data]
        x_vals =[p[1] for p in shape_data]
        return (int((min(y_vals) + max(y_vals)) / 2),
                int((min(x_vals) + max(x_vals)) / 2))

    # ------------------------------------------------------------------
    #  Shape-layer sanitiser – prevent cross-layer drawing accidents
    # ------------------------------------------------------------------
    def _sanitize_shapes_layers(self):
        """Remove shapes that were accidentally drawn in the wrong layer.

        Rules
        -----
        * 1_MAIN_CROP  / 2_ZOOM_ROI  – only keep shapes with 4 vertices
          (rectangles).  Any 2-vertex lines are removed.
        * 3_LINE_PROFILE – only keep shapes with exactly 2 vertices
          (lines).  Any 4-vertex rectangles are removed.
        """
        def _filter(layer, keep_npts):
            if layer is None or len(layer.data) == 0:
                return
            cleaned = [s for s in layer.data if len(s) == keep_npts]
            removed = len(layer.data) - len(cleaned)
            if removed > 0:
                print(f"   >> WARNING: Removed {removed} stray shape(s) "
                      f"from '{layer.name}' (wrong vertex count).")
                layer.data = cleaned if cleaned else layer.data[:0]

        _filter(self.shapes_layer, 4)
        if self.cfg.do_zoom:
            _filter(self.zoom_shapes_layer, 4)
        if self.cfg.do_intensity_profile:
            _filter(self.line_shapes_layer, 2)

    # ------------------------------------------------------------------
    #  Image loading into Napari
    # ------------------------------------------------------------------
    def load_current_set(self):
        if self.index >= len(self.experiments):
            print(">> ALL DATASETS PROCESSED.")
            self.finalize_analysis()
            return

        self.current_base_path = self.experiments[self.index]
        self.current_name      = os.path.basename(self.current_base_path)

        files     = sorted(glob.glob(f"{self.current_base_path}_ch*.tif"))

        print(f"\n[{self.index+1}/{len(self.experiments)}] {self.current_name}")

        self.map  = self.detector.assign_roles(files)
        self.mode = "TRIPLE" if self.map['cyan'] else "DUAL"

        self.viewer.layers.select_all()
        self.viewer.layers.remove_selected()

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

        self.shapes_layer = self.viewer.add_shapes(
            name='1_MAIN_CROP',
            edge_color='yellow', face_color='transparent', edge_width=3)
        if self.cfg.do_zoom:
            self.zoom_shapes_layer = self.viewer.add_shapes(
                name='2_ZOOM_ROI',
                edge_color='cyan', face_color='transparent', edge_width=2)
        if self.cfg.do_intensity_profile:
            self.line_shapes_layer = self.viewer.add_shapes(
                name='3_LINE_PROFILE',
                shape_type='line',
                edge_color='white', face_color='transparent', edge_width=2)

        self.add_default_boxes()

    # ------------------------------------------------------------------
    #  Key-binding callbacks
    # ------------------------------------------------------------------
    def skip_and_next(self, viewer):
        print(f"   >> Skipping {self.current_name}...")
        self.index += 1
        self.load_current_set()

    def process_and_next(self, viewer):
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

            lim_c  = self.auto_contrast_limits(raw_crop_c)
            lim_g  = self.auto_contrast_limits(raw_crop_g)
            lim_m  = self.auto_contrast_limits(raw_crop_m)
            lim_bf = self.auto_contrast_limits(raw_crop_bf)

            norm_c  = ensure_shape(self.auto_contrast(raw_crop_c, limits=lim_c))
            norm_g  = ensure_shape(self.auto_contrast(raw_crop_g, limits=lim_g))
            norm_m  = ensure_shape(self.auto_contrast(raw_crop_m, limits=lim_m))
            norm_bf = ensure_shape(self.auto_contrast(raw_crop_bf, limits=lim_bf))

            c_rgb     = self.make_rgb(norm_c, self.cfg.MIX_CYAN)
            g_rgb     = self.make_rgb(norm_g, self.cfg.MIX_GREEN)
            m_rgb     = self.make_rgb(norm_m, self.cfg.MIX_MAGENTA)
            merge_rgb = np.clip(c_rgb + g_rgb + m_rgb, 0, 1)
            bf_rgb    = np.clip(np.stack((norm_bf,)*3, axis=-1), 0, 1).astype(np.float32)

            zoom_rgb            = None
            zoom_coords_in_main = None

            if (self.cfg.do_zoom and self.zoom_shapes_layer is not None
                    and len(self.zoom_shapes_layer.data) > 0):

                zyc, zxc = self.get_box_center(self.zoom_shapes_layer.data[-1])
                zs       = self.cfg.zoom_size    
                zh       = zs // 2
                zy1, zy2 = zyc - zh, zyc + zh
                zx1, zx2 = zxc - zh, zxc + zh
                target_px = self.cfg.zoom_upscaled_size  

                def zoom_crop_upscale_contrast(img_data, limits):
                    crop = self.safe_crop(img_data, zy1, zy2, zx1, zx2)
                    if crop is None or crop.shape != (zs, zs):
                        return None
                    upscaled = self.upscale_channel(crop, target_px, target_px)
                    return self.auto_contrast(upscaled, limits=limits)

                zn_c = zoom_crop_upscale_contrast(img_c, lim_c)
                zn_g = zoom_crop_upscale_contrast(img_g, lim_g)
                zn_m = zoom_crop_upscale_contrast(img_m, lim_m)

                z_rgb = np.zeros((target_px, target_px, 3))
                if zn_c is not None: z_rgb += self.make_rgb(zn_c, self.cfg.MIX_CYAN)
                if zn_g is not None: z_rgb += self.make_rgb(zn_g, self.cfg.MIX_GREEN)
                if zn_m is not None: z_rgb += self.make_rgb(zn_m, self.cfg.MIX_MAGENTA)

                zoom_rgb            = np.clip(z_rgb, 0, 1)
                zoom_coords_in_main = (zy1 - y1, zx1 - x1)

                print(f"   >> ZOOM: {zs}×{zs} → {target_px}×{target_px} "
                      f"({self.cfg.zoom_mag}× bicubic upscale)")

            # --- INTENSITY LINE PROFILE -----------------------------------
            intensity_data   = None
            line_coords_crop = None
            _raw_profiles    = None

            if (self.cfg.do_intensity_profile
                    and self.line_shapes_layer is not None
                    and len(self.line_shapes_layer.data) > 0):

                # Sanitize: remove any accidental shapes from wrong layers
                self._sanitize_shapes_layers()

                line_pts = self.line_shapes_layer.data[-1]  # [[y0,x0],[y1,x1]]
                p0 = (int(line_pts[0][0]) - y1, int(line_pts[0][1]) - x1)
                p1 = (int(line_pts[1][0]) - y1, int(line_pts[1][1]) - x1)
                line_coords_crop = (p0, p1)

                def _profile(arr_2d):
                    """Extract intensity along line, normalise 0-1."""
                    prof = profile_line(arr_2d, p0, p1, order=1, mode='constant')
                    pmin, pmax = prof.min(), prof.max()
                    if pmax > pmin:
                        return (prof - pmin) / (pmax - pmin)
                    return np.zeros_like(prof)

                prof_dict = {'Distance_px': np.arange(len(_profile(norm_g)))}
                if self.mode == "TRIPLE":
                    prof_dict[self.cfg.lbl_c] = _profile(norm_c)
                prof_dict[self.cfg.lbl_g] = _profile(norm_g)
                prof_dict[self.cfg.lbl_m] = _profile(norm_m)

                intensity_data = pd.DataFrame(prof_dict)

                # Store raw arrays for relabeling later
                _raw_profiles = {'distance': prof_dict['Distance_px']}
                if self.mode == "TRIPLE":
                    _raw_profiles['cyan'] = prof_dict[self.cfg.lbl_c]
                _raw_profiles['green'] = prof_dict[self.cfg.lbl_g]
                _raw_profiles['mag']   = prof_dict[self.cfg.lbl_m]

                # Save CSV
                save_dir = os.path.join(self.cfg.output_dir, self.current_name)
                os.makedirs(save_dir, exist_ok=True)
                intensity_data.to_csv(
                    os.path.join(save_dir, "Intensity_Profile.csv"), index=False)
                print(f"   >> INTENSITY PROFILE saved ({len(intensity_data)} pts)")

            self.save_individual_and_local_panel(
                self.current_name,
                c_rgb, g_rgb, m_rgb, merge_rgb, bf_rgb,
                zoom_rgb, zoom_coords_in_main, mask_g, mask_m,
                intensity_data=intensity_data,
                line_coords_crop=line_coords_crop)

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
                'intensity_data':  intensity_data,
                'line_coords':     line_coords_crop,
                'raw_profiles':    _raw_profiles if intensity_data is not None else None,
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
        if img is None:
            return None
        h, w = img.shape[:2]
        y1, x1 = max(0, y1), max(0, x1)
        y2, x2 = min(h, y2), min(w, x2)
        if y1 >= y2 or x1 >= x2:
            return None
        return img[y1:y2, x1:x2]

    def upscale_channel(self, crop, target_h, target_w):
        if crop is None:
            return None
        return resize(
            crop.astype(np.float64), (target_h, target_w),
            order=3, anti_aliasing=False, preserve_range=True)

    # ------------------------------------------------------------------
    #  Figure layout helpers
    # ------------------------------------------------------------------
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
        axes  =[]
        for r in range(n_rows):
            row_top = r * (row_h + self.cfg.spacing_inch)
            bottom  = 1.0 - (row_top + row_h) / total_h
            h_frac  = row_h / total_h

            row_axes =[]
            x = 0.0
            for ci in range(n_cols):
                left   = x / total_w
                w_frac = col_w[ci] / total_w
                ax     = fig.add_axes([left, bottom, w_frac, h_frac])
                row_axes.append(ax)
                x += col_w[ci] + self.cfg.spacing_inch
            axes.append(row_axes)
        return axes

    # ------------------------------------------------------------------
    #  Intensity plot renderer (shared by individual panels & montage)
    # ------------------------------------------------------------------
    def _render_intensity_axes(self, ax, intensity_data, mode='DUAL',
                                linewidth=0.8, tick_scale=1.0):
        """Draw channel intensity traces on *ax*. Transparent background,
        only left & bottom spines, no x-tick labels, plot fills full panel."""
        ax.set_facecolor('none')  # transparent background
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

        # Remove top & right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_edgecolor('black')
        ax.spines['left'].set_linewidth(0.5 * tick_scale)
        ax.spines['bottom'].set_edgecolor('black')
        ax.spines['bottom'].set_linewidth(0.5 * tick_scale)

        ts = self.cfg.font_size * 0.7 * tick_scale
        # Y-axis: keep tick labels.  X-axis: ticks but no labels.
        ax.tick_params(axis='y', colors='black', labelsize=ts,
                       length=2 * tick_scale, pad=1)
        ax.tick_params(axis='x', colors='black', labelsize=ts,
                       length=2 * tick_scale, pad=1,
                       labelbottom=False)  # hide x-axis numbers
        ax.yaxis.set_label_position('left')
        ax.yaxis.tick_left()

        # Adjust position: small left margin for y-tick labels,
        # then let the plot area fill the rest of the panel.
        pos = ax.get_position()
        left_margin = pos.width * 0.10   # space for y-axis ticks
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

    # ------------------------------------------------------------------
    #  Output saving
    # ------------------------------------------------------------------
    def save_individual_and_local_panel(self, name, c, g, m, merge, bf,
                                        zoom, zoom_coords, mask_g, mask_m,
                                        intensity_data=None,
                                        line_coords_crop=None):
        save_path = os.path.join(self.cfg.output_dir, name)
        os.makedirs(save_path, exist_ok=True)

        def to_uint8(arr):
            return (arr * 255).astype(np.uint8)

        if self.mode == "TRIPLE":
            tifffile.imwrite(os.path.join(save_path, "1_Cyan.tif"),    to_uint8(c))
        tifffile.imwrite(os.path.join(save_path, "2_Green.tif"),       to_uint8(g))
        tifffile.imwrite(os.path.join(save_path, "3_Magenta.tif"),     to_uint8(m))
        tifffile.imwrite(os.path.join(save_path, "4_Merge.tif"),       to_uint8(merge))
        if zoom is not None:
            tifffile.imwrite(os.path.join(save_path, "5_Zoom.tif"),    to_uint8(zoom))
        if self.cfg.do_bf:
            tifffile.imwrite(os.path.join(save_path, "6_BF.tif"),      to_uint8(bf))

        if mask_g is not None and mask_m is not None:
            qc_path = os.path.join(save_path, "QC_Masks")
            os.makedirs(qc_path, exist_ok=True)
            tifffile.imwrite(
                os.path.join(qc_path, "Mask_Green_Otsu.tif"), to_uint8(mask_g))
            tifffile.imwrite(
                os.path.join(qc_path, "Mask_Mag_Otsu.tif"),   to_uint8(mask_m))

        plot_list =[]
        if self.mode == "TRIPLE":
            plot_list.append((c,     self.cfg.lbl_c))
        plot_list.append((g,         self.cfg.lbl_g))
        plot_list.append((m,         self.cfg.lbl_m))
        plot_list.append((merge,     self.cfg.lbl_merge))
        if zoom is not None:
            plot_list.append((zoom,  self.cfg.lbl_zoom))
        if self.cfg.do_bf:
            plot_list.append((bf,    self.cfg.lbl_bf))

        # Append intensity profile as the rightmost panel
        has_profile = (intensity_data is not None and not intensity_data.empty)
        if has_profile:
            plot_list.append((None, '__intensity__'))  # sentinel

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

            # --- Intensity profile panel (rightmost) -----------------------
            if label == '__intensity__':
                self._render_intensity_axes(ax, intensity_data, mode=self.mode,
                                            linewidth=0.8, tick_scale=1.0)
                continue

            # --- Normal image panel ----------------------------------------
            ax.imshow(img_data, aspect='auto', interpolation='none')
            ax.axis('off')
            ax.text(0.05, 0.95, label, transform=ax.transAxes, color='white',
                    fontproperties=self.cfg.font_prop, ha='left', va='top')
            if i == 0:
                self.add_scale_bar_patch(ax)

            if label == self.cfg.lbl_merge:
                # Draw zoom ROI rectangle
                if zoom_coords is not None:
                    ry, rx = zoom_coords
                    ax.add_patch(patches.Rectangle(
                        (rx, ry), self.cfg.zoom_size, self.cfg.zoom_size,
                        linewidth=1, edgecolor='white',
                        facecolor='none', linestyle='--'))
                # Draw intensity line
                if line_coords_crop is not None:
                    (ly0, lx0), (ly1, lx1) = line_coords_crop
                    ax.plot([lx0, lx1], [ly0, ly1],
                            color='white', linewidth=0.8, linestyle='--')

        plt.savefig(os.path.join(save_path, "Panel_View.pdf"), dpi=self.cfg.dpi)
        plt.savefig(os.path.join(save_path, "Panel_View.png"), dpi=self.cfg.dpi)
        plt.close(fig)

    # ------------------------------------------------------------------
    #  Finalisation
    # ------------------------------------------------------------------
    def finalize_analysis(self):
        if self.gallery:
            self.build_global_montage()

        if self.cfg.do_quant and self.stats_log:
            df       = pd.DataFrame(self.stats_log)
            csv_path = os.path.join(self.cfg.output_dir,
                                    f"{self.cfg.exp_name}_QUANTIFICATION.csv")
            df.to_csv(csv_path, index=False)
            print(f">> STATS SAVED: {csv_path}")

        # Save session record for future replay
        save_session_record(self.cfg, self.gallery, self.stats_log)

        # Build date-level mega-montage
        gap_px = round(self.cfg.spacing_pt / 72.0 * self.cfg.dpi)
        build_date_summary(self.cfg.date_dir, self.cfg.date_folder, gap_px=gap_px)

    # ------------------------------------------------------------------
    #  Regeneration helpers (used by replay)
    # ------------------------------------------------------------------
    def _regenerate_all_panels(self):
        """Re-draw and overwrite every individual Panel_View using current labels."""
        for item in self.gallery:
            name = item['name']
            self.mode = item['mode']
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
                intensity_data=item.get('intensity_data'),
                line_coords_crop=item.get('line_coords'),
            )
        print(f"   >> {len(self.gallery)} panel view(s) regenerated.")

    def _regenerate_intensity_csvs(self):
        """Re-save intensity CSVs so column headers reflect new labels."""
        for item in self.gallery:
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
        if any(item.get('raw_profiles') for item in self.gallery):
            print("   >> Intensity CSVs updated with new labels.")

    # ------------------------------------------------------------------
    #  Summary Montage
    # ------------------------------------------------------------------
    def build_global_montage(self):
        print("\n>>> GENERATING SUMMARY MONTAGES (PDF + PNG)...")

        n_rows   = len(self.gallery)
        has_cyan = any(x['mode'] == "TRIPLE"    for x in self.gallery)
        has_zoom = any(x['zoom'] is not None     for x in self.gallery)

        has_intensity = any(x.get('intensity_data') is not None for x in self.gallery)

        col_keys   =[]
        col_labels =[]
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
                ax = axes[r][ci]

                # --- Intensity profile column --------------------------------
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

                # --- Normal image column -------------------------------------
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
                    # Draw intensity line on merge panel in montage
                    if item.get('line_coords') is not None:
                        (ly0, lx0), (ly1, lx1) = item['line_coords']
                        ax.plot([lx0, lx1], [ly0, ly1],
                                color='white', linewidth=0.5, linestyle='--')

        base_name = os.path.join(self.cfg.output_dir,
                                 f"{self.cfg.exp_name}_SUMMARY_MONTAGE")
        plt.savefig(f"{base_name}.pdf", dpi=self.cfg.dpi)
        plt.savefig(f"{base_name}.png", dpi=self.cfg.dpi)
        print(f">> MONTAGE SAVED: {base_name}.pdf")
        plt.close(fig)


# ================================================================
#  SESSION RECORD  (pickle save / load)
# ================================================================
def save_session_record(cfg, gallery, stats_log):
    """Pickle the session state so it can be replayed later."""
    record = {
        'cfg': cfg,
        'gallery': gallery,
        'stats_log': stats_log,
    }
    pkl_path = os.path.join(cfg.output_dir,
                            f"{cfg.exp_name}_session_record.pkl")
    with open(pkl_path, 'wb') as f:
        pickle.dump(record, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f">> SESSION RECORD SAVED: {pkl_path}")


def load_session_record(pkl_path):
    """Load a previously saved session record."""
    with open(pkl_path, 'rb') as f:
        record = pickle.load(f)
    return record['cfg'], record['gallery'], record['stats_log']


# ================================================================
#  DATE-LEVEL MEGA-MONTAGE
# ================================================================
def build_date_summary(date_dir, date_folder, gap_px=20):
    """Collect all *_SUMMARY_MONTAGE.png under *date_dir* and arrange
    them in a two-column layout.

    Each sub-montage is scaled to a uniform column width while preserving
    its original aspect ratio.  Images are placed into whichever column
    is currently shorter (masonry packing) to keep columns balanced.
    """
    montage_paths = sorted(glob.glob(
        os.path.join(date_dir, "**", "*_SUMMARY_MONTAGE.png"), recursive=True))
    if len(montage_paths) < 1:
        return

    images = [Image.open(p) for p in montage_paths]
    gap = max(gap_px, 1)

    # Scale every image to the same column width, preserving aspect ratio
    max_w = max(img.width for img in images)
    col_w = max_w  # each column gets this width
    scaled = []
    for img in images:
        if img.width != col_w:
            new_h = round(img.height * col_w / img.width)
            img = img.resize((col_w, new_h), Image.LANCZOS)
        scaled.append(img)

    # Masonry packing: place each image in the shorter column
    col_heights = [0, 0]  # running height for each column
    placements = []       # (col_index, y_offset) per image
    for img in scaled:
        ci = 0 if col_heights[0] <= col_heights[1] else 1
        placements.append((ci, col_heights[ci]))
        col_heights[ci] += img.height + gap

    total_w = col_w * 2 + gap
    total_h = max(col_heights)  # strip trailing gap
    if total_h > gap:
        total_h -= gap

    mega = Image.new('RGBA', (total_w, total_h), (0, 0, 0, 0))
    for img, (ci, y) in zip(scaled, placements):
        x = ci * (col_w + gap)
        mega.paste(img, (x, y))

    base = os.path.join(date_dir, f"{date_folder}_MEGA_SUMMARY")
    mega.save(f"{base}.png", dpi=(300, 300))
    mega.save(f"{base}.pdf", dpi=(300, 300))
    print(f">> MEGA-MONTAGE SAVED: {base}.png / .pdf")


# ================================================================
#  REPLAY  (load record → tweak parameters → regenerate)
# ================================================================
def replay_session():
    """Interactive replay: load a pickle, let user tweak labels/layout, re-render."""
    pkl_path = input("Path to session_record.pkl: ").strip()
    if not os.path.isfile(pkl_path):
        print(f"File not found: {pkl_path}")
        return

    cfg, gallery, stats_log = load_session_record(pkl_path)

    print(f"\nLoaded session: {cfg.exp_name}")
    print(f"  Images: {len(gallery)}")
    print(f"  Output: {cfg.output_dir}")

    # --- Let user tweak any display parameter ---
    print("\n--- ADJUST PARAMETERS (Enter to keep current) ---")
    cfg.lbl_c     = get_user_input(f"Cyan Label",        cfg.lbl_c)
    cfg.lbl_g     = get_user_input(f"Green Label",       cfg.lbl_g)
    cfg.lbl_m     = get_user_input(f"Magenta Label",     cfg.lbl_m)
    cfg.lbl_merge = get_user_input(f"Merge Label",       cfg.lbl_merge)
    if cfg.do_zoom:
        cfg.lbl_zoom = get_user_input(f"Enlargement Label", cfg.lbl_zoom)
    if cfg.do_bf:
        cfg.lbl_bf = get_user_input(f"Brightfield Label", cfg.lbl_bf)

    new_font = get_user_input(f"Font Size (pt)", str(cfg.font_size))
    cfg.font_size = int(new_font)
    cfg.font_prop = fm.FontProperties(family='Arial', size=cfg.font_size)

    new_pw = get_user_input(f"Panel Width (mm)", str(cfg.panel_w_mm))
    cfg.panel_w_mm = float(new_pw)
    cfg.panel_w_inch = cfg.panel_w_mm / 25.4
    if cfg.crop_w > 0:
        cfg.scale = cfg.panel_w_inch / cfg.crop_w
        cfg.panel_h_inch = cfg.crop_h * cfg.scale

    new_dpi = get_user_input(f"DPI", str(cfg.dpi))
    cfg.dpi = int(new_dpi)

    # --- Create a temporary ColocManager-like object to call regeneration ---
    mgr = ColocManager.__new__(ColocManager)
    mgr.cfg = cfg
    mgr.gallery = gallery
    mgr.stats_log = stats_log
    mgr.mode = gallery[0]['mode'] if gallery else 'DUAL'

    # Regenerate intensity CSVs with updated labels
    mgr._regenerate_intensity_csvs()

    # Regenerate all panels
    print("\n>>> Regenerating all panel views...")
    mgr._regenerate_all_panels()

    # Regenerate summary montage
    if gallery:
        print(">>> Regenerating summary montage...")
        mgr.build_global_montage()

    # Save updated session record
    save_session_record(cfg, gallery, stats_log)

    # Rebuild date-level mega-montage
    if hasattr(cfg, 'date_dir') and hasattr(cfg, 'date_folder'):
        gap_px = round(cfg.spacing_pt / 72.0 * cfg.dpi)
        build_date_summary(cfg.date_dir, cfg.date_folder, gap_px=gap_px)

    print("\n>>> Replay complete.")


# ================================================================
#  ENTRY POINT
# ================================================================
if __name__ == "__main__":
    choice = get_user_input(
        "\nNew session [n] / Replay from record [r]", "n").lower()
    if choice == 'r':
        replay_session()
    else:
        ColocManager()