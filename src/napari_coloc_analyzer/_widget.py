"""
_widget.py
==========
Napari dock widget that replaces the terminal wizard and key bindings
from the standalone ``coloc_analyzer.py``.
"""

import os
import glob
import json
import traceback

import napari
import numpy as np
import pandas as pd
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QLabel, QLineEdit, QSpinBox, QDoubleSpinBox,
    QCheckBox, QPushButton, QTextEdit, QFileDialog,
    QScrollArea, QSizePolicy, QProgressBar,
)
from qtpy.QtCore import Qt
from skimage.measure import profile_line

from napari_coloc_analyzer._core import (
    BASE_DIR,
    SessionConfig,
    XMLMetadataDetector,
    Quantifier,
    scan_directory,
    load_flat,
    auto_contrast,
    auto_contrast_limits,
    make_rgb,
    get_box_center,
    safe_crop,
    upscale_channel,
    FigureBuilder,
)


# ================================================================
#  STATE CONSTANTS
# ================================================================
_UNCONFIGURED = "UNCONFIGURED"
_VIEWING = "VIEWING"
_FINALIZED = "FINALIZED"


# ================================================================
#  DOCK WIDGET
# ================================================================
class ColocWidget(QWidget):
    """Single dock widget that drives the full colocalization workflow."""

    _SETTINGS_PATH = os.path.join(
        os.path.expanduser("~"), ".napari_coloc_settings.json")

    def __init__(self, viewer: napari.Viewer):
        super().__init__()
        self.viewer = viewer

        # Processing state
        self.cfg = None
        self.detector = XMLMetadataDetector()
        self.fig_builder = None
        self.experiments = []
        self.index = 0
        self.gallery = []
        self.stats_log = []
        self.current_base_path = None
        self.current_name = None
        self.mode = "DUAL"

        # Shape layer references
        self.shapes_layer = None
        self.zoom_shapes_layer = None
        self.line_shapes_layer = None

        self._build_ui()
        self._load_settings()          # (E) Restore last-used settings
        self._bind_keyboard_shortcuts() # (D) Keyboard shortcuts
        self._set_state(_UNCONFIGURED)

    # ----------------------------------------------------------------
    #  UI CONSTRUCTION
    # ----------------------------------------------------------------
    def _build_ui(self):
        # -- Outer scroll wrapper --
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        inner = QWidget()
        layout = QVBoxLayout(inner)

        # 1. SESSION SETUP
        grp_session = QGroupBox("Session Setup")
        sl = QVBoxLayout()

        row_folder = QHBoxLayout()
        row_folder.addWidget(QLabel("Input Folder:"))
        self.input_dir_edit = QLineEdit(os.path.join(BASE_DIR, "Input_Raw"))
        self.input_dir_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        row_folder.addWidget(self.input_dir_edit)
        self.browse_btn = QPushButton("Browse")
        self.browse_btn.clicked.connect(self._on_browse)
        row_folder.addWidget(self.browse_btn)
        sl.addLayout(row_folder)

        row_exp = QHBoxLayout()
        row_exp.addWidget(QLabel("Experiment:"))
        self.exp_name_edit = QLineEdit()
        self.exp_name_edit.setPlaceholderText("Experiment_YYYY-MM-DD")
        row_exp.addWidget(self.exp_name_edit)
        sl.addLayout(row_exp)

        row_crop = QHBoxLayout()
        row_crop.addWidget(QLabel("Crop (px):  W:"))
        self.crop_w_spin = QSpinBox(); self.crop_w_spin.setRange(10, 9999); self.crop_w_spin.setValue(360)
        row_crop.addWidget(self.crop_w_spin)
        row_crop.addWidget(QLabel(" H:"))
        self.crop_h_spin = QSpinBox(); self.crop_h_spin.setRange(10, 9999); self.crop_h_spin.setValue(360)
        row_crop.addWidget(self.crop_h_spin)
        sl.addLayout(row_crop)

        grp_session.setLayout(sl)
        layout.addWidget(grp_session)
        self.grp_session = grp_session

        # 2. FEATURES
        grp_feat = QGroupBox("Features")
        fl = QVBoxLayout()

        self.bf_check = QCheckBox("Brightfield Panel")
        fl.addWidget(self.bf_check)

        row_zoom = QHBoxLayout()
        self.zoom_check = QCheckBox("Zoom / Enlargement")
        row_zoom.addWidget(self.zoom_check)
        row_zoom.addWidget(QLabel("Magnification:"))
        self.zoom_mag_spin = QSpinBox(); self.zoom_mag_spin.setRange(2, 10); self.zoom_mag_spin.setValue(3)
        self.zoom_mag_spin.setEnabled(False)
        row_zoom.addWidget(self.zoom_mag_spin)
        fl.addLayout(row_zoom)
        self.zoom_check.toggled.connect(self.zoom_mag_spin.setEnabled)

        self.profile_check = QCheckBox("Intensity Line Profile")
        fl.addWidget(self.profile_check)

        self.quant_check = QCheckBox("Quantification (Pearson / Manders)")
        self.quant_check.setChecked(True)
        fl.addWidget(self.quant_check)

        grp_feat.setLayout(fl)
        layout.addWidget(grp_feat)
        self.grp_feat = grp_feat

        # 3. CHANNEL LABELS
        grp_labels = QGroupBox("Channel Labels")
        ll = QVBoxLayout()

        self.lbl_c_edit = self._label_row(ll, "Cyan:", "Protein-Cyan")
        self.lbl_g_edit = self._label_row(ll, "Green:", "Protein-GFP")
        self.lbl_m_edit = self._label_row(ll, "Magenta:", "Protein-mCherry")
        self.lbl_merge_edit = self._label_row(ll, "Merge:", "Merged")
        self.lbl_zoom_edit = self._label_row(ll, "Zoom:", "3X enlarged")
        self.lbl_bf_edit = self._label_row(ll, "BF:", "BF")

        grp_labels.setLayout(ll)
        layout.addWidget(grp_labels)
        self.grp_labels = grp_labels

        # 4. EXPORT SETTINGS
        grp_export = QGroupBox("Export Settings (Illustrator)")
        el = QVBoxLayout()

        row_e1 = QHBoxLayout()
        row_e1.addWidget(QLabel("Panel Width (mm):"))
        self.panel_w_spin = QDoubleSpinBox(); self.panel_w_spin.setRange(1, 200); self.panel_w_spin.setValue(20.0); self.panel_w_spin.setDecimals(1)
        row_e1.addWidget(self.panel_w_spin)
        row_e1.addWidget(QLabel("Spacing (pt):"))
        self.spacing_spin = QDoubleSpinBox(); self.spacing_spin.setRange(0, 50); self.spacing_spin.setValue(3.0); self.spacing_spin.setDecimals(1)
        row_e1.addWidget(self.spacing_spin)
        el.addLayout(row_e1)

        row_e2 = QHBoxLayout()
        row_e2.addWidget(QLabel("Font Size (pt):"))
        self.font_spin = QDoubleSpinBox(); self.font_spin.setRange(1, 72); self.font_spin.setValue(5.0); self.font_spin.setDecimals(1)
        row_e2.addWidget(self.font_spin)
        el.addLayout(row_e2)

        # Scale bar row: µm-based (primary) + pixel size display
        row_sb = QHBoxLayout()
        row_sb.addWidget(QLabel("Scale Bar (µm):"))
        self.sb_um_spin = QDoubleSpinBox()
        self.sb_um_spin.setRange(0.1, 1000)
        self.sb_um_spin.setValue(10.0)
        self.sb_um_spin.setDecimals(1)
        row_sb.addWidget(self.sb_um_spin)

        row_sb.addWidget(QLabel("Pixel Size (µm/px):"))
        self.px_size_spin = QDoubleSpinBox()
        self.px_size_spin.setRange(0, 100)
        self.px_size_spin.setValue(0.0)
        self.px_size_spin.setDecimals(4)
        self.px_size_spin.setToolTip(
            "Auto-detected from Leica XML metadata.\n"
            "Set to 0 to use the manual mm fallback.\n"
            "You can override this value manually.")
        row_sb.addWidget(self.px_size_spin)
        el.addLayout(row_sb)

        # Fallback mm row (used when pixel size is unknown)
        row_sb_mm = QHBoxLayout()
        row_sb_mm.addWidget(QLabel("Scale Bar mm (fallback):"))
        self.sb_spin = QDoubleSpinBox()
        self.sb_spin.setRange(0.1, 50)
        self.sb_spin.setValue(2.646)
        self.sb_spin.setDecimals(3)
        self.sb_spin.setToolTip(
            "Used only when pixel size is not available (0).\n"
            "Specifies scale bar as a physical print-panel length in mm.")
        row_sb_mm.addWidget(self.sb_spin)

        self.sb_px_label = QLabel("")
        self.sb_px_label.setToolTip("Computed scale bar width in pixels")
        row_sb_mm.addWidget(self.sb_px_label)
        el.addLayout(row_sb_mm)

        # Connect signals to update the pixel preview label
        self.sb_um_spin.valueChanged.connect(self._update_sb_preview)
        self.px_size_spin.valueChanged.connect(self._update_sb_preview)
        self._update_sb_preview()

        grp_export.setLayout(el)
        layout.addWidget(grp_export)
        self.grp_export = grp_export

        # === LOAD BUTTON ===
        self.load_btn = QPushButton("LOAD EXPERIMENT")
        self.load_btn.setStyleSheet("font-weight: bold; padding: 8px;")
        self.load_btn.clicked.connect(self._on_load_experiment)
        layout.addWidget(self.load_btn)

        # 5. NAVIGATION
        grp_nav = QGroupBox("Navigation")
        nl = QVBoxLayout()

        self.nav_label = QLabel("No experiment loaded")
        nl.addWidget(self.nav_label)

        # (B) Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setFormat("%v / %m")
        self.progress_bar.setValue(0)
        nl.addWidget(self.progress_bar)

        row_btns = QHBoxLayout()
        # (A) Go-back button
        self.back_btn = QPushButton("\u2190 Back")
        self.back_btn.setToolTip("Return to the previous image (undo last crop)  [B]")
        self.back_btn.clicked.connect(self._on_go_back)
        row_btns.addWidget(self.back_btn)

        self.process_btn = QPushButton("Process && Next")
        self.process_btn.setStyleSheet("font-weight: bold; padding: 6px;")
        self.process_btn.setToolTip("Process current crop and advance  [Enter]")
        self.process_btn.clicked.connect(self._on_process_and_next)
        row_btns.addWidget(self.process_btn)

        self.skip_btn = QPushButton("Skip")
        self.skip_btn.setToolTip("Skip this image without processing  [S]")
        self.skip_btn.clicked.connect(self._on_skip)
        row_btns.addWidget(self.skip_btn)
        nl.addLayout(row_btns)

        grp_nav.setLayout(nl)
        layout.addWidget(grp_nav)
        self.grp_nav = grp_nav

        # 6. STATUS LOG
        grp_log = QGroupBox("Status Log")
        ll2 = QVBoxLayout()
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumHeight(160)
        ll2.addWidget(self.log_text)

        # (C) Export log button
        self.export_log_btn = QPushButton("Save Log to File")
        self.export_log_btn.setToolTip(
            "Export the entire status log to a .txt file in the output folder.")
        self.export_log_btn.clicked.connect(self._on_export_log)
        ll2.addWidget(self.export_log_btn)

        grp_log.setLayout(ll2)
        layout.addWidget(grp_log)

        # === FINALIZE BUTTON ===
        self.finalize_btn = QPushButton("FINALIZE")
        self.finalize_btn.setStyleSheet("font-weight: bold; padding: 8px;")
        self.finalize_btn.clicked.connect(self._on_finalize)
        layout.addWidget(self.finalize_btn)

        # 7. RELABEL
        grp_relabel = QGroupBox("Relabel Channels")
        rl = QVBoxLayout()

        self.relbl_c_edit = self._label_row(rl, "Cyan:", "Protein-Cyan")
        self.relbl_g_edit = self._label_row(rl, "Green:", "Protein-GFP")
        self.relbl_m_edit = self._label_row(rl, "Magenta:", "Protein-mCherry")
        self.relbl_merge_edit = self._label_row(rl, "Merge:", "Merged")
        self.relbl_zoom_edit = self._label_row(rl, "Zoom:", "3X enlarged")
        self.relbl_bf_edit = self._label_row(rl, "BF:", "BF")

        self.relabel_btn = QPushButton("APPLY RELABELING")
        self.relabel_btn.setStyleSheet("font-weight: bold; padding: 6px;")
        self.relabel_btn.clicked.connect(self._on_relabel)
        rl.addWidget(self.relabel_btn)

        grp_relabel.setLayout(rl)
        layout.addWidget(grp_relabel)
        self.grp_relabel = grp_relabel

        # === NEW EXPERIMENT (RESET) BUTTON ===
        self.reset_btn = QPushButton("\u21bb  NEW EXPERIMENT")
        self.reset_btn.setStyleSheet(
            "font-weight: bold; padding: 10px; background-color: #335577; color: white;")
        self.reset_btn.setToolTip(
            "Clear all layers and results, then open a new input folder.")
        self.reset_btn.clicked.connect(self._on_reset)
        layout.addWidget(self.reset_btn)

        # Finalise scroll area
        layout.addStretch()
        scroll.setWidget(inner)
        outer.addWidget(scroll)

    def _update_sb_preview(self, _=None):
        """Refresh the scale bar pixel-width preview label."""
        px_size = self.px_size_spin.value()
        if px_size > 0:
            sb_px = self.sb_um_spin.value() / px_size
            self.sb_px_label.setText(f"= {sb_px:.1f} px")
            self.sb_spin.setEnabled(False)
        else:
            self.sb_px_label.setText("(using mm fallback)")
            self.sb_spin.setEnabled(True)

    @staticmethod
    def _label_row(parent_layout, label_text, default):
        row = QHBoxLayout()
        row.addWidget(QLabel(label_text))
        edit = QLineEdit(default)
        row.addWidget(edit)
        parent_layout.addLayout(row)
        return edit

    # ----------------------------------------------------------------
    #  STATE MANAGEMENT
    # ----------------------------------------------------------------
    def _set_state(self, state):
        self._state = state
        is_unconfigured = (state == _UNCONFIGURED)
        is_viewing = (state == _VIEWING)
        is_finalized = (state == _FINALIZED)

        # Setup sections: editable only when unconfigured
        for w in (self.grp_session, self.grp_feat, self.grp_labels, self.grp_export):
            w.setEnabled(is_unconfigured or is_finalized)

        self.load_btn.setEnabled(is_unconfigured or is_finalized)
        self.grp_nav.setEnabled(is_viewing)
        self.back_btn.setEnabled(is_viewing and self.index > 0)
        self.finalize_btn.setEnabled(is_finalized)
        self.grp_relabel.setEnabled(is_finalized)
        self.reset_btn.setEnabled(is_finalized)
        self.export_log_btn.setEnabled(is_viewing or is_finalized)

    # ----------------------------------------------------------------
    #  LOG HELPER
    # ----------------------------------------------------------------
    def _log(self, msg):
        self.log_text.append(msg)
        self.log_text.verticalScrollBar().setValue(
            self.log_text.verticalScrollBar().maximum()
        )

    # ----------------------------------------------------------------
    #  BROWSE
    # ----------------------------------------------------------------
    def _on_browse(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Input Folder")
        if folder:
            self.input_dir_edit.setText(folder)

    # ----------------------------------------------------------------
    #  LOAD EXPERIMENT
    # ----------------------------------------------------------------
    def _on_load_experiment(self):
        input_dir = self.input_dir_edit.text().strip()
        if not os.path.isdir(input_dir):
            self._log(f"ERROR: Folder not found: {input_dir}")
            return

        # Scan for experiments first so we can detect pixel size
        experiments = scan_directory(input_dir)
        if not experiments:
            self._log("ERROR: No TIF files found in the input folder.")
            return

        # Auto-detect pixel size from the first experiment set's XML
        first_files = sorted(glob.glob(f"{experiments[0]}_ch*.tif"))
        pixel_size = self.detector.extract_pixel_size(first_files, log=self._log)
        if pixel_size is not None:
            self.px_size_spin.setValue(pixel_size)
            self._update_sb_preview()

        # Determine which scale bar path to use
        px_val = self.px_size_spin.value()
        pixel_size_um = px_val if px_val > 0 else None

        # Build SessionConfig from widget values
        self.cfg = SessionConfig(
            input_dir=input_dir,
            exp_name=self.exp_name_edit.text().strip() or None,
            crop_w=self.crop_w_spin.value(),
            crop_h=self.crop_h_spin.value(),
            do_bf=self.bf_check.isChecked(),
            do_zoom=self.zoom_check.isChecked(),
            zoom_mag=self.zoom_mag_spin.value(),
            lbl_c=self.lbl_c_edit.text(),
            lbl_g=self.lbl_g_edit.text(),
            lbl_m=self.lbl_m_edit.text(),
            lbl_merge=self.lbl_merge_edit.text(),
            lbl_zoom=self.lbl_zoom_edit.text() or None,
            lbl_bf=self.lbl_bf_edit.text(),
            do_quant=self.quant_check.isChecked(),
            do_intensity_profile=self.profile_check.isChecked(),
            panel_w_mm=self.panel_w_spin.value(),
            spacing_pt=self.spacing_spin.value(),
            font_size=self.font_spin.value(),
            sb_len_mm=self.sb_spin.value(),
            sb_um=self.sb_um_spin.value(),
            pixel_size_um=pixel_size_um,
        )
        self.fig_builder = FigureBuilder(self.cfg)

        self.experiments = experiments

        self.index = 0
        self.gallery = []
        self.stats_log = []

        # (B) Initialise progress bar
        self.progress_bar.setMaximum(len(self.experiments))
        self.progress_bar.setValue(0)

        self._log(f"Found {len(self.experiments)} experiment set(s).")
        self._save_settings()          # (E) Persist settings for next launch
        self._set_state(_VIEWING)
        self._load_current_set()

    # ----------------------------------------------------------------
    #  IMAGE LOADING INTO NAPARI
    # ----------------------------------------------------------------
    def _load_current_set(self):
        if self.index >= len(self.experiments):
            self._log(">> ALL DATASETS PROCESSED.")
            self._on_finalize()
            return

        self.current_base_path = self.experiments[self.index]
        self.current_name = os.path.basename(self.current_base_path)

        files = sorted(glob.glob(f"{self.current_base_path}_ch*.tif"))

        # (B) Update progress bar + back-button state
        self.progress_bar.setValue(self.index)
        self.back_btn.setEnabled(self.index > 0)

        self._log(f"\n[{self.index+1}/{len(self.experiments)}] {self.current_name}")

        self.map = self.detector.assign_roles(files, log=self._log)
        self.mode = "TRIPLE" if self.map['cyan'] else "DUAL"

        # Clear existing layers
        self.viewer.layers.select_all()
        self.viewer.layers.remove_selected()

        if self.map['cyan']:
            self.viewer.add_image(load_flat(self.map['cyan']),
                name='Cyan_Ch', colormap='cyan', blending='additive')
        if self.map['green']:
            self.viewer.add_image(load_flat(self.map['green']),
                name='Green_Ch', colormap='green', blending='additive')
        if self.map['mag']:
            self.viewer.add_image(load_flat(self.map['mag']),
                name='Mag_Ch', colormap='magenta', blending='additive')
        if self.cfg.do_bf and self.map['bf']:
            self.viewer.add_image(load_flat(self.map['bf']),
                name='BF_Ch', colormap='gray',
                blending='translucent', opacity=0.4)

        self.shapes_layer = self.viewer.add_shapes(
            name='1_MAIN_CROP',
            edge_color='yellow', face_color='transparent', edge_width=3)
        if self.cfg.do_zoom:
            self.zoom_shapes_layer = self.viewer.add_shapes(
                name='2_ZOOM_ROI',
                edge_color='cyan', face_color='transparent', edge_width=2)
        else:
            self.zoom_shapes_layer = None
        if self.cfg.do_intensity_profile:
            self.line_shapes_layer = self.viewer.add_shapes(
                name='3_LINE_PROFILE',
                shape_type='line',
                edge_color='white', face_color='transparent', edge_width=2)
        else:
            self.line_shapes_layer = None

        self._add_default_boxes()
        self.nav_label.setText(
            f"Image {self.index+1} / {len(self.experiments)}: {self.current_name}")

    # ----------------------------------------------------------------
    #  SHAPES HELPERS
    # ----------------------------------------------------------------
    def _add_default_boxes(self):
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

    def _sanitize_shapes_layers(self):
        def _filter(layer, keep_npts):
            if layer is None or len(layer.data) == 0:
                return
            cleaned = [s for s in layer.data if len(s) == keep_npts]
            removed = len(layer.data) - len(cleaned)
            if removed > 0:
                self._log(f"   >> WARNING: Removed {removed} stray shape(s) "
                          f"from '{layer.name}' (wrong vertex count).")
                layer.data = cleaned if cleaned else layer.data[:0]

        _filter(self.shapes_layer, 4)
        if self.cfg.do_zoom:
            _filter(self.zoom_shapes_layer, 4)
        if self.cfg.do_intensity_profile:
            _filter(self.line_shapes_layer, 2)

    # ----------------------------------------------------------------
    #  PROCESS & NEXT
    # ----------------------------------------------------------------
    def _on_process_and_next(self):
        if self.shapes_layer is None or len(self.shapes_layer.data) == 0:
            self._log("   >> No crop box found!")
            return

        yc, xc = get_box_center(self.shapes_layer.data[-1])
        rh, rw = self.cfg.crop_h // 2, self.cfg.crop_w // 2
        y1, y2 = yc - rh, yc + rh
        x1, x2 = xc - rw, xc + rw

        try:
            def get_layer_data(name):
                return self.viewer.layers[name].data if name in self.viewer.layers else None

            img_c  = get_layer_data('Cyan_Ch')
            img_g  = get_layer_data('Green_Ch')
            img_m  = get_layer_data('Mag_Ch')
            img_bf = get_layer_data('BF_Ch')

            raw_crop_c  = safe_crop(img_c,  y1, y2, x1, x2)
            raw_crop_g  = safe_crop(img_g,  y1, y2, x1, x2)
            raw_crop_m  = safe_crop(img_m,  y1, y2, x1, x2)
            raw_crop_bf = safe_crop(img_bf, y1, y2, x1, x2)

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
                    self.current_name, stats_g, stats_m, log=self._log)
                self.stats_log.append(stats)
                self._log(f"   >> STATS: R={stats['Pearson_R']}, "
                          f"M1={stats['Manders_M1 (G_in_M)']}")
            else:
                self._log("   >> Visuals only")

            lim_c  = auto_contrast_limits(raw_crop_c)
            lim_g  = auto_contrast_limits(raw_crop_g)
            lim_m  = auto_contrast_limits(raw_crop_m)
            lim_bf = auto_contrast_limits(raw_crop_bf)

            norm_c  = ensure_shape(auto_contrast(raw_crop_c, limits=lim_c))
            norm_g  = ensure_shape(auto_contrast(raw_crop_g, limits=lim_g))
            norm_m  = ensure_shape(auto_contrast(raw_crop_m, limits=lim_m))
            norm_bf = ensure_shape(auto_contrast(raw_crop_bf, limits=lim_bf))

            c_rgb     = make_rgb(norm_c, self.cfg.MIX_CYAN)
            g_rgb     = make_rgb(norm_g, self.cfg.MIX_GREEN)
            m_rgb     = make_rgb(norm_m, self.cfg.MIX_MAGENTA)
            merge_rgb = np.clip(c_rgb + g_rgb + m_rgb, 0, 1)
            bf_rgb    = np.clip(np.stack((norm_bf,)*3, axis=-1), 0, 1).astype(np.float32)

            zoom_rgb = None
            zoom_coords_in_main = None

            if (self.cfg.do_zoom and self.zoom_shapes_layer is not None
                    and len(self.zoom_shapes_layer.data) > 0):

                zyc, zxc = get_box_center(self.zoom_shapes_layer.data[-1])
                zs = self.cfg.zoom_size
                zh = zs // 2
                zy1, zy2 = zyc - zh, zyc + zh
                zx1, zx2 = zxc - zh, zxc + zh
                target_px = self.cfg.zoom_upscaled_size

                def zoom_crop_upscale_contrast(img_data, limits):
                    crop = safe_crop(img_data, zy1, zy2, zx1, zx2)
                    if crop is None or crop.shape != (zs, zs):
                        return None
                    upscaled = upscale_channel(crop, target_px, target_px)
                    return auto_contrast(upscaled, limits=limits)

                zn_c = zoom_crop_upscale_contrast(img_c, lim_c)
                zn_g = zoom_crop_upscale_contrast(img_g, lim_g)
                zn_m = zoom_crop_upscale_contrast(img_m, lim_m)

                z_rgb = np.zeros((target_px, target_px, 3))
                if zn_c is not None: z_rgb += make_rgb(zn_c, self.cfg.MIX_CYAN)
                if zn_g is not None: z_rgb += make_rgb(zn_g, self.cfg.MIX_GREEN)
                if zn_m is not None: z_rgb += make_rgb(zn_m, self.cfg.MIX_MAGENTA)

                zoom_rgb = np.clip(z_rgb, 0, 1)
                zoom_coords_in_main = (zy1 - y1, zx1 - x1)

                self._log(f"   >> ZOOM: {zs}x{zs} -> {target_px}x{target_px} "
                          f"({self.cfg.zoom_mag}x bicubic upscale)")

            # --- INTENSITY LINE PROFILE ---
            intensity_data = None
            line_coords_crop = None
            _raw_profiles = None

            if (self.cfg.do_intensity_profile
                    and self.line_shapes_layer is not None
                    and len(self.line_shapes_layer.data) > 0):

                self._sanitize_shapes_layers()

                line_pts = self.line_shapes_layer.data[-1]
                p0 = (int(line_pts[0][0]) - y1, int(line_pts[0][1]) - x1)
                p1 = (int(line_pts[1][0]) - y1, int(line_pts[1][1]) - x1)
                line_coords_crop = (p0, p1)

                def _profile(arr_2d):
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

                _raw_profiles = {'distance': prof_dict['Distance_px']}
                if self.mode == "TRIPLE":
                    _raw_profiles['cyan'] = prof_dict[self.cfg.lbl_c]
                _raw_profiles['green'] = prof_dict[self.cfg.lbl_g]
                _raw_profiles['mag'] = prof_dict[self.cfg.lbl_m]

                save_dir = os.path.join(self.cfg.output_dir, self.current_name)
                os.makedirs(save_dir, exist_ok=True)
                intensity_data.to_csv(
                    os.path.join(save_dir, "Intensity_Profile.csv"), index=False)
                self._log(f"   >> INTENSITY PROFILE saved ({len(intensity_data)} pts)")

            self.fig_builder.save_individual_and_local_panel(
                self.current_name,
                c_rgb, g_rgb, m_rgb, merge_rgb, bf_rgb,
                zoom_rgb, zoom_coords_in_main, mask_g, mask_m,
                mode=self.mode,
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
                'intensity_data': intensity_data,
                'line_coords':    line_coords_crop,
                'raw_profiles':   _raw_profiles if intensity_data is not None else None,
            })

            self.index += 1
            self._load_current_set()

        except Exception as e:
            self._log(f"ERROR: {e}")
            self._log(traceback.format_exc())

    # ----------------------------------------------------------------
    #  SKIP
    # ----------------------------------------------------------------
    def _on_skip(self):
        self._log(f"   >> Skipping {self.current_name}...")
        self.index += 1
        self._load_current_set()

    # ----------------------------------------------------------------
    #  (A) GO BACK / UNDO
    # ----------------------------------------------------------------
    def _on_go_back(self):
        """Return to the previous image, discarding the last processed entry."""
        if self.index <= 0:
            self._log("   >> Already at the first image.")
            return

        # If the last gallery entry corresponds to the previous index, remove it
        prev_base = self.experiments[self.index - 1]
        prev_name = os.path.basename(prev_base)
        if self.gallery and self.gallery[-1]['name'] == prev_name:
            self.gallery.pop()
            if self.stats_log and self.stats_log[-1].get('Image') == prev_name:
                self.stats_log.pop()
            self._log(f"   << Undid processing for {prev_name}")
        else:
            self._log(f"   << Going back to {prev_name} (was skipped)")

        self.index -= 1
        self._load_current_set()

    # ----------------------------------------------------------------
    #  (C) EXPORT LOG
    # ----------------------------------------------------------------
    def _on_export_log(self):
        """Save the status log to a text file."""
        # Determine save location
        if self.cfg and self.cfg.output_dir:
            default_path = os.path.join(self.cfg.output_dir, "session_log.txt")
        else:
            default_path = os.path.join(BASE_DIR, "session_log.txt")

        path, _ = QFileDialog.getSaveFileName(
            self, "Save Log File", default_path, "Text files (*.txt)")
        if not path:
            return
        try:
            with open(path, 'w') as f:
                f.write(self.log_text.toPlainText())
            self._log(f"   >> Log exported to {path}")
        except Exception as e:
            self._log(f"ERROR saving log: {e}")

    # ----------------------------------------------------------------
    #  (D) KEYBOARD SHORTCUTS
    # ----------------------------------------------------------------
    def _bind_keyboard_shortcuts(self):
        """Register viewer key bindings for fast navigation."""
        @self.viewer.bind_key('Enter', overwrite=True)
        def _key_process(viewer):
            if self._state == _VIEWING:
                self._on_process_and_next()

        @self.viewer.bind_key('s', overwrite=True)
        def _key_skip(viewer):
            if self._state == _VIEWING:
                self._on_skip()

        @self.viewer.bind_key('b', overwrite=True)
        def _key_back(viewer):
            if self._state == _VIEWING and self.index > 0:
                self._on_go_back()

    # ----------------------------------------------------------------
    #  (E) REMEMBER LAST SETTINGS
    # ----------------------------------------------------------------
    def _save_settings(self):
        """Persist current widget values to a JSON file."""
        data = {
            'input_dir':    self.input_dir_edit.text(),
            'exp_name':     self.exp_name_edit.text(),
            'crop_w':       self.crop_w_spin.value(),
            'crop_h':       self.crop_h_spin.value(),
            'do_bf':        self.bf_check.isChecked(),
            'do_zoom':      self.zoom_check.isChecked(),
            'zoom_mag':     self.zoom_mag_spin.value(),
            'lbl_c':        self.lbl_c_edit.text(),
            'lbl_g':        self.lbl_g_edit.text(),
            'lbl_m':        self.lbl_m_edit.text(),
            'lbl_merge':    self.lbl_merge_edit.text(),
            'lbl_zoom':     self.lbl_zoom_edit.text(),
            'lbl_bf':       self.lbl_bf_edit.text(),
            'do_quant':     self.quant_check.isChecked(),
            'do_profile':   self.profile_check.isChecked(),
            'panel_w_mm':   self.panel_w_spin.value(),
            'spacing_pt':   self.spacing_spin.value(),
            'font_size':    self.font_spin.value(),
            'sb_um':        self.sb_um_spin.value(),
            'sb_len_mm':    self.sb_spin.value(),
        }
        try:
            with open(self._SETTINGS_PATH, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception:
            pass   # Non-critical; silently ignore write errors

    def _load_settings(self):
        """Restore widget values from the saved JSON file, if it exists."""
        if not os.path.isfile(self._SETTINGS_PATH):
            return
        try:
            with open(self._SETTINGS_PATH) as f:
                d = json.load(f)
        except Exception:
            return

        # Apply values defensively (missing keys are simply skipped)
        if d.get('input_dir'):    self.input_dir_edit.setText(d['input_dir'])
        if d.get('exp_name'):     self.exp_name_edit.setText(d['exp_name'])
        if 'crop_w' in d:        self.crop_w_spin.setValue(d['crop_w'])
        if 'crop_h' in d:        self.crop_h_spin.setValue(d['crop_h'])
        if 'do_bf' in d:         self.bf_check.setChecked(d['do_bf'])
        if 'do_zoom' in d:       self.zoom_check.setChecked(d['do_zoom'])
        if 'zoom_mag' in d:      self.zoom_mag_spin.setValue(d['zoom_mag'])
        if d.get('lbl_c'):       self.lbl_c_edit.setText(d['lbl_c'])
        if d.get('lbl_g'):       self.lbl_g_edit.setText(d['lbl_g'])
        if d.get('lbl_m'):       self.lbl_m_edit.setText(d['lbl_m'])
        if d.get('lbl_merge'):   self.lbl_merge_edit.setText(d['lbl_merge'])
        if d.get('lbl_zoom'):    self.lbl_zoom_edit.setText(d['lbl_zoom'])
        if d.get('lbl_bf'):      self.lbl_bf_edit.setText(d['lbl_bf'])
        if 'do_quant' in d:      self.quant_check.setChecked(d['do_quant'])
        if 'do_profile' in d:    self.profile_check.setChecked(d['do_profile'])
        if 'panel_w_mm' in d:    self.panel_w_spin.setValue(d['panel_w_mm'])
        if 'spacing_pt' in d:    self.spacing_spin.setValue(d['spacing_pt'])
        if 'font_size' in d:     self.font_spin.setValue(d['font_size'])
        if 'sb_um' in d:         self.sb_um_spin.setValue(d['sb_um'])
        if 'sb_len_mm' in d:     self.sb_spin.setValue(d['sb_len_mm'])

    # ----------------------------------------------------------------
    #  FINALIZE
    # ----------------------------------------------------------------
    def _on_finalize(self):
        if self.fig_builder and self.gallery:
            self.fig_builder.finalize_analysis(
                self.gallery, self.stats_log, log=self._log)
        self._log("All datasets processed. Montage and CSV saved.")
        self._set_state(_FINALIZED)

        # Pre-fill relabel fields with current labels
        self.relbl_c_edit.setText(self.cfg.lbl_c)
        self.relbl_g_edit.setText(self.cfg.lbl_g)
        self.relbl_m_edit.setText(self.cfg.lbl_m)
        self.relbl_merge_edit.setText(self.cfg.lbl_merge)
        self.relbl_zoom_edit.setText(self.cfg.lbl_zoom)
        self.relbl_bf_edit.setText(self.cfg.lbl_bf)

    # ----------------------------------------------------------------
    #  RELABEL
    # ----------------------------------------------------------------
    def _on_relabel(self):
        self.cfg.lbl_c = self.relbl_c_edit.text()
        self.cfg.lbl_g = self.relbl_g_edit.text()
        self.cfg.lbl_m = self.relbl_m_edit.text()
        self.cfg.lbl_merge = self.relbl_merge_edit.text()
        self.cfg.lbl_zoom = self.relbl_zoom_edit.text()
        self.cfg.lbl_bf = self.relbl_bf_edit.text()

        self._log(">>> Regenerating with new labels...")
        self.fig_builder.regenerate_intensity_csvs(self.gallery, log=self._log)
        self.fig_builder.regenerate_all_panels(self.gallery, log=self._log)
        if self.gallery:
            self.fig_builder.build_global_montage(self.gallery, log=self._log)
        self._log(">>> Relabeling complete.")

    # ----------------------------------------------------------------
    #  RESET (NEW EXPERIMENT)
    # ----------------------------------------------------------------
    def _on_reset(self):
        """Clear all layers, results, and state so the user can load a
        new input folder without restarting Napari."""
        # 1. Remove all viewer layers
        self.viewer.layers.select_all()
        self.viewer.layers.remove_selected()

        # 2. Reset processing state
        self.cfg = None
        self.fig_builder = None
        self.experiments = []
        self.index = 0
        self.gallery = []
        self.stats_log = []
        self.current_base_path = None
        self.current_name = None
        self.mode = "DUAL"
        self.shapes_layer = None
        self.zoom_shapes_layer = None
        self.line_shapes_layer = None

        # 3. Reset pixel size (will be re-detected from new folder)
        self.px_size_spin.setValue(0.0)
        self._update_sb_preview()

        # 4. Clear UI feedback
        self.nav_label.setText("No experiment loaded")
        self.log_text.clear()
        self._log("Session reset. Ready for a new experiment folder.")

        # 5. Return to unconfigured state
        self._set_state(_UNCONFIGURED)
