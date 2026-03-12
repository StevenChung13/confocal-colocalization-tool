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
    QScrollArea, QSizePolicy, QProgressBar, QComboBox,
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
    unified_contrast_limits,
    make_rgb,
    get_box_center,
    safe_crop,
    upscale_channel,
    denoise_image,
    deconvolve_image,
    generate_gaussian_psf,
    load_psf,
    FigureBuilder,
    save_session_record,
    load_session_record,
    build_date_summary,
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

        row_date = QHBoxLayout()
        row_date.addWidget(QLabel("Date Folder:"))
        self.date_folder_edit = QLineEdit()
        import datetime as _dt
        self.date_folder_edit.setPlaceholderText(_dt.datetime.now().strftime("%Y%m%d"))
        row_date.addWidget(self.date_folder_edit)
        sl.addLayout(row_date)

        row_exp = QHBoxLayout()
        row_exp.addWidget(QLabel("Experiment:"))
        self.exp_name_edit = QLineEdit()
        self.exp_name_edit.setPlaceholderText("Experiment_YYYYMMDD")
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

        # --- Denoise controls ---
        self.denoise_check = QCheckBox("Denoise")
        fl.addWidget(self.denoise_check)

        row_dn = QHBoxLayout()
        row_dn.addWidget(QLabel("Method:"))
        self.denoise_method_combo = QComboBox()
        self.denoise_method_combo.addItems([
            "Wavelet", "Non-Local Means", "Total Variation",
            "Median", "Gaussian", "Bilateral"])
        _dn_tips = {
            0: "Wavelet — Best all-round. Preserves edges and fine"
               " structure.  Use for most confocal images.",
            1: "Non-Local Means — Averages similar patches."
               "  Best when noise is moderate and textures repeat.",
            2: "Total Variation — Reduces noise while keeping"
               " sharp edges.  Good for piecewise-constant signals.",
            3: "Median — Fast rank filter.  Removes salt-and-pepper"
               " / shot noise without blurring edges.",
            4: "Gaussian — Simple smoothing.  Use for mild,"
               " uniform noise when speed matters.",
            5: "Bilateral — Edge-aware smoothing.  Preserves"
               " boundaries; good for mixed-texture samples.",
        }
        for idx, tip in _dn_tips.items():
            self.denoise_method_combo.setItemData(idx, tip, Qt.ToolTipRole)
        self.denoise_method_combo.setEnabled(False)
        row_dn.addWidget(self.denoise_method_combo)
        row_dn.addWidget(QLabel("Strength:"))
        self.denoise_strength_spin = QDoubleSpinBox()
        self.denoise_strength_spin.setRange(0.1, 10.0)
        self.denoise_strength_spin.setValue(1.0)
        self.denoise_strength_spin.setSingleStep(0.1)
        self.denoise_strength_spin.setDecimals(1)
        self.denoise_strength_spin.setEnabled(False)
        row_dn.addWidget(self.denoise_strength_spin)
        fl.addLayout(row_dn)
        self.denoise_check.toggled.connect(self.denoise_method_combo.setEnabled)
        self.denoise_check.toggled.connect(self.denoise_strength_spin.setEnabled)

        # --- Deconvolution controls ---
        self.deconv_check = QCheckBox("Deconvolution")
        fl.addWidget(self.deconv_check)

        row_dc = QHBoxLayout()
        row_dc.addWidget(QLabel("Method:"))
        self.deconv_method_combo = QComboBox()
        self.deconv_method_combo.addItems([
            "Richardson-Lucy", "Wiener", "Unsupervised Wiener"])
        _dc_tips = {
            0: "Richardson-Lucy — Iterative ML method.  Best for"
               " Poisson (photon) noise typical of confocal/"
               "fluorescence data.  More iterations = sharper.",
            1: "Wiener — Single-step frequency-domain filter."
               "  Very fast; good when noise level is roughly"
               " known and speed is a priority.",
            2: "Unsupervised Wiener — Auto-estimates noise power."
               "  Use when the noise level is unknown; no manual"
               " tuning needed.",
        }
        for idx, tip in _dc_tips.items():
            self.deconv_method_combo.setItemData(idx, tip, Qt.ToolTipRole)
        self.deconv_method_combo.setEnabled(False)
        row_dc.addWidget(self.deconv_method_combo)
        row_dc.addWidget(QLabel("Iterations:"))
        self.deconv_iter_spin = QSpinBox()
        self.deconv_iter_spin.setRange(1, 200)
        self.deconv_iter_spin.setValue(15)
        self.deconv_iter_spin.setEnabled(False)
        row_dc.addWidget(self.deconv_iter_spin)
        fl.addLayout(row_dc)
        self.deconv_check.toggled.connect(self.deconv_method_combo.setEnabled)
        self.deconv_check.toggled.connect(self.deconv_iter_spin.setEnabled)

        row_psf = QHBoxLayout()
        row_psf.addWidget(QLabel("PSF:"))
        self.psf_source_combo = QComboBox()
        self.psf_source_combo.addItems(["Synthetic", "Load from file"])
        self.psf_source_combo.setEnabled(False)
        row_psf.addWidget(self.psf_source_combo)
        fl.addLayout(row_psf)
        self.deconv_check.toggled.connect(self.psf_source_combo.setEnabled)

        # Synthetic PSF parameters
        row_psf_params = QHBoxLayout()
        row_psf_params.addWidget(QLabel("NA:"))
        self.na_spin = QDoubleSpinBox()
        self.na_spin.setRange(0.1, 1.7)
        self.na_spin.setValue(1.4)
        self.na_spin.setDecimals(2)
        self.na_spin.setSingleStep(0.05)
        self.na_spin.setEnabled(False)
        row_psf_params.addWidget(self.na_spin)
        row_psf_params.addWidget(QLabel("λ (nm):"))
        self.wavelength_spin = QDoubleSpinBox()
        self.wavelength_spin.setRange(300, 800)
        self.wavelength_spin.setValue(520.0)
        self.wavelength_spin.setDecimals(0)
        self.wavelength_spin.setEnabled(False)
        row_psf_params.addWidget(self.wavelength_spin)
        fl.addLayout(row_psf_params)
        self.deconv_check.toggled.connect(self.na_spin.setEnabled)
        self.deconv_check.toggled.connect(self.wavelength_spin.setEnabled)

        # PSF file path (shown only when "Load from file" is selected)
        row_psf_file = QHBoxLayout()
        row_psf_file.addWidget(QLabel("PSF File:"))
        self.psf_path_edit = QLineEdit()
        self.psf_path_edit.setPlaceholderText("Path to measured PSF TIFF")
        self.psf_path_edit.setEnabled(False)
        row_psf_file.addWidget(self.psf_path_edit)
        self.psf_browse_btn = QPushButton("Browse")
        self.psf_browse_btn.setEnabled(False)
        self.psf_browse_btn.clicked.connect(self._on_browse_psf)
        row_psf_file.addWidget(self.psf_browse_btn)
        fl.addLayout(row_psf_file)

        def _toggle_psf_source(idx):
            is_file = (idx == 1)
            is_on = self.deconv_check.isChecked()
            self.psf_path_edit.setEnabled(is_file and is_on)
            self.psf_browse_btn.setEnabled(is_file and is_on)
            self.na_spin.setEnabled(not is_file and is_on)
            self.wavelength_spin.setEnabled(not is_file and is_on)
        self.psf_source_combo.currentIndexChanged.connect(_toggle_psf_source)

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

        # === LOAD BUTTONS ===
        row_load = QHBoxLayout()
        self.load_btn = QPushButton("LOAD EXPERIMENT")
        self.load_btn.setStyleSheet("font-weight: bold; padding: 8px;")
        self.load_btn.clicked.connect(self._on_load_experiment)
        row_load.addWidget(self.load_btn)

        self.load_record_btn = QPushButton("LOAD RECORD")
        self.load_record_btn.setStyleSheet("font-weight: bold; padding: 8px;")
        self.load_record_btn.setToolTip(
            "Load a previous session record (.pkl) to tweak settings and re-generate.")
        self.load_record_btn.clicked.connect(self._on_load_record)
        row_load.addWidget(self.load_record_btn)
        layout.addLayout(row_load)

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

        # 7. SESSION RECORD (save only — load button is above with LOAD EXPERIMENT)
        grp_record = QGroupBox("Session Record")
        rl = QVBoxLayout()

        self.save_record_btn = QPushButton("SAVE SESSION RECORD")
        self.save_record_btn.setStyleSheet("font-weight: bold; padding: 6px;")
        self.save_record_btn.setToolTip(
            "Save cfg + gallery + stats as a pickle for future replay.")
        self.save_record_btn.clicked.connect(self._on_save_record)
        rl.addWidget(self.save_record_btn)

        grp_record.setLayout(rl)
        layout.addWidget(grp_record)
        self.grp_record = grp_record

        # === NEW EXPERIMENT (RESET) BUTTON ===
        self.reset_btn = QPushButton("\u21bb  NEW EXPERIMENT")
        self.reset_btn.setStyleSheet(
            "font-weight: bold; padding: 10px; background-color: #335577; color: white;")
        self.reset_btn.setToolTip(
            "Clear all layers and results, then open a new input folder.")
        self.reset_btn.clicked.connect(self._on_reset)
        layout.addWidget(self.reset_btn)

        # --- Collapsible sections ---
        self._collapse_btns = {}
        for grp, collapsed in [
            (grp_session, False),
            (grp_feat,    False),
            (grp_labels,  True),
            (grp_export,  True),
            (grp_nav,     False),
            (grp_log,     False),
            (grp_record,  True),
        ]:
            self._make_collapsible(grp, layout, collapsed)

        # Finalise scroll area
        layout.addStretch()
        scroll.setWidget(inner)
        outer.addWidget(scroll)

    # ----------------------------------------------------------------
    #  COLLAPSIBLE HELPER
    # ----------------------------------------------------------------
    def _make_collapsible(self, grp, parent_layout, collapsed=False):
        """Insert a flat toggle button before *grp* that hides/shows it."""
        title = grp.title()
        arrow = "\u25B6" if collapsed else "\u25BC"  # ▶ / ▼
        btn = QPushButton(f"{arrow}  {title}")
        btn.setFlat(True)
        btn.setStyleSheet(
            "text-align: left; font-weight: bold; padding: 3px 6px;"
            "color: palette(text);")
        btn.setCursor(Qt.PointingHandCursor)
        btn.setFocusPolicy(Qt.NoFocus)

        idx = parent_layout.indexOf(grp)
        parent_layout.insertWidget(idx, btn)

        grp.setTitle("")            # title now lives on the button
        grp.setFlat(True)           # remove the box frame for a cleaner look
        grp.setVisible(not collapsed)

        def _toggle():
            vis = not grp.isVisible()
            grp.setVisible(vis)
            btn.setText(f"{'\u25BC' if vis else '\u25B6'}  {title}")

        btn.clicked.connect(_toggle)
        self._collapse_btns[grp] = btn

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
            # Keep collapse toggle clickable regardless of state
            if w in self._collapse_btns:
                self._collapse_btns[w].setEnabled(True)

        self.load_btn.setEnabled(is_unconfigured or is_finalized)
        self.load_record_btn.setEnabled(is_unconfigured or is_finalized)
        self.grp_nav.setEnabled(is_viewing)
        self.back_btn.setEnabled(is_viewing and self.index > 0)
        self.finalize_btn.setEnabled(is_finalized)
        self.grp_record.setEnabled(is_finalized)
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

    def _on_browse_psf(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select PSF TIFF", "", "TIFF files (*.tif *.tiff);;All files (*)")
        if path:
            self.psf_path_edit.setText(path)

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

        # Map combo display names to internal keys
        _denoise_methods = {
            "Wavelet": "wavelet", "Non-Local Means": "nlm",
            "Total Variation": "tv", "Median": "median",
            "Gaussian": "gaussian", "Bilateral": "bilateral",
        }
        _deconv_methods = {
            "Richardson-Lucy": "richardson_lucy",
            "Wiener": "wiener",
            "Unsupervised Wiener": "unsupervised_wiener",
        }
        _psf_sources = {"Synthetic": "synthetic", "Load from file": "file"}

        # Build SessionConfig from widget values
        self.cfg = SessionConfig(
            input_dir=input_dir,
            date_folder=self.date_folder_edit.text().strip() or None,
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
            do_denoise=self.denoise_check.isChecked(),
            denoise_method=_denoise_methods.get(
                self.denoise_method_combo.currentText(), "wavelet"),
            denoise_strength=self.denoise_strength_spin.value(),
            do_deconvolve=self.deconv_check.isChecked(),
            deconv_method=_deconv_methods.get(
                self.deconv_method_combo.currentText(), "richardson_lucy"),
            deconv_iterations=self.deconv_iter_spin.value(),
            psf_source=_psf_sources.get(
                self.psf_source_combo.currentText(), "synthetic"),
            psf_path=self.psf_path_edit.text().strip() or None,
            na=self.na_spin.value(),
            wavelength_nm=self.wavelength_spin.value(),
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

            # --- DENOISE (before deconvolution) ---
            if self.cfg.do_denoise:
                meth = self.cfg.denoise_method
                stren = self.cfg.denoise_strength
                self._log(f"   >> DENOISE: {meth} (strength={stren})")
                if raw_crop_c is not None:
                    raw_crop_c = denoise_image(raw_crop_c, meth, strength=stren)
                if raw_crop_g is not None:
                    raw_crop_g = denoise_image(raw_crop_g, meth, strength=stren)
                if raw_crop_m is not None:
                    raw_crop_m = denoise_image(raw_crop_m, meth, strength=stren)

            # --- DECONVOLUTION (after denoise) ---
            psf = None
            if self.cfg.do_deconvolve:
                if self.cfg.psf_source == "file" and self.cfg.psf_path:
                    psf = load_psf(self.cfg.psf_path)
                    self._log(f"   >> PSF loaded from file: {self.cfg.psf_path}")
                else:
                    psf = generate_gaussian_psf(
                        self.cfg.na, self.cfg.wavelength_nm,
                        self.cfg.pixel_size_um)
                    self._log(f"   >> Synthetic PSF (NA={self.cfg.na}, "
                              f"λ={self.cfg.wavelength_nm}nm)")

                dm = self.cfg.deconv_method
                nit = self.cfg.deconv_iterations
                self._log(f"   >> DECONVOLVE: {dm} ({nit} iter)")
                if raw_crop_c is not None:
                    raw_crop_c = deconvolve_image(raw_crop_c, psf, dm, num_iter=nit)
                if raw_crop_g is not None:
                    raw_crop_g = deconvolve_image(raw_crop_g, psf, dm, num_iter=nit)
                if raw_crop_m is not None:
                    raw_crop_m = deconvolve_image(raw_crop_m, psf, dm, num_iter=nit)

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

            # --- CONTRAST LIMITS (unified when denoise/deconv active) ---
            if self.cfg.do_denoise or self.cfg.do_deconvolve:
                shared_lim = unified_contrast_limits(
                    raw_crop_c, raw_crop_g, raw_crop_m)
                lim_c = lim_g = lim_m = shared_lim
                self._log(f"   >> UNIFIED CONTRAST: {shared_lim}")
            else:
                lim_c = auto_contrast_limits(raw_crop_c)
                lim_g = auto_contrast_limits(raw_crop_g)
                lim_m = auto_contrast_limits(raw_crop_m)
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
                    if self.cfg.do_denoise:
                        crop = denoise_image(
                            crop, self.cfg.denoise_method,
                            strength=self.cfg.denoise_strength)
                    if self.cfg.do_deconvolve and psf is not None:
                        crop = deconvolve_image(
                            crop, psf, self.cfg.deconv_method,
                            num_iter=self.cfg.deconv_iterations)
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
            'date_folder':  self.date_folder_edit.text(),
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
            # Denoise / Deconvolution
            'do_denoise':       self.denoise_check.isChecked(),
            'denoise_method':   self.denoise_method_combo.currentText(),
            'denoise_strength': self.denoise_strength_spin.value(),
            'do_deconvolve':    self.deconv_check.isChecked(),
            'deconv_method':    self.deconv_method_combo.currentText(),
            'deconv_iterations': self.deconv_iter_spin.value(),
            'psf_source':       self.psf_source_combo.currentText(),
            'psf_path':         self.psf_path_edit.text(),
            'na':               self.na_spin.value(),
            'wavelength_nm':    self.wavelength_spin.value(),
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
        if d.get('date_folder'):  self.date_folder_edit.setText(d['date_folder'])
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
        # Denoise / Deconvolution
        if 'do_denoise' in d:       self.denoise_check.setChecked(d['do_denoise'])
        if d.get('denoise_method'):
            idx = self.denoise_method_combo.findText(d['denoise_method'])
            if idx >= 0: self.denoise_method_combo.setCurrentIndex(idx)
        if 'denoise_strength' in d: self.denoise_strength_spin.setValue(d['denoise_strength'])
        if 'do_deconvolve' in d:    self.deconv_check.setChecked(d['do_deconvolve'])
        if d.get('deconv_method'):
            idx = self.deconv_method_combo.findText(d['deconv_method'])
            if idx >= 0: self.deconv_method_combo.setCurrentIndex(idx)
        if 'deconv_iterations' in d: self.deconv_iter_spin.setValue(d['deconv_iterations'])
        if d.get('psf_source'):
            idx = self.psf_source_combo.findText(d['psf_source'])
            if idx >= 0: self.psf_source_combo.setCurrentIndex(idx)
        if d.get('psf_path'):       self.psf_path_edit.setText(d['psf_path'])
        if 'na' in d:               self.na_spin.setValue(d['na'])
        if 'wavelength_nm' in d:    self.wavelength_spin.setValue(d['wavelength_nm'])

    # ----------------------------------------------------------------
    #  FINALIZE
    # ----------------------------------------------------------------
    def _on_finalize(self):
        if not self.cfg or not self.gallery:
            self._log("Nothing to finalize — no processed data.")
            return

        # --- Sync current UI values back into cfg -----------------------
        self.cfg.lbl_c     = self.lbl_c_edit.text()
        self.cfg.lbl_g     = self.lbl_g_edit.text()
        self.cfg.lbl_m     = self.lbl_m_edit.text()
        self.cfg.lbl_merge = self.lbl_merge_edit.text()
        self.cfg.lbl_zoom  = self.lbl_zoom_edit.text()
        self.cfg.lbl_bf    = self.lbl_bf_edit.text()

        new_pw = self.panel_w_spin.value()
        if new_pw != self.cfg.panel_w_mm:
            self.cfg.panel_w_mm = new_pw
            self.cfg.panel_w_inch = new_pw / 25.4
            if self.cfg.crop_w > 0:
                self.cfg.scale = self.cfg.panel_w_inch / self.cfg.crop_w
                self.cfg.panel_h_inch = self.cfg.crop_h * self.cfg.scale
                self.cfg.dpi = self.cfg.crop_w / self.cfg.panel_w_inch

        new_fs = self.font_spin.value()
        if new_fs != self.cfg.font_size:
            self.cfg.font_size = new_fs
            self.cfg.font_prop.set_size(new_fs)

        self.cfg.spacing_pt = self.spacing_spin.value()
        self.cfg.spacing_inch = self.cfg.spacing_pt / 72.0

        # Rebuild FigureBuilder with updated cfg
        self.fig_builder = FigureBuilder(self.cfg)

        # --- Regenerate individual panels & CSVs ------------------------
        self._log(">>> Regenerating individual panels with current settings...")
        self.fig_builder.regenerate_intensity_csvs(self.gallery, log=self._log)
        self.fig_builder.regenerate_all_panels(self.gallery, log=self._log)

        # --- Montage, CSV, session record, mega-montage -----------------
        self.fig_builder.finalize_analysis(
            self.gallery, self.stats_log, log=self._log)
        self._log("All outputs regenerated. Montage, CSV, and session record saved.")
        self._set_state(_FINALIZED)

    # ----------------------------------------------------------------
    #  SESSION RECORD
    # ----------------------------------------------------------------
    def _on_save_record(self):
        if not self.cfg or not self.gallery:
            self._log("Nothing to save — no processed data.")
            return
        save_session_record(self.cfg, self.gallery, self.stats_log)
        self._log(f"Session record saved to {self.cfg.output_dir}")

    def _on_load_record(self):
        pkl_path, _ = QFileDialog.getOpenFileName(
            self, "Select Session Record", "",
            "Pickle files (*.pkl);;All files (*)")
        if not pkl_path:
            return

        cfg, gallery, stats_log = load_session_record(pkl_path)

        # Populate widget fields from loaded config
        self.cfg = cfg
        self.gallery = gallery
        self.stats_log = stats_log
        self.fig_builder = FigureBuilder(cfg)

        self.date_folder_edit.setText(getattr(cfg, 'date_folder', ''))
        self.exp_name_edit.setText(cfg.exp_name)
        self.lbl_c_edit.setText(cfg.lbl_c)
        self.lbl_g_edit.setText(cfg.lbl_g)
        self.lbl_m_edit.setText(cfg.lbl_m)
        self.lbl_merge_edit.setText(cfg.lbl_merge)
        self.lbl_zoom_edit.setText(cfg.lbl_zoom)
        self.lbl_bf_edit.setText(cfg.lbl_bf)
        self.panel_w_spin.setValue(cfg.panel_w_mm)
        self.font_spin.setValue(cfg.font_size)
        self.spacing_spin.setValue(cfg.spacing_pt)

        self._log(f"Loaded record: {cfg.exp_name} ({len(gallery)} images)")
        self._log("Edit labels / export settings above, then click FINALIZE to regenerate.")
        self._set_state(_FINALIZED)

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
