# Confocal Colocalization Tool

An interactive, semi-automated tool for **fluorescence colocalization analysis** of multi-channel confocal microscopy images.  
Built on [Napari](https://napari.org/), it lets you visually select regions of interest, crop them, apply pseudocolour overlays, draw intensity line profiles, and compute **Pearson's R** and **Manders' coefficients** — all with publication-ready PDF/PNG panel output fully compatible with **Adobe Illustrator**.

---

## Features

| Feature | Description |
|---|---|
| **Interactive ROI selection** | Drag a yellow box in Napari to select your region of interest |
| **Zoom panel** | Optional cyan box to define an enlarged inset with bicubic upscaling |
| **Dual / Triple channel** | Automatically detects 2- or 3-fluorescence-channel datasets |
| **Brightfield support** | Optional BF overlay panel (true-colour RGB) |
| **Fluorescent intensity line profile** | Draw a line in Napari; normalised intensity CSV and integrated plot panel are generated automatically |
| **Colocalization quantification** | Pearson R and Manders M1/M2 with **Maximum Entropy** (Kapur) thresholding |
| **Publication-ready output** | PDF panels with `interpolation='none'`, Type 42 vector fonts, configurable print size (mm) — optimised for Adobe Illustrator |
| **Summary montage** | Auto-assembled multi-row montage PDF + PNG including intensity profile column |
| **Post-processing channel relabelling** | Rename channel labels after Napari closes; all panels, CSVs, and montage are regenerated instantly |
| **Shape layer safety** | Automatic vertex-count sanitiser prevents accidental cross-layer drawing |
| **CSV export** | Colocalization statistics and intensity profiles saved per image |

---

## Requirements

- Python ≥ 3.9
- See [`requirements.txt`](requirements.txt) for all dependencies

### Recommended setup

> It is strongly recommended to use a dedicated **conda** or **venv** environment.

```bash
# conda (recommended)
conda create -n coloc python=3.11
conda activate coloc
pip install -r requirements.txt
```

```bash
# venv alternative
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

**VS Code users:** Run *Python: Select Interpreter* → choose the environment, then *Developer: Reload Window* if import warnings persist.

---

## Input File Format

The tool expects **single-channel TIFF files** named with a `_ch<N>.tif` suffix:

```
MyExperiment_s001_ch00.tif   ← Cyan / first fluorophore
MyExperiment_s001_ch01.tif   ← Green / second fluorophore
MyExperiment_s001_ch02.tif   ← Magenta / third fluorophore
MyExperiment_s001_ch03.tif   ← Brightfield (optional)
```

- Files belonging to the **same field of view** share a common prefix before `_ch`.
- Sub-directories are scanned **recursively**.
- Compatible with the default export format of **Zeiss ZEN** and **Leica LAS X**.

---

## Folder Structure

By default the tool reads from and writes to `~/Desktop/Confocal_Workflow/` (override at runtime):

```
Confocal_Workflow/
├── Input_Raw/                          ← Place raw .tif files here
└── Output_Coloc/
    └── <ExperimentName>/
        ├── <ImageBasename>/
        │   ├── 1_Cyan.tif
        │   ├── 2_Green.tif
        │   ├── 3_Magenta.tif
        │   ├── 4_Merge.tif
        │   ├── 5_Zoom.tif              (if zoom enabled)
        │   ├── 6_BF.tif                (if BF enabled)
        │   ├── Intensity_Profile.csv    (if intensity profile enabled)
        │   ├── Panel_View.pdf
        │   ├── Panel_View.png
        │   └── QC_Masks/
        │       ├── Mask_Green_MaxEnt.tif
        │       └── Mask_Mag_MaxEnt.tif
        ├── <ExperimentName>_SUMMARY_MONTAGE.pdf
        ├── <ExperimentName>_SUMMARY_MONTAGE.png
        └── <ExperimentName>_QUANTIFICATION.csv
```

---

## Usage

```bash
python coloc_analyzer.py
```

### Setup Wizard

The script opens an **interactive wizard** in the terminal:

1. **Input folder** — path to your raw TIFF directory
2. **Experiment name** — used for naming output files
3. **Crop geometry** — width × height in pixels
4. **Extra features** — brightfield panel, zoom/inset panel, intensity line profile
5. **Channel labels** — text labels printed on each panel
6. **Illustrator export settings** — physical panel size (mm), font size, scale bar length
7. **Quantification** — enable/disable Pearson + Manders calculation

### Napari Controls

| Key | Action |
|---|---|
| `c` | Crop current ROI, save outputs, advance to next image |
| `s` | Skip current image |

### Interactive Layers

| Layer (colour) | Purpose |
|---|---|
| **1_MAIN_CROP** (yellow) | Main crop region — drag to position |
| **2_ZOOM_ROI** (cyan) | Zoom inset region — must be inside the yellow box |
| **3_LINE_PROFILE** (white) | Intensity line — draw a 2-point line for the fluorescent intensity profile |

> The shape-layer sanitiser automatically removes accidental shapes drawn in the wrong layer so you don't have to worry about cross-layer mistakes.

### Post-processing Relabelling

After closing Napari, the tool offers a **channel relabelling loop**:

```
--- Channel relabelling (enter blank line to finish) ---
Current labels:
  Cyan  = DAPI
  Green = GFP
  Mag   = mCherry

New Cyan label [DAPI]:
```

Type a new label or press Enter to keep the current one. All panels (PDF/PNG), intensity CSVs, and the summary montage are regenerated automatically with the updated labels — **no need to re-run the entire pipeline**.

---

## Output Description

### Per-image

| File | Description |
|---|---|
| `1_Cyan.tif` … `4_Merge.tif` | Pseudocolour **uint8 RGB** TIFF crops |
| `5_Zoom.tif` | Bicubic-upscaled inset (if zoom enabled) |
| `6_BF.tif` | Brightfield RGB (if BF enabled) |
| `Intensity_Profile.csv` | Normalised fluorescence intensity along the drawn line (if enabled) |
| `Panel_View.pdf` | Vector PDF panel — Illustrator-ready (`interpolation='none'`, Type 42 fonts) |
| `Panel_View.png` | Raster PNG at calculated DPI |
| `QC_Masks/` | Maximum Entropy threshold binary masks for QC |

### Global

| File | Description |
|---|---|
| `*_SUMMARY_MONTAGE.pdf/png` | Multi-row montage of all processed images (includes intensity profile column) |
| `*_QUANTIFICATION.csv` | Pearson R, Manders M1/M2, threshold values per image |

---

## Colocalization Metrics

| Metric | Description |
|---|---|
| **Pearson R** | Linear correlation of pixel intensities between the green and magenta channels |
| **Manders M1** | Fraction of green signal overlapping with magenta mask |
| **Manders M2** | Fraction of magenta signal overlapping with green mask |

Thresholds for Manders coefficients are determined by **Maximum Entropy thresholding** (Kapur et al., 1985) applied independently to each channel crop. This method is more robust than Otsu for sparse fluorescent puncta common in confocal images.

---

## Intensity Line Profile

When enabled, an additional **white line layer** (`3_LINE_PROFILE`) appears in Napari. Draw a two-point line across the structure of interest:

- A **normalised intensity CSV** is saved per image with columns for each channel and distance in pixels.
- An **intensity plot panel** is appended as the rightmost column in both the per-image `Panel_View` and the summary montage.
- A **dashed white line** is drawn on the merge panel to indicate the profiled region.

The plot uses minimal styling: left/bottom axes only (L-shaped), no x-axis numbers, and a transparent background to integrate cleanly with your figure layout.

---

## Adobe Illustrator Compatibility

All PDF output is optimised for Adobe Illustrator:

- **Type 42 fonts** (`pdf.fonttype = 42`, `ps.fonttype = 42`) — text remains editable in Illustrator
- **`interpolation='none'`** on all image panels — raw pixel data is embedded, avoiding resampling artefacts
- **True-colour RGB** arrays for every channel (including brightfield) — consistent colour model across all panels

---

## Configuration Tips

- **Scale bar**: `sb_len_mm` should match the physical length you want to display. Calculate from your microscope's pixel size × crop width.
- **DPI**: Automatically calculated from crop width (px) ÷ panel width (mm / 25.4). Typical values: 300–600 DPI for print.
- **Panel width**: 20 mm is standard for a single journal figure column. Adjust to your journal's requirements.

---

## License

This project is released under the [MIT License](LICENSE).

---

## Citation

If you use this tool in published research, please cite this repository:

```
StevenChung13. Confocal Colocalization Tool. GitHub. https://github.com/StevenChung13/confocal-colocalization-tool
```
