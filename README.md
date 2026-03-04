# Confocal Colocalization Tool

An interactive, semi-automated tool for **fluorescence colocalization analysis** of multi-channel confocal microscopy images.  
Built on [Napari](https://napari.org/), it lets you visually select regions of interest, crops them, applies pseudocolour overlays, and computes **Pearson's R** and **Manders' coefficients** — all with publication-ready PDF/PNG panel output.

---

## Features

| Feature | Description |
|---|---|
| Interactive ROI selection | Drag a yellow box in Napari to select your region of interest |
| Zoom panel | Optional cyan box to define an enlarged inset |
| Dual/Triple channel | Automatically detects cyan, green, and magenta channels |
| Brightfield support | Optional BF overlay panel |
| Colocalization quantification | Pearson R and Manders M1/M2 with Otsu thresholding |
| Publication-ready output | PDF panels at configurable print size (mm), vector fonts (Type 42) |
| Summary montage | Auto-assembled multi-row montage PDF + PNG |
| CSV export | All colocalization statistics saved per experiment |

---

## Requirements

- Python >= 3.9
- See [`requirements.txt`](requirements.txt) for all dependencies

### Recommended local setup (VS Code)

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Then in VS Code:

1. Run **Python: Select Interpreter**
2. Choose `.venv/bin/python`
3. Run **Developer: Reload Window** if `Import "napari" could not be resolved` persists

Install dependencies with:

```bash
pip install -r requirements.txt
```

> **Tip:** It is strongly recommended to use a dedicated conda or venv environment.

```bash
conda create -n coloc python=3.11
conda activate coloc
pip install -r requirements.txt
```

---

## Input File Format

The tool expects **single-channel TIFF files** named with a `_ch<N>.tif` suffix:

```
MyExperiment_s001_ch00.tif   <- Cyan / first fluorophore
MyExperiment_s001_ch01.tif   <- Green / second fluorophore
MyExperiment_s001_ch02.tif   <- Magenta / third fluorophore
MyExperiment_s001_ch03.tif   <- Brightfield (optional)
```

- Files belonging to the **same field of view** share a common prefix before `_ch`.
- Sub-directories are scanned **recursively**.
- This naming convention is the default export format of **Zeiss ZEN** and **Leica LAS X** software.

---

## Folder Structure

By default, the tool expects (and creates) the following structure under `~/Desktop/Confocal_Workflow/`:

```
Confocal_Workflow/
├── Input_Raw/          <- Place your raw .tif files here
└── Output_Coloc/
    └── <ExperimentName>/
        ├── <ImageBasename>/
        │   ├── 1_Cyan.tif
        │   ├── 2_Green.tif
        │   ├── 3_Magenta.tif
        │   ├── 4_Merge.tif
        │   ├── 5_Zoom.tif          (if zoom enabled)
        │   ├── 6_BF.tif            (if BF enabled)
        │   ├── Panel_View.pdf
        │   ├── Panel_View.png
        │   └── QC_Masks/
        │       ├── Mask_Green_Otsu.tif
        │       └── Mask_Mag_Otsu.tif
        ├── <ExperimentName>_SUMMARY_MONTAGE.pdf
        ├── <ExperimentName>_SUMMARY_MONTAGE.png
        └── <ExperimentName>_QUANTIFICATION.csv
```

You can override the input/output paths at runtime through the setup wizard.

---

## Usage

```bash
python coloc_analyzer.py
```

The script opens an **interactive wizard** in the terminal:

1. **Input folder** — path to your raw TIFF directory  
2. **Experiment name** — used for naming output files  
3. **Crop geometry** — width x height in pixels  
4. **Extra features** — brightfield panel, zoom/inset panel  
5. **Channel labels** — text labels printed on each panel  
6. **Illustrator export settings** — physical panel size (mm), font size, scale bar length  
7. **Quantification** — enable/disable Pearson + Manders calculation

### Napari Controls

| Key | Action |
|---|---|
| `c` | Crop current ROI, save outputs, advance to next image |
| `s` | Skip current image |

**Zoom workflow:**  
If zoom is enabled, two boxes appear:  
- **Yellow** box -> main crop region  
- **Cyan** box -> zoom inset region (must be fully inside the yellow box)

---

## Output Description

### Per-image
| File | Description |
|---|---|
| `1_Cyan.tif` ... `4_Merge.tif` | Pseudocolour uint8 TIFF crops |
| `5_Zoom.tif` | Bicubic-upscaled inset (if zoom enabled) |
| `Panel_View.pdf` | Vector PDF panel for Adobe Illustrator |
| `Panel_View.png` | Raster PNG at configured DPI |
| `QC_Masks/` | Otsu threshold binary masks for QC |

### Global
| File | Description |
|---|---|
| `*_SUMMARY_MONTAGE.pdf/png` | Multi-row summary of all processed images |
| `*_QUANTIFICATION.csv` | Pearson R, Manders M1/M2, Otsu thresholds |

---

## Colocalization Metrics

| Metric | Description |
|---|---|
| **Pearson R** | Linear correlation of pixel intensities between channels |
| **Manders M1** | Fraction of green signal overlapping with magenta mask |
| **Manders M2** | Fraction of magenta signal overlapping with green mask |

Thresholds for Manders coefficients are determined by **Otsu's method** applied independently to each channel crop.

---

## Configuration Tips

- **Scale bar**: `sb_len_mm` should match the physical length you want to display. Calculate it from your microscope's pixel size x crop width.  
- **DPI**: Automatically calculated from crop width (px) divided by panel width (mm / 25.4). Typical values: 300–600 DPI for print.  
- **Panel width**: 20 mm is standard for a journal figure column. Adjust to your journal's requirements.

---

## License

This project is released under the [MIT License](LICENSE).

---

## Citation

If you use this tool in published research, please cite this repository:

```
StevenChung13. Confocal Colocalization Tool. GitHub. https://github.com/StevenChung13/confocal-colocalization-tool
```
