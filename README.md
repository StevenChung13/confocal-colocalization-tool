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
| **Date folder organisation** | Group experiments by date (`YYYYMMDD`) for cleaner project structure |
| **Session record (pickle)** | Save your entire session (settings + processed data) to a `.pkl` file; reload later to change labels, font, layout — no reprocessing needed |
| **Mega-montage** | Automatically combines all experiment montages from a date folder into a single two-column overview (PNG + PDF) |
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

By default the tool reads from and writes to `~/Desktop/Confocal_Workflow/`:

```
Confocal_Workflow/
├── Input_Raw/                              ← Place raw .tif files here
└── Output_Coloc/
    └── <DateFolder>/                       ← e.g. 20260312
        ├── <ExperimentName>/
        │   ├── <ImageBasename>/
        │   │   ├── 1_Cyan.tif
        │   │   ├── 2_Green.tif
        │   │   ├── 3_Magenta.tif
        │   │   ├── 4_Merge.tif
        │   │   ├── 5_Zoom.tif              (if zoom enabled)
        │   │   ├── 6_BF.tif                (if BF enabled)
        │   │   ├── Intensity_Profile.csv    (if intensity profile enabled)
        │   │   ├── Panel_View.pdf
        │   │   ├── Panel_View.png
        │   │   └── QC_Masks/
        │   │       ├── Mask_Green_MaxEnt.tif
        │   │       └── Mask_Mag_MaxEnt.tif
        │   ├── <Experiment>_SUMMARY_MONTAGE.pdf
        │   ├── <Experiment>_SUMMARY_MONTAGE.png
        │   ├── <Experiment>_QUANTIFICATION.csv
        │   └── <Experiment>_session_record.pkl
        ├── <DateFolder>_MEGA_SUMMARY.png    ← all experiments combined
        └── <DateFolder>_MEGA_SUMMARY.pdf
```

---

## Quick-Start Workflow (Step-by-Step)

This guide assumes you have never used the tool before. Follow each step in order.

### 1. Prepare your images

1. Export your confocal images as **single-channel TIFFs** (one file per channel).
2. Place all TIFFs in `~/Desktop/Confocal_Workflow/Input_Raw/` (or any folder you like).

### 2. Launch the plugin

```bash
# Napari widget (recommended)
python -m napari_coloc_analyzer
```

Napari opens with the **Colocalization Analyzer** panel docked on the right.

### 3. Configure your session

Fill in the settings in the dock widget:

| Setting | What to enter | Example |
|---|---|---|
| **Input Folder** | Click **Browse** → navigate to your TIFF folder | `~/Desktop/Confocal_Workflow/Input_Raw` |
| **Date Folder** | The date of your experiment in `YYYYMMDD` format (defaults to today) | `20260312` |
| **Experiment** | A descriptive name for this batch | `ATG5_KO_Starvation` |
| **Crop (px)** | Width × Height of your crop region | `360 × 360` |

### 4. Enable the features you need

| Checkbox | When to enable |
|---|---|
| **Brightfield Panel** | You have a transmitted-light / DIC channel |
| **Zoom / Enlargement** | You want a magnified inset (set magnification, e.g. 3×) |
| **Intensity Line Profile** | You want a fluorescence intensity plot across a structure |
| **Quantification** | You want Pearson R and Manders M1/M2 (enabled by default) |

### 5. Set channel labels

Enter the names of your fluorescent proteins / dyes (e.g. "ATG5", "LC3B-GFP", "LAMP1-mCherry"). These appear on every panel and montage.

### 6. Adjust export settings

| Setting | Recommendation |
|---|---|
| **Panel Width (mm)** | 20 mm for a single journal column; 40 mm for double |
| **Font Size (pt)** | 5–7 pt for most journals |
| **Scale Bar (µm)** | 10 µm is typical for cell biology |
| **Pixel Size (µm/px)** | Auto-detected from Leica XML metadata; adjust manually if needed |

### 7. Click **LOAD EXPERIMENT**

The plugin scans your folder and loads the first image set into Napari with coloured shape layers:

| Layer (colour) | What to do |
|---|---|
| **Yellow box** (1_MAIN_CROP) | Drag to position over your cell / structure of interest |
| **Cyan box** (2_ZOOM_ROI) | Drag inside the yellow box to define the zoom inset (if enabled) |
| **White line** (3_LINE_PROFILE) | Drag endpoints across the structure for the intensity profile (if enabled) |

### 8. Process each image

For each field of view:

1. **Position the yellow box** over the area you want to crop
2. **Position the cyan box** and **white line** if applicable
3. Click **Process & Next** (or press **Enter**) → the tool crops, generates panels, computes stats, and advances to the next image
4. Click **Skip** (or press **S**) to skip an image
5. Click **← Back** (or press **B**) to return to the previous image

The progress bar shows how many images remain.

### 9. Click **FINALIZE**

When all images are done, click **FINALIZE**. This:

- Generates the **summary montage** (PDF + PNG)
- Saves the **quantification CSV**
- Saves a **session record** (`.pkl` file)
- Builds the **date-level mega-montage** combining all experiments from the same date folder

### 10. Done!

Your outputs are in `Output_Coloc/<DateFolder>/<ExperimentName>/`. Open the PDFs in Illustrator for final figure assembly.

---

## Reloading a Session Record (Change Labels / Settings Without Reprocessing)

This is the key time-saving feature. If you need to change channel labels, font size, panel width, or spacing **after** you have already processed and finalized:

### Step-by-step

1. **Launch the plugin** (or click **↻ NEW EXPERIMENT** if already open)
2. Click **LOAD RECORD** (next to "LOAD EXPERIMENT")
3. In the file dialog, navigate to your experiment folder and select the `*_session_record.pkl` file
4. The plugin populates **all settings** from the saved session — labels, font, spacing, etc.
5. **Edit whatever you want** in the dock widget:
   - Change channel labels (e.g. rename "Protein-GFP" → "LC3B")
   - Adjust font size or panel width
   - Change spacing
6. Click **FINALIZE**
7. **All outputs are regenerated** — individual Panel_View images, intensity CSVs, the summary montage, the session record, and the mega-montage — using your new settings

> **No need to re-select ROIs or reprocess images.** The session record stores all cropped image data, so regeneration is instant.

### Standalone script (terminal)

```bash
python coloc_analyzer.py
# Choose:  New session [n] / Replay from record [r]
# Press 'r', then paste the path to your .pkl file
# Follow the prompts to edit labels, font, DPI, etc.
```

---

## Date Folder & Mega-Montage

### Date folder

Every experiment is now saved under a **date folder** (`YYYYMMDD`) inside `Output_Coloc/`. This keeps your daily work organised:

```
Output_Coloc/
├── 20260312/
│   ├── ATG5_KO_Starvation/
│   ├── ATG5_WT_Control/
│   └── 20260312_MEGA_SUMMARY.png
└── 20260313/
    └── ...
```

### Mega-montage

After finalization, all `*_SUMMARY_MONTAGE.png` files in the date folder are automatically combined into a **two-column mega-montage**:

- Each sub-montage is scaled to the same column width **without distortion** (aspect ratio preserved)
- Images are placed in a balanced two-column masonry layout
- The gap between images matches your configured spacing
- Saved as both PNG (with transparent background) and PDF

This gives you a single overview of all experiments done on that date.

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

### Per-experiment

| File | Description |
|---|---|
| `*_SUMMARY_MONTAGE.pdf/png` | Multi-row montage of all processed images (includes intensity profile column) |
| `*_QUANTIFICATION.csv` | Pearson R, Manders M1/M2, threshold values per image |
| `*_session_record.pkl` | Saved session for future replay (labels, settings, all image data) |

### Per-date

| File | Description |
|---|---|
| `*_MEGA_SUMMARY.png` | Two-column overview of all experiment montages (transparent background) |
| `*_MEGA_SUMMARY.pdf` | Same mega-montage in PDF format |

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

- **Scale bar**: Enter the desired scale bar length in **µm** (e.g., 10). The pixel size is **auto-detected** from Leica LAS X XML metadata (`Voxel` attribute or `Length / NumberOfElements`). If metadata is unavailable, the tool falls back to a manual **mm** value specifying the bar's physical print-panel length.
- **Pixel size override**: The auto-detected pixel size (µm/px) is shown in the widget and can be manually adjusted if needed. Setting it to 0 reverts to the mm-based fallback.
- **DPI**: Automatically calculated from crop width (px) ÷ panel width (mm / 25.4). Typical values: 300–600 DPI for print.
- **Panel width**: 20 mm is standard for a single journal figure column. Adjust to your journal's requirements.

---

## Common Problems & Solutions

### Installation issues

| Problem | Solution |
|---|---|
| `ModuleNotFoundError: No module named 'imagecodecs'` | Run `pip install imagecodecs`. This is required for reading LZW/Deflate TIFFs exported by Leica LAS X and Zeiss ZEN. |
| `ModuleNotFoundError: No module named 'napari'` | Make sure you are in the correct conda/venv environment: `conda activate coloc` (or `source .venv/bin/activate`), then `pip install -r requirements.txt`. |
| Napari does not open / crashes immediately | Try `pip install "napari[all]"` to install all optional Napari dependencies. On macOS, ensure you have the latest PyQt5: `pip install --upgrade PyQt5`. |
| `pip install` fails on Apple Silicon (M1/M2/M3) | Use `conda install` for NumPy and SciPy first: `conda install numpy scipy scikit-image`, then `pip install -r requirements.txt`. |

### Image loading issues

| Problem | Solution |
|---|---|
| "No image sets found" after clicking LOAD EXPERIMENT | Check that your TIFFs follow the `_ch00.tif` / `_ch01.tif` naming convention. The tool matches files by the common prefix before `_ch`. |
| Only 2 channels detected (expected 3) | Ensure the third channel file (`_ch02.tif`) exists and shares the same prefix. If you have only 2 fluorescent channels, the tool automatically runs in **Dual** mode — this is normal. |
| Images look black in Napari | Your raw TIFF may have very low signal. The tool applies auto-contrast (0.02–99.98 percentile) during processing — the final output will look correct even if the Napari preview appears dark. |
| "Could not read TIFF" / `tifffile.TiffFileError` | Your TIFF may use an unsupported compression. Re-export from your microscopy software as **uncompressed** or **LZW-compressed** TIFF. |

### Processing issues

| Problem | Solution |
|---|---|
| Yellow box is missing or won't move | Make sure the **1_MAIN_CROP** layer is selected (click it in the layer list). Use the **pointer tool** (arrow icon, top left of Napari), not the drawing tool. |
| Accidentally drew extra shapes | The shape layer sanitiser automatically removes extra vertices. If you still see unexpected shapes, select the layer → press **A** (select all) → **Delete**, then re-draw. |
| Zoom box appears outside the yellow crop box | The zoom inset must be **inside** the yellow box. If placed outside, the zoom panel will show the wrong region. Drag the cyan box to be fully contained within the yellow box. |
| Intensity line profile is flat | The line must cross a region with fluorescent signal. If drawn over background, the profile will be flat. Redraw across a bright structure. |
| Output panels look pixelated | Increase the **Panel Width (mm)** value — this increases DPI. At 20 mm with a 360 px crop, DPI ≈ 457. For higher resolution, use a larger crop size. |

### Session record / replay issues

| Problem | Solution |
|---|---|
| "Nothing to finalize" when clicking FINALIZE | You need to either process images first (LOAD EXPERIMENT → Process) or load a session record (LOAD RECORD). |
| Labels didn't change after FINALIZE with loaded record | This was a known bug — **now fixed**. FINALIZE reads all current UI values before regenerating. Make sure you are running the latest version. |
| `.pkl` file is very large | Session records store all cropped image arrays. A 20-image experiment with 360×360 crops uses ~50–100 MB. This is normal and ensures instant replay without reprocessing. |
| `LOAD RECORD` button is greyed out | The button should be available at startup (UNCONFIGURED state) and after finalization. If greyed out, click **↻ NEW EXPERIMENT** to reset to startup state. |

### Output issues

| Problem | Solution |
|---|---|
| PDF text is not editable in Illustrator | Ensure you haven't changed the font settings. The tool uses Type 42 fonts by default. If you modified `matplotlib` rcParams globally, reset them. |
| Scale bar is wrong length | Check the **Pixel Size (µm/px)** field. If it shows 0, the tool uses the manual mm fallback. Enter your microscope's pixel size manually (check your acquisition software's metadata). |
| Mega-montage is missing | The mega-montage is only generated when there are `*_SUMMARY_MONTAGE.png` files in the date folder. Make sure at least one experiment has been finalized. |
| Mega-montage shows only one column | If you have only 1–2 experiments, one column may be empty. The two-column layout balances images by height — with few images, they may all fit in one column. |

### General tips

- **Save your work often**: The session record is automatically saved on FINALIZE. You can also click **SAVE SESSION RECORD** at any time after finalization.
- **Use descriptive experiment names**: They appear on your output files and montages. E.g., `ATG5_KO_EBSS_2h` is more useful than `Experiment_1`.
- **Process one condition at a time**: Run the tool once per experimental condition (e.g., WT control, KO treatment). Each becomes a separate experiment folder under the same date folder.
- **Keyboard shortcuts**: **Enter** = Process & Next, **S** = Skip, **B** = Go Back. These work when the Napari viewer is focused.
- **If something goes wrong**: Click **↻ NEW EXPERIMENT** to reset everything and start fresh. Your previously saved outputs are not affected.

---

## License

This project is released under the [MIT License](LICENSE).

---

## Citation

If you use this tool in published research, please cite this repository:

```
StevenChung13. Confocal Colocalization Tool. GitHub. https://github.com/StevenChung13/confocal-colocalization-tool
```
