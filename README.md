# Akdeniz Tool for HU Calculation

Version 1.0.0

Desktop tool for measuring Hounsfield units (HU) in a three-dimensional region of interest (3D ROI)
placed on a spinal CT series. It reports the mean HU, median HU, the Gaussian-based FWHM interval and
FWHM mean HU of the ROI voxels, and exports a summary figure suitable for publication.

Developed by Muhammad Ateya, Department of Neurosurgery, Akdeniz University.

> \*\*Intended use.\*\* Research tool. It has not been validated for clinical decision-making.

## Features

* Loads a folder of DICOM slices and converts stored values to HU using RescaleSlope and RescaleIntercept.
* Axial, coronal and sagittal views with linked crosshair, drawn with voxel-spacing-aware proportions.
* 3D ROI defined by eight points placed on sagittal slices.
* ROI statistics: voxel count, mean HU, median HU, Gaussian-based FWHM interval, FWHM mean HU, minimum and maximum HU.
* Results are copied to the clipboard as tab-separated text for pasting into a spreadsheet.
* Optional 3D point-cloud preview of the ROI region (visualization only).
* Publication figure export (PNG and PDF): sagittal and coronal slices with the ROI outline, and the ROI HU distribution.

## Requirements

* Python 3 with Tkinter (the figure export was tested with Python 3.13 and matplotlib 3.11; other versions have not been tested).
* numpy, pydicom, Pillow, matplotlib, scipy (see `requirements.txt`).
* Some compressed DICOM transfer syntaxes need additional pydicom pixel-data handlers (for example
`pylibjpeg` or `python-gdcm`); these are not installed by `requirements.txt`.

## Installation

```bash
git clone https://github.com/mateya86/akdeniz-tool-for-hu-calculation.git
cd akdeniz-tool-for-hu-calculation
pip install -r requirements.txt
```

## Usage

```bash
python hounsfield\_calculator.py
```

1. Click **Open Spinal DICOM Folder** and select a folder containing a single CT series.
2. Scroll over a view to change slice. Left-click a view to move the crosshair.
3. In the **sagittal** view (right panel), right-click to place eight points in the order shown on screen:
four points on the first sagittal slice (Posterior Superior, Anterior Superior, Anterior Inferior,
Posterior Inferior), then scroll to another sagittal slice and place the same four on that slice.
The step labels read "LEFT" for the first four and "RIGHT" for the last four.
4. Press **Enter** to compute the ROI statistics. The result is copied to the clipboard and shown in a new window.
5. Press **P** to export the publication figure (PNG and PDF).

### Keyboard and mouse controls

|Input|Action|
|-|-|
|Left-click|Move crosshair|
|Right-click (sagittal view)|Place ROI point|
|Scroll wheel|Change slice in the view under the cursor|
|Backspace or U|Undo last point|
|C|Clear all points|
|V|Toggle 3D preview|
|Z|Toggle zoom mode (points cannot be placed while active)|
|Enter|Analyze (requires eight points)|
|P|Export publication figure (requires eight points)|
|N|Load a new subject|

## Methods

**HU conversion.** Stored pixel values are converted to HU: HU = RescaleSlope × pixel value + RescaleIntercept.

**Voxel spacing.** In-plane spacing comes from PixelSpacing. Slice spacing is the median difference between
consecutive ImagePositionPatient z-positions, because SliceThickness can differ from the reconstruction
increment. Spacing is used for display proportions and scale bars only; the ROI is defined in voxel indices.

**3D ROI.** For each instrumented vertebra, four boundary points are placed on each of two sagittal slices. The corresponding vertices are linearly interpolated across the intervening slices to define a three-dimensional region of interest (ROI), and all enclosed voxels are included in the HU calculations. The boundary points should be manually positioned to include as much accessible trabecular bone as possible while avoiding the cortical shell, endplates, posterior vertebral wall, focal lesions, osteophytes, and imaging artifacts.

**Gaussian-based FWHM interval.** The Gaussian-based FWHM interval is defined as the overall ROI mean ±1.177 standard deviations, corresponding to mean ±sqrt(2 ln 2) × SD. **FWHM mean HU** is defined as the mean of ROI voxels within this interval. This operational definition does not require the observed HU distribution to follow a Gaussian distribution.

**Density cloud.** The 3D point cloud is for visualization only. The displayed crop is regularly subsampled
and each HU bin is randomly thinned, so the cloud can differ slightly between redraws. This does not affect
any reported statistic.

## Outputs

**Clipboard (tab-separated):** Voxel Count, Overall Mean, Overall Median, FWHM Mean, FWHM Min, FWHM Max, Min, Max
(FWHM Min and FWHM Max are the limits of the Gaussian-based FWHM interval).

**Publication figure** (300 dpi PNG, plus a PDF with the same name): white background, Arial (Liberation Sans or
DejaVu Sans if Arial is not installed), text of at least 8 pt, colorblind-safe colors. Panel A shows sagittal and
coronal CT slices (display window -200 to 1000 HU) with the 3D ROI outline and a 10 mm scale bar. Panel B shows one
histogram of the ROI voxel HU with lines for mean, median and FWHM mean HU and the shaded Gaussian-based FWHM interval.
The voxel count and HU range are printed below the plot. No DICOM header fields are read or drawn, so no patient
identifier appears in the figure.

## Notes

* Optional: place a `background.png` next to the script to show a logo on the welcome window.
* Folder selection loads every readable DICOM file in the folder; use a folder with a single series.

## Citation

Please cite this software using the DOI below (see also `CITATION.cff`).

```
Ateya M. Akdeniz Tool for HU Calculation, version 1.0.0. Zenodo. https://doi.org/10.5281/zenodo.24043853
```


## License

Licensed under the Apache License, Version 2.0. See `LICENSE` and `NOTICE`.

