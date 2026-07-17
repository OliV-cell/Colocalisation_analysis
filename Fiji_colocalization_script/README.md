# 🔬 Fiji Colocalization Analysis Macros

Welcome to the Fiji Colocalization Analysis Macro suite. This toolkit provides semi-automated, robust scripts for computing colocalization metrics (2D and 3D) within Fiji/ImageJ, specifically tailored for fluorescence microscopy.

## 📺 Video Tutorial

A comprehensive step-by-step tutorial for the basic 2D colocalization macro is available on YouTube: 

[![Colocalization Analysis Tutorial](https://img.youtube.com/vi/CHfANNardjs/0.jpg)](https://youtu.be/CHfANNardjs)

> **Note:** The tutorial demonstrates the default pipeline. You are encouraged to modify the macro parameters to suit your specific experimental conditions and imaging modalities.

## 🚀 Getting Started

### Prerequisites
- **Software:** Fiji/ImageJ with macro scripting support enabled.
- **Data Format:** Image files in `.tif` or any bio-formats compatible format.
- **Channels:** At least two fluorescence channels to perform colocalization analysis.

### Running the Macro
1. Launch **Fiji/ImageJ**.
2. Load your image file.
3. Navigate to `Plugins` → `Macros` → `Edit...` and open your desired `.ijm` macro file.
4. Execute the script by clicking `Run` (or `Macros` → `Run Macro`).

### Batch Processing (Optional)
For high-throughput analysis, the macro can be executed in batch mode to significantly improve performance. 
Locate the following line in the script:
```java
// setBatchMode(true);
```
Uncomment it (remove the `//`) to disable screen updates during processing:
```java
setBatchMode(true);
```
*This is disabled by default to allow visual feedback during the tutorial.*

## 🧬 Macro Versions

### 2D Colocalization Macro (`Macro_colocalisation.ijm`)
Designed for standard two-dimensional images. Calculates the following colocalization metrics:
- **Pearson's Correlation Coefficient (PCC)**
- **Manders' Overlap Coefficients (M1 & M2)**
- **Intensity Correlation Quotient (ICQ)**

### 3D Colocalization Macro (Z-Stack)
The 3D macros (`Macro_analyse_colocalization_forZstack.ijm` & `Macro_analyse_colocalisation_forZstack_only.ijm`) operate on volumetric (Z-stack) data. They extend the 2D logic across additional Z-planes.

## 📁 Output & File Structure

The macros will generate output across three main directories. 
> **Important:** For 3D analysis, you will need separate directories for 2D (MIP or single slice) and 3D results for each category below.

- 📊 **`Coloc_score/`**: Quantitative metrics stored as text/spreadsheet files.
- 📈 **`Coloc_graph/`**: Graphical representations of colocalization patterns (scatter plots, etc.).
- 🖼️ **`Tiff_bin/`**: Binary thresholded `.tif` images used during the calculation.

## ⚙️ Customization

These macros are open-source and fully customizable:
- **Thresholds:** Adjust thresholding algorithms (e.g., Renyi's entropy) to match your specific dyes and signal-to-noise ratios.
- **Naming Conventions:** Modify output file names for better integration into your data management pipelines.
- **Additional Parameters:** Include more morphological or intensity-based measurements as needed.

## ⚠️ Important Notes
- Ensure **consistent image preprocessing** (e.g., background subtraction, channel alignment/registration, deconvolution) before running the analysis.
- The results heavily depend on **proper channel registration** and appropriate imaging parameters during acquisition.
