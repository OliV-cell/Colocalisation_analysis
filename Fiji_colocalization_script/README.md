# Colocalization Macro Tutorial

## Video Tutorial

A comprehensive tutorial for the basic 2D colocalization macro is available here: [Colocalization Analysis Tutorial](https://youtu.be/CHfANNardjs)

This tutorial demonstrates the default pipeline. Feel free to modify the macro to suit your specific experimental needs.

## Getting Started

### Prerequisites
- Fiji/ImageJ with macro scripting support
- Image files in TIFF or compatible format
- Two or more fluorescence channels for colocalization analysis

### Running the Macro

1. Open Fiji/ImageJ
2. Load your image file
3. Go to `Macros` → `Edit Macros` (or open the macro file directly)
4. Run the macro using `Macros` → `Run Macro`

### Batch Processing (Optional)

The macro can be run in batch mode for improved performance:

```
setBatchMode(true);
```

This line is commented out by default to allow visualization during the tutorial, but you can uncomment it to disable screen updates and speed up processing when running multiple images.

## Macro Versions

### 2D Colocalization Macro

The basic 2D macro processes two-dimensional images and calculates colocalization metrics including:
- Pearson correlation coefficient
- Manders' coefficients
- Intensity correlation quotient

**Output files required:**
- `Coloc_score/` – Numerical colocalization scores
- `Coloc_graph/` – Visualization graphs
- `Tiff_bin/` – Binary thresholded images

### 3D Colocalization Macro

The 3D macro operates on the same logic as the 2D version but processes volumetric (Z-stack) data. The key difference is the handling of additional Z-plane data throughout the analysis pipeline.

**Additional requirements for 3D analysis:**
- Three additional files per category to store Z-stack image data:
  - `Coloc_score/` – Numerical colocalization scores for the stack
  - `Coloc_graph/` – Visualization graphs for the stack
  - `Tiff_bin/` – Binary thresholded images for each Z-plane

**File structure:** You will need 2 directories (one for 2D results and one for 3D results) for each analysis category (Coloc_score, Coloc_graph, Tiff_bin).

## Output

The macros generate:
- **Coloc_score**: Quantitative colocalization metrics stored as text files or spreadsheets
- **Coloc_graph**: Graphical representations of colocalization patterns
- **Tiff_bin**: Tiff image data

## Customization

Both 2D and 3D macros are fully customizable:
- Adjust threshold values for your specific dyes and imaging conditions
- Modify output file naming conventions
- Include additional image analysis parameters as needed

## Notes

- The macros are designed to work with standard fluorescence microscopy data
- Ensure consistent image preprocessing (background subtraction, alignment, etc.) before running the analysis
- Results depend on proper channel registration and imaging parameters
