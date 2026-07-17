# 🔬 Colocalization Analysis of SatDNA Transcripts in *Drosophila virilis* Oocytes

<p align="center">
  <img width="336" height="329" alt="Dvir48_CTACcy3_green_CTATcy5_magenta_09112025 lif - Image105-1" src="https://github.com/user-attachments/assets/9b4b60af-32ae-4f23-aea3-016340934b6a" style="border-radius: 8px; margin-right: 10px;" />
  <img width="336" height="329" alt="CTACcy5green_CTATcy3magenta_sensboth_03102025 lif - Image104-1" src="https://github.com/user-attachments/assets/0589e3f7-a803-43cc-8bcb-a0a1424747e2" style="border-radius: 8px;" />
</p>

## 📖 Overview

This repository contains a comprehensive suite of scripts for analyzing the **2D and 3D colocalization** of transcripts within oocytes across various *Drosophila* species. 
The analysis pipeline integrates the [JACoP2](link-to-documentation) plugin and employs [Renyi's entropy](https://www.sciencedirect.com/science/article/pii/S0031320396000659) as a threshold method for robust and unbiased colocalization detection.

> 💡 **Note:** While these scripts are specifically tailored for the analyses presented in our article, they can be easily adapted for semi-automated colocalization studies on other fluorescence image datasets.

## 🛠️ Colocalization Scripts

The primary macro scripts calculate 2D and 3D colocalization metrics for transcript localization within oocytes. The pipeline explicitly relies on **Renyi's entropy** for threshold-setting, selected to provide statistically rigorous and reproducible colocalization measurements across complex biological samples.

## 🗂️ Repository Structure

```text
📦 Colocalisation_analysis
├── 📂 Fiji_colocalization_script/
│   ├── 📄 Macro_analyse_colocalisation_forZstack_only.ijm
│   ├── 📄 Macro_analyse_colocalization_forZstack.ijm
│   ├── 📄 Macro_colocalisation.ijm
│   └── 📄 README.md (Usage instructions for Fiji macros)
|
├── 📂 Script_colocalisation/
│   ├── 📄 Analyse_colocalisation_R_ovaire_3D.R
│   ├── 📄 README.md
│   ├── 📄 Traitement_data_colocalisation_sondeantisens_Dvir48_CTAC_CTAT.R
│   └── 📄 Traitement_data_colocalisation_Dvir48_sondesens_CTAC_CTAT.R
|
├── 📂 Script_score_occurence/
│   ├── 📄 Analyse__Dvir48_combinaison.R
│   └── 📄 Comparaison_CTAC.R
|
└── 📄 README.md
```

## 📊 Statistical Analysis

All subsequent statistical analyses of the image outputs are included in this repository. These were performed using:
- 🔵 **RStudio** for R-based statistical testing and data wrangling.

Feel free to review the code or adapt the analyses for your own experimental studies!

## 💾 Data Availability

Original image data are available upon request *(link to be provided)*. This ensures full reproducibility of our analyses. 
> ⚠️ **Important:** These scripts are optimized for specific image file naming conventions. You may need to modify the regular expressions or string parsing in the macros/scripts to accommodate different file naming schemes.

## 🤖 AI Assistance Disclaimer

This project leveraged **Gemini Pro 3** as an assistant for identifying critical errors and generating a specific key function within the scripts. Approximately **95%** of the codebase was developed independently by the authors, with AI assistance accounting for ~**5%** of the final implementation.

## 📝 Citation

If you use these scripts, methodologies, or data in your research, please cite our corresponding article:

> **Expression of AAACTAC satellite repeats as a long noncoding RNA in the early oocyte of *Drosophila virilis***

---

💬 For any questions, issues, or suggestions, please feel free to **open an issue** on this repository.
