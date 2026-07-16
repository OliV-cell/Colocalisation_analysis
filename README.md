# Colocalization Analysis of SatDNA Transcripts in *Drosophila virilis* Oocytes

<p align="center">
<img width="336" height="329" alt="Dvir48_CTACcy3_green_CTATcy5_magenta_09112025 lif - Image105-1" src="https://github.com/user-attachments/assets/9b4b60af-32ae-4f23-aea3-016340934b6a" />
<img width="336" height="329" alt="CTACcy5green_CTATcy3magenta_sensboth_03102025 lif - Image104-1" src="https://github.com/user-attachments/assets/0589e3f7-a803-43cc-8bcb-a0a1424747e2" />
</p>

## Overview

This repository contains a comprehensive suite of scripts for analyzing the 2D and 3D colocalization of transcripts within oocytes across various *Drosophila* species. The analysis pipeline integrates the [JACoP2](link-to-documentation) plugin and employs [Renyi's entropy](https://www-sciencedirect-com.acces.bibl.ulaval.ca/science/article/pii/S0031320396000659) as a threshold method for robust colocalization detection.

While these scripts are specifically tailored for the analyses presented in our article, they can be adapted for semi-automated colocalization studies on other image datasets.

## Colocalization Scripts

The main colocalization analysis script calculates 2D and 3D colocalization metrics for transcript localization within oocytes. The pipeline uses Renyi's entropy as a threshold-setting method, which was selected to provide statistically rigorous colocalization measurements.

## Reperotory structure


## Data Availability

Image data are available upon request (link to be provided). This ensures full reproducibility of the analyses. Please note that these scripts are optimized for specific image file naming conventions; they may require modification to work with different file naming schemes.

## Statistical Analysis

All statistical analyses are included in this repository. The analyses were performed using:
- **RStudio** for R-based statistical analyses
- **Julia** with the Antigravity IDE for advanced computational analyses

Feel free to review the code or adapt the analyses for your own studies.

## AI Assistance Disclaimer

This project used Gemini Pro 3 as an assistant for identifying critical errors and generating a key function within the scripts. Approximately 95% of the codebase was developed independently, with AI assistance accounting for ~5% of the final implementation.

## Citation

If you use these scripts or data in your research, please cite:

**Expression of AAACTAC satellite repeats as a long noncoding RNA in the early oocyte of *Drosophila virilis***

---

For questions or issues, please feel free to open an issue on this repository.
