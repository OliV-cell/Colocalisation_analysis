# Colocalization analysis of SatDNA transcript in Drosophila virilis oocyte

<p align="center">
<img width="336" height="329" alt="Dvir48_CTACcy3_green_CTATcy5_magenta_09112025 lif - Image105-1" src="https://github.com/user-attachments/assets/9b4b60af-32ae-4f23-aea3-016340934b6a" />
<img width="336" height="329" alt="CTACcy5green_CTATcy3magenta_sensboth_03102025 lif - Image104-1" src="https://github.com/user-attachments/assets/0589e3f7-a803-43cc-8bcb-a0a1424747e2" />
</p>

## Colocalization script
The script for the colocalization analysis has been developed to calculate 2D and 3D colocalization of transcripts within the oocyte across various Drosophila species. Furthermore, the [JACoP2](https://imagej.net/plugins/jacop) plugin has been used as the core tool for the analysis. We have also used the [DiAna](https://imagej.net/plugins/distance-analysis) plugin to measure, for example, the centroide-centroide distance of our 3D object, using a [CLIJ2](https://clij.github.io/) implementation with GPU acceleration during the labelling of our object. 
Finally, all colocalization analyses use the [Renyi's entropy](https://www-sciencedirect-com.acces.bibl.ulaval.ca/science/article/pii/S0031320396000659) as a threshold; indeed, we chose this method because it best captured our object of study. 
Those scripts are highly specialized for this article (Put the name); nevertheless, they can be used to implement semi-automated colocalization for any image data.

## Data availability
The image data are available here (link provided), so anyone can reproduce the output data themselves. It is important to note that these scripts are highly adapted to my image files' names; therefore, a tutorial on the colocalisation script is available here https://youtu.be/CHfANNardjs

## Statistical analysis
All executed statistical analyses have also been given. Feel free to use them if needed or to check the code. For the record, the analysis was done in RStudio and in Julia with the Antigravity IDE. 

## IA use disclaimer
I have been using Gemini Pro 3 as an assistant to find critical errors and generate a key function in my script. Therefore, only 5% of this script has been written by an IA assistant. 

## Usage of this material
If you find those scripts and data useful, please cite from: Expression of AAACTAC satellite repeats as a long noncoding RNA in the early oocyte of Drosophila virilis

