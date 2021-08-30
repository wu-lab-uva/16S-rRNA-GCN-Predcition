# 16S-rRNA-GCN-Predcition
This site contains all the necessary scripts and a minimal set of data to reproduce the results (with different random numbers generated) in

`Accounting for 16S rRNA copy number prediction uncertainty and its implications in microbial diversity analyses`, in prep.

The exact results generated in the study can be found on Dryad once the manuscript goes to publication:

https://doi.org/10.5061/dryad.2rbnzs7p5

(Before that happens, a private link can be provided upon request.)

## System requirements
The scripts have been tested on Windows 10, macOS 11.4 (Big Sur) and Ubuntu 18.04.1 LTS with R 3.6.3. 
The following package is required to run all the scripts: `RasperGade16S`,`RasperGade`,`castor`,`ape`,`phyloseq`,`ggplot2`,`ggpubr`,`ggtree`

## Reproduction of results in the preprint
To reproduce the results in the preprint on macOS or Linux systems, follow the instrucitons in 

`Microbiome/Instruction_Ubuntu.sh`

For Windows users, the reproduction of all figures can be done by sourcing the scripts under Scripts/Figures_and_Tables/ in RStudio.
