#!/bin/bash
# This script lists the commands to reproduce the key analyses done in the manuscript
# A minimal set of data to reproduce the analyses in part 1 and 2 is provided in GitHub
# Random numbers generated may be different
# The exact results generated by the authors are provided in Dryad repository (see GitHub site's description file)

# Part 1. Fitting pulsed evolution models and preparing reference for prediction
# Fit homogeneous BM and PE models
# Starting data:  Reference/reference.tre
#                 Reference/16S_GCN.txt
# Expected output:  Reference/homogeneous_data.RDS
#                   Reference/homogeneous_models.RDS
Rscript fit_homogeneous_models.R Reference/reference.tre Reference/16S_GCN.txt Reference/homogeneous

# Binary partition of the phylogeny
Rscript binary_partition_by_AIC.R Reference/homogeneous_data.RDS Reference/homogeneous_models.RDS Reference/binary_partition.RDS Reference/rescaled_data_model.RDS
# Prepare prediction reference
Rscript prepare_prediction_data.R Reference/rescaled_data_model.RDS Reference/prepared_reference.RDS

# Part 2. Cross-validation at different NSTD cutoff
# Starting data: results from Part 1.
# Expected output:  CV/GCN.PE.CV.RDS
#                   CV/GCN.BM.CV.RDS
#                   CV/GCN.MP_EMP.CV.RDS

# Part 3. Predicting GCN for SILVA and HMP

# Part 4. Microbiome analyses with GCN correction

# Part 5. Figures and Tables
# Figure 1 and 2
# Starting data:  Reference/homogeneous_data.RDS
#                 CV/GCN.PE.CV.RDS
#                 CV/GCN.BM.CV.RDS
#                 CV/GCN.MP_EMP.CV.RDS
#                 CV/Baseline.trait.RDS
#                 CV/Baseline.PE.CV.RDS
#                 CV/Baseline.BM.CV.RDS
#                 CV/Baseline.MP_EMP.CV.RDS
#                 CV/CV.NSTD.RDS
# Expected output:  Fig_1.png Fig_1.pdf Fig_2.png Fig_2.pdf
Rscript Figure_1_Expectation.R
Rscript Figure_2_GCN_Classification.R

# Figure 3
# Starting data:
# Expected output:  
Rscript Figure_3_Abundance.R
# Figure 4
Rscript Figure_4_Beta_diversity_Bray.R

# Figure 5
# Starting data: EBI/EBI.adj.NSTI.txt
# Expected output: Fig_5.png Fig_5.pdf
Rscript Figure_5_EBI_NSTI_Biomes.R

# Supplementary Figure 1
# Starting data:  Reference/homogeneous_data.RDS
#                 Reference/taxids.RDS
#                 Reference/lineage_table.RDS
# Expected output Fig_S1.png Fig_S1.pdf
Rscript Figure_S1_Rate_heterogeneity.R

# Supplementary Figure 2
# Starting data:  Reference/prepared_reference.RDS
#                 Reference/taxids.RDS
#                 Reference/lineage_table.RDS
# Expected output Fig_S2.pdf
Rscript Figure_S2_Reference_insertion.R

# Supplementary Figure 3
# Starting data:  Reference/reference.tre
#                 Reference/taxids.RDS
#                 Reference/lineage_table.RDS
# Expected output Fig_S3.png Fig_S3.pdf
Rscript Figure_S3_distance_between_taxonomic_group.R

# Table 1

# Table 2

