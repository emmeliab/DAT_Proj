# DAT_Proj

This repository contains the data and code to process 32 DAT and 27 steady-state A/Ci curves from 13 trees at K67 in the Tapajos National Forest in Para, Brazil in August 2022.

## Data
* A/Ci data were collected from 13 trees in the Tapajos National Forest accessible from the walkup tower at K67; 33 DAT curves and 28 SS curves (note that one pair of curves were not a matched set and were excluded from the analyses, totaling 32 DAT curves and 27 SS curves)
* Data were collected using the LI-6800 Portable Photosynthesis System
  * DAT curves were collected with the DAT_continuous program and SS curves were collected manually at each point
  * DAT curves were run right before SS curves on the same leaves; two leaves per tree (typically)
  * Of note, MNG revisited MACA1 on 8/16 and measured new, different leaves than on the first day (8/6?). MNG recorded a Leaf_number of 2, which had already been recorded in the data. In a manual recording of data, CDS gave a unique identifier to each leaf, with the second leaf 2 recorded as leaf 8.

## File Structure
### 0_Archive
Contains results and scripts from our first pass-over with the data. Contains redundant data files and MNG scripts for processing curves. Files in here are likely not needed but are being kept in case.

### 1_Raw_data
Contains raw .xlsx files from the LI-6800, with added Data_QC column, which was inputted manually

### 2_Cleaning_Scripts
Contains the script (1_Tapajos_DAT_data_assembly.R) used to clean and assemble the raw data files. 

- Takes raw .xlsx files as input and compiles to 3_Clean_data/raw_combined_aci_data.csv
- Removes outliers in the raw data and saves as 3_Clean_data/clean_aci_noOutliers.csv

### 3_Clean_data
Contains .csv of cleaned data as well as tree identification

- clean_aci_noOutliers.csv: All A/Ci data points with erroneous points and outliers removed
- clean_aci_with_uniquecode.csv : Contains the same data as clean_aci_noOutliers, but with manual fixes to the IDs
- id_codebook: contains the species names associated with each tree ID, as well as 
- raw_combined_aci_data.csv: composite of all the data in the .xlsx files with added Data_QC column, but no points removed
- rel_canopy_pos.csv: contains the relative canopy position for each tree
- unique_ids.csv : contains the unique tree ID, K67 ID, and species code for internal use

### 4_Fitting_stats_figs_scripts
Contains scripts for fitting curves, performing statistical analyses, and building figures for the manuscripts

- 2_Tapajos_Fit_ACi.R : script for fitting the A/Ci curves with 'photosynthesis'. Takes clean_aci_with_uniquecode.csv as inputs and outputs temperature-corrected fitted parameters and .RData files of the fits
- 3_Photo_Stat_Analysis.R: scripts for running statistical analyses on the fitted parameters
- 4_Figs_and_tables.R: scripts for making all in-text and supplementary tables and figures
- bootstrapping_estimates.R : function for running bootstrapping on linear mixed models; developed by Dusty Gannon

### 5_Results
Contains fitted parameters as well as .RData files of curve fits

- .RData files:
     - DAT_photo_fits_TPU.RData: fits for DAT curves with TPU enabled
     - DAT_photo_fits_noTPU.RData: fits for DAT curves without TPU enabled
     - SS_photo_fits_TPU.RData: fits for SS curves with TPU enabled
     - SS_photo_fits_noTPU.RData: fits for SS curves without TPU enabled

- .csv files
     - boot_res.csv : results from bootstrapping linear mixed models. Contains output from original model and bootstrapped results
     - DAT_photo_pars_crct_TPU.csv: the fitted parameters for the DAT curves with TPU enabled
     - DAT_photo_pars_crct_noTPU.csv: the fitted parameters for the DAT curves without TPU enabled
     - estimates_boot.csv : the output from linear mixed model bootstrapping from 500 iterations
     - lf_diffs_summ_TPU.csv: difference in SS and DAT fitted parameters with TPU enabled at the leaf level
     - lf_diffs_summ_noTPU.csv: difference in SS and DAT fitted parameters without TPU enabled at the leaf level
     - null_boot.csv : the null model output from bootstrapping over 500 iterations
     - pars_ecophys_noTPU.csv: the fitted parameters on the DAT curves without TPU fitting using 'plantecophys.' Note that the package was especially sensitive to noise in the data, and we elected not to proceed with this package.
     - pho_nd_stat.csv: relevant columns of DAT and SS curves fitted without TPU and without curves that displayed overshoot
     - pho_nd_stat_tpu.csv: relevant columns of DAT and SS curves fitted with TPU and without curves that display overshoot
     - pho_stat.csv: relevant columns of DAT and SS curves fitted without TPU
     - pho_stat_TPU.csv: relevant columns of DAT and SS curves fitted with TPU
     - SS_photo_pars_crct_noTPU.csv: the fitted parameters for the SS curves without TPU enabled
     - SS_photo_pars_crct_TPU.csv: the fitted parameters for the SS curves with TPU enabled
     - tree_diffs_summary_noTPU.csv: difference in SS and DAT fitted parameters without TPU enabled at the tree level
     - tree_diffs_summary_TPU.csv: difference in SS and DAT fitted parameters with TPU enabled at the tree level
     - tree_nOS_diffs_summary_TPU.csv: difference in SS and DAT fitted parameters with TPU enabled and without the curves with an overshoot at the leaf level
     - tree_nOS_diffs_summary_noTPU.csv: difference in SS and DAT fitted parameters without TPU enables and without the curves with overshoot at the leaf level
     - wilcox_table.csv : results from the tree-level Wilcoxon sign-rank tests
### 6_Figures
Contains plots of cleaned data, fitted figures, and figures for the manuscript
