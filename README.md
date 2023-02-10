# DAT_Proj

## Purpose
This repository contains the data and analysis of the DAT/ACi curve data we collected in the Tapajos in Aug 2022

## Data
* ACi data were collected on 13 trees in the Tapajos National Forest in Brazil accessible from the walkup tower at K67
* Data were collected using the LI-6800
* DAT curves were ran right before Traditional Curves on the same leaves; two leaves per tree (typically)
* Of note, MNG revisted MACA1 on 8/16 and measured new, different leaves than on the first day (8/6?). MNG recorded a Leaf_number of 2, which had already been recorded in the data. In a manual recoding of data, CDS gave a unique identifier to each leaf, with the second leaf 2 recorded as leaf 8.

## Files
### Raw_data
Contains raw data files from the 6800, with added Data_QC column

### Inputs
Contains: 
* Fully assembled data file titled: clean_aci_data_one_file.csv. 
* A table of unique identifiers: unique_ids
* Fully Assembled data file + unique identifiers: clean_aci_with_uniquecode.csv
* DAT data that filters out the backwards points and three curves that did not fit with fitacis: filtered_DAT_data.csv
* All aci data with outliers removed: Aci_no_out.csv

### Figures
Contains:
* Scatterplots of DAT and Trad curves
* Curve fits from plantecophys package
* Curve fits from photosynthesis package

### Results
Contains results of A-ci.adj.MG.R
* Fitted curves for each curve without back-filtering: A_ci_fit_DAT_Tapajos2.pdf
* Parameters for the curve fits without back-filtering: A_ci_fit_DAT_Tapajos2.csv
* Fitted curves for each DAT curve with back-filtering: A_ci_fit_DAT_Tapajos2_noback.pdf
* Parameters for the DAT curve fits with back-filtering: A_ci_fit_DAT_Tapajos2_noback.csv
* Compliled parameters of non-filtered and back-filtered data: curve_fitting_MG_out.csv

Contains results of fits from photosynthesis package
* Parameters of Trad curves: trad_fits_photo_pars.csv
* Parameters of DAT curves: dat_fits_photo_pars.csv

### Scripts
#### Tapajos_Data_Assembly.R
This R file assembles all of the individual clean data files (in the Inputs Folder), filters by quality-controlled data, and outputs in to a data file titled clean_aci_data_one_file.csv

#### Tapajos_Fit_Aci.R
Cleans outliers and contains script for back-filtering data. Fits curves with plantecophys and photosynthesis packages. Takes clean_aci_with_uniquecode.csv and Aci_no_out.csv as inputs. Outputs plots of fitted curves into Figures Folder and csvs of parameters in the Results folder. Previously charlie_playing.R

#### A-ci.adj.MG.R
Maquelle's code for fitting aci curves. Outputs pdf of plots and csv file of values into Results folder

#### MG_fit_stats_and_figs.R
Charlie's code for comparing non-filtered and back-filtered data, and plots for comparing Vcmax, Jmax, and TPU between traditional and DAT curves, as calculated by Maquelle's code

#### Tapajos_ACi_stat_analysis.R
Contains normality tests and preliminary stats for analyzing parameters calculated through plantecophys (Tapajos_Fit_Aci.R)
