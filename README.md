# DAT_Proj

## Purpose
This repository contains the data and analysis of the DAT/ACi curve data we collected in the Tapajos in Aug 2022

## Data
* ACi data were collected on 13 trees in the Tapajos National Forest in Brazil accessible from the walkup tower at K67
* Data were collected using the LI-6800
* DAT curves were ran right before Traditional Curves on the same leaves; two leaves per tree (typically)
* Of note, MNG revisted MACA1 on 8/16 and measured new, different leaves than on the first day (8/6?). MNG recorded a Leaf_number of 2, which had already been recorded in the data. In a manual recoding of data, CDS gave a unique identifier to each leaf, with the second leaf 2 recorded as leaf 8.

## Files
<<<<<<< HEAD
=======
### Raw_data
Contains raw data files from the 6800, with added Data_QC column

>>>>>>> b67e019a731fc6b4722a3ff7b2875ccb271bb7f9
### Inputs
Contains: 
* Fully assembled data file titled: clean_aci_data_one_file.csv. 
* A table of unique identifiers: unique_ids
* Fully Assembled data file + unique identifiers: clean_aci_with_uniquecode.csv
* DAT data that filters out the backwards points and three curves that did not fit with fitacis: filtered_DAT_data.csv
* All aci data with outliers removed: Aci_no_out.csv

### Figures
<<<<<<< HEAD
Contains plots of A/Ci curves

### Results
Contains a pdf of plots and csv file of the values resulting from A-ci.adj.MG.R
* Curve fit of each plot: A_ci_fit_DAT_Tapajos.pdf
* Output values from curve fitting: A_ci_fit_DAT_Tapajos.csv
=======
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
>>>>>>> b67e019a731fc6b4722a3ff7b2875ccb271bb7f9

### Scripts
#### Tapajos_Data_Assembly.R
This R file assembles all of the individual clean data files (in the Inputs Folder), filters by quality-controlled data, and outputs in to a data file titled clean_aci_data_one_file.csv

#### Tapajos_Fit_Aci.R
<<<<<<< HEAD
Takes clean_aci_data_one_file.csv as inputs and outputs plots of aci curves from each species (Outputs folder) as well as fitting data for each leaf (redundant)

#### charlie_playing.R
This has now become our 'main' script for running stats. Plans to merge this into Tapajos_Fit_Aci later. Takes clean_aci_with_onecode.csv as inputs. Removes outliers and the "backwards" portion of the curve and fits curves with plantecophys. Outputs plots of the A/Ci curves to the Figures folder and csv files with outliers removed to the Inputs folder. 

#### A-ci.adj.MG.R
Maquelle's code for fitting A/Ci curves. Outputs a pdf of the plots (A_ci_fit_DAT_Tapajos.pdf) and a csv file of the values (A_ci_fit_DAT_Tapajos.txt) to the Results folder
=======
Cleans outliers and contains script for back-filtering data. Fits curves with plantecophys and photosynthesis packages. Takes clean_aci_with_uniquecode.csv and Aci_no_out.csv as inputs. Outputs plots of fitted curves into Figures Folder and csvs of parameters in the Results folder. Previously charlie_playing.R

#### A-ci.adj.MG.R
Maquelle's code for fitting aci curves. Outputs pdf of plots and csv file of values into Results folder

#### MG_fit_stats_and_figs.R
Charlie's code for comparing non-filtered and back-filtered data, and plots for comparing Vcmax, Jmax, and TPU between traditional and DAT curves, as calculated by Maquelle's code

#### Tapajos_ACi_stat_analysis.R
Contains normality tests and preliminary stats for analyzing parameters calculated through plantecophys (Tapajos_Fit_Aci.R)
>>>>>>> b67e019a731fc6b4722a3ff7b2875ccb271bb7f9
