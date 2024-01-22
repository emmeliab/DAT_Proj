# DAT_Proj

## Purpose
This repository contains the data and code to process 32 DAT and 27 steady-state A/Ci curves during a field campaign to the Tapajos National Forest in Para, Brazil during August of 2022.

## Data
* A/Ci data were collected on 13 trees in the Tapajos National Forest accessible from the walkup tower at K67; 33 DAT curves and 28 SS curves (note that one pair of curves were not a matched set and were excluded from the analyses, totalling 32 DAT curves and 27 SS curves)
* Data were collected using the LI-6800
* DAT curves were collected with the DAT_continuous program and SS curves collected manually at each point
* DAT curves were ran right before SS curves on the same leaves; two leaves per tree (typically)
* Of note, MNG revisted MACA1 on 8/16 and measured new, different leaves than on the first day (8/6?). MNG recorded a Leaf_number of 2, which had already been recorded in the data. In a manual recoding of data, CDS gave a unique identifier to each leaf, with the second leaf 2 recorded as leaf 8.

## File Structure
### 0_Archive
Contains results and scripts from our first pass-over with the data. Contains redundant data files and MNG scripts for processing curves. Files in here are likely not needed, but are being kept in case.

### 1_Raw_data
Contains raw data files from the 6800, with added Data_QC column

### 2_Cleaning_Scripts
Contains the script used to clean and assemble the raw data files

### 3_Clean_data
Contains .csv of cleaned data as well as tree identification

### 4_Fitting_stats_figs_scripts
Contains scripts for fitting curves, performing statistical analyses, and building figures for the manuscripts

### 5_Results
Contains fitted parameters as well as .RData files of curve fits

### 6_Figures
Contains plots of cleaned data, fitted figures, and figures for the manuscript



# Old stuff - delete later
### Results
Contains results of A-ci.adj.MG.R
* Fitted curves for each curve without back-filtering: A_ci_fit_DAT_Tapajos2.pdf
* Parameters for the curve fits without back-filtering: A_ci_fit_DAT_Tapajos2.csv
* Fitted curves for each DAT curve with back-filtering: A_ci_fit_DAT_Tapajos2_noback.pdf
* Parameters for the DAT curve fits with back-filtering: A_ci_fit_DAT_Tapajos2_noback.csv
* Compliled parameters of non-filtered and back-filtered data: curve_fitting_MG_out.csv

Contains results of fits from photosynthesis package
* Parameters of Trad curves with TPU: trad_fits_photo_pars.csv
* Parameters of Trad curves without TPU: trad_fits_photo_pars_no_TPU.csv
* Parameters of Trad fits without TPU with temperature correction: trad_fits_photo_pars_correct_no_TPU.csv
* Parameters of DAT curves with TPU: dat_fits_photo_pars.csv
* Parameters of DAT curves with backwards point removed: dat_fit_ex_photo_pars.csv
* Parameters of DAT curves without TPU: dat_fits_photo_pars_filt_no_TPU.csv
* Parameters of DAT fits without TPU with temperature correction: dat_fits_photo_pars_filt_correct_no_TPU.csv

Contains results of fits from plantecophys package
* Parameters with TPU: params_ecophys.csv
* Parameters without TPU: params_ecophys_no_TPU.csv

### Scripts
#### 1_Tapajos_Data_Assembly.R
This R file assembles all of the individual clean data files (in the Inputs Folder), filters by quality-controlled data, and outputs in to a data file titled clean_aci_data_one_file.csv

#### 2_Tapajos_Fit_Aci.R
Cleans outliers and contains script for back-filtering data. Fits curves with plantecophys and photosynthesis packages. Performs a temperature correction on photosynthesis results. Takes clean_aci_with_uniquecode.csv and Aci_no_out.csv as inputs. Outputs plots of fitted curves into Figures Folder and csvs of parameters in the Results folder. Previously charlie_playing.R

#### 3_Tapajos_ACi_stat_analysis.R
Takes parameters from MG code, photosynthesis, and plantecophys as inputs. Compares Vcmax, Jmax, and TPU results between packages. Performs statistical comparisons between DAT and Trad curves for each package and does preliminary visualizatiom. Outputs figures to Figures folder.

#### A_ci_oldVersion.R
Maquelle's A/Ci fitting code containing fitting for TPU. Takes Aci_no_out.csv as inputs. Outputs a pdf of fitted curves to Figures folder and a csv of parameters to Results folder.

#### A_ci_in_progress.R
Maquelle's A/Ci fitting code with TPU fitting removed. Takes Aci_no_out.csv as inputs. Outputs a pdf of fitted curves to Figures folder and a csv of parameters to Results folder.

#### temp_correct_photo.R
Charlie's code for performing temperature correction on photosynthesis results. This has been incorporated into 2_Tapajos_Fit_Aci.R. Plans to delete this later.

