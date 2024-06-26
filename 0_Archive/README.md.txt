# Some files in the archive:

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

