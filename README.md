# DAT_Proj

## Purpose
This repository contains the data and analysis of the DAT/ACi curve data we collected in the Tapajos in Aug 2022

## Data
* ACi data were collected on 13 trees in the Tapajos National Forest in Brazil accessible from the walkup tower at K67
* Data were collected using the LI-6800
* DAT curves were ran right before Traditional Curves on the same leaves; two leaves per tree (typically)
* Of note, MNG revisted MACA1 on 8/16 and measured new, different leaves than on the first day (8/6?). MNG recorded a Leaf_number of 2, which had already been recorded in the data. In a manual recoding of data, CDS gave a unique identifier to each leaf, with the second leaf 2 recorded as leaf 8.

## Files
### Inputs
Contains: 
* Clean data files from the LI-6800
* Fully assembled data file titled: clean_aci_data_one_file.csv. 
* A table of unique identifiers: unique_ids
* Fully Assembled data file + unique identifiers: clean_aci_with_uniquecode.csv
* DAT data that filters out the backwards points and three curves that did not fit with fitacis: filtered_DAT_data.csv
* All aci data with outliers removed: Aci_no_out.csv

### Figures
Contains plots of A/Ci curves

### Results
Contains a pdf of plots and csv file of the values resulting from A-ci.adj.MG.R
* Curve fit of each plot: A_ci_fit_DAT_Tapajos.pdf
* Output values from curve fitting: A_ci_fit_DAT_Tapajos.csv

### Scripts
#### Tapajos_Data_Assembly.R
This R file assembles all of the individual clean data files (in the Inputs Folder), filters by quality-controlled data, and outputs in to a data file titled clean_aci_data_one_file.csv

#### Tapajos_Fit_Aci.R
Takes clean_aci_data_one_file.csv as inputs and outputs plots of aci curves from each species (Outputs folder) as well as fitting data for each leaf (redundant)

#### charlie_playing.R
This has now become our 'main' script for running stats. Plans to merge this into Tapajos_Fit_Aci later. Takes clean_aci_with_onecode.csv as inputs. Removes outliers and the "backwards" portion of the curve and fits curves with plantecophys. Outputs plots of the A/Ci curves to the Figures folder and csv files with outliers removed to the Inputs folder. 

#### A-ci.adj.MG.R
Maquelle's code for fitting A/Ci curves. Outputs a pdf of the plots (A_ci_fit_DAT_Tapajos.pdf) and a csv file of the values (A_ci_fit_DAT_Tapajos.txt) to the Results folder
