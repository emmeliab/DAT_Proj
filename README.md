# DAT_Proj

## Purpose
### This repository contains the data and analysis of the DAT/ACi curve data we collected in the Tapajos in Aug 2022

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
Contains plots of aci curves for each curve fitted with plantecophys

### Results
Contains results of A-ci.adj.MG.R
* Fitted curves for each curve: A_ci_fit_DAT_Tapajos2.pdf
* Values for the curve fits: A_ci_fit_DAT_Tapajos2.csv

### Tapajos_Data_Assembly.R
This R file assembles all of the individual clean data files (in the Inputs Folder), filters by quality-controlled data, and outputs in to a data file titled clean_aci_data_one_file.csv

### Tapajos_Fit_Aci
Takes clean_aci_data_one_file.csv as inputs and outputs plots of aci curves from each species (Outputs folder) as well as fitting data for each leaf (superseeded)

### charlie_playing.R
Cleans, plots, and fits aci curves with plantecophys. This has now become our 'main' script for running stats. Plans to merge this into Tapajos_Fit_Aci later. Takes clean_aci_with_uniquecode.csv as inputs. Filters out outliers and "backwards" datapoints. Outputs curve fit plots into Figures Folder.

### A-ci.adj.MG.R
Maquelle's code for fitting aci curves. Outputs pdf of plots and csv file of values into Results folder
