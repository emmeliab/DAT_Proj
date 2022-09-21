# DAT_Proj

# Purpose
## This repository contains the data and analysis of the DAT/ACi curve data we collected in the Tapajos in Aug 2022

# Data
* ACi data were collected on 13 trees in the Tapajos National Forest in Brazil accessible from the walkup tower at K67
* Data were collected using the LI-6800
* DAT curves were ran right before Traditional Curves on the same leaves; two leaves per tree (typically)
* Of note, MNG revisted MACA1 on 8/16 and measured new, different leaves than on the first day (8/6?). MNG recorded a Leaf_number of 2, which had already been recorded in the data. In a manual recoding of data, CDS gave a unique identifier to each leaf, with the second leaf 2 recorded as leaf 8.

# Files
## Inputs
Contains Clean data files from the LI-6800 as well as the fully assembled data file titled: clean_aci_data_one_file.csv

## Tapajos_Data_Assembly.R
This R file assembles all of the individual clean data files (in the Inputs Folder), filters by quality-controlled data, and outputs in to a data file titled clean_aci_data_one_file.csv

## Tapajos_Fit_Aci
Takes clean_aci_data_one_file.csv as inputs and outputs plots of aci curves from each species (Outputs folder) as well as fitting data for each leaf (in progress)

## charlie_playing.R
This is Charlie messing with our DAT data to get real numbers :)