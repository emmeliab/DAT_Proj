# Load Packages
library(tidyverse)
library(readxl)
library(here)


# Load Data from Directory ------------------------------------------------
aci_files_nm <- list.files(path = here("1_Raw_data"),
                           pattern = "walkup_aci_clean")
print(aci_files_nm)

for (i in 1:length(aci_files_nm)) {
  files1 <- read_xlsx(path = paste0(here("1_Raw_data"), "/", aci_files_nm[i]),
                      sheet = "Measurements", skip = 14)
  sliced <- slice(files1, -(1))
  assign(paste0("data_", aci_files_nm[i]), sliced)                     
}





# Assemble into one data file ---------------------------------------------

# ## Make sure all the data have the same columns (in this case, the 6th and 7th are missing
# ## some columns)
# Reduce(setdiff, list(colnames(`data_2022-08-06-1020_walkup_aci_clean.xlsx`), 
#         colnames(`data_2022-08-07-0941_walkup_aci_clean.xlsx`),
#         colnames(`data_2022-08-11-0904_walkup_aci_clean.xlsx`),
#         colnames(`data_2022-08-11-0949_walkup_aci_clean.xlsx`),
#         colnames(`data_2022-08-11-1050_walkup_aci_clean.xlsx`),
#         colnames(`data_2022-08-12-0928_walkup_aci_clean.xlsx`),
#         colnames(`data_2022-08-13-1007_walkup_aci_clean.xlsx`),
#         colnames(`data_2022-08-14-0906_walkup_aci_clean.xlsx`),
#         colnames(`data_2022-08-08-0854_walkup_aci_clean.xlsx`),
#         colnames(`data_2022-08-09-0911_walkup_aci_clean.xlsx`),
#         colnames(`data_2022-08-10-0854_walkup_aci_clean.xlsx`),
#         colnames(`data_2022-08-16-walkup_aci_clean.xlsx`))) #finds the unique columns


## Make a function to remove all the stability criteria columns, since they are uneven 
## across files
rmv_stbl <- function(file_name) {
  file_name <- file_name[,!grepl(":", colnames(file_name))]
}


## Rename the different columns and remove stability columns
`data_2022-08-06-1020_walkup_aci_clean.xlsx` <- `data_2022-08-06-1020_walkup_aci_clean.xlsx` %>% 
  rename("Leaf_position" = "Leaf_Position") %>% 
  rename("Data_QC" = "...244") %>% 
  rmv_stbl()
`data_2022-08-07-0941_walkup_aci_clean.xlsx` <- `data_2022-08-07-0941_walkup_aci_clean.xlsx` %>% 
  rename("Tree_Identifier" = "Tree Identifier") %>% 
  rename("Data_QC" = "...244") %>% 
  rmv_stbl()
`data_2022-08-11-0904_walkup_aci_clean.xlsx` <- `data_2022-08-11-0904_walkup_aci_clean.xlsx`%>% 
  rename("Tree_Identifier" = "Tree Identifier") %>% 
  rmv_stbl()
`data_2022-08-11-0949_walkup_aci_clean.xlsx` <- `data_2022-08-11-0949_walkup_aci_clean.xlsx`%>% 
  rename("Tree_Identifier" = "Tree Identifier") %>% 
  rmv_stbl()
`data_2022-08-11-1050_walkup_aci_clean.xlsx` <- `data_2022-08-11-1050_walkup_aci_clean.xlsx`%>% 
  rename("Tree_Identifier" = "Tree Identifier") %>% 
  rename("Data_QC" = "...248") %>% 
  rmv_stbl()
`data_2022-08-12-0928_walkup_aci_clean.xlsx` <- `data_2022-08-12-0928_walkup_aci_clean.xlsx`%>% 
  rename("Tree_Identifier" = "Tree Identifier") %>% 
  rename("Data_QC" = "...248") %>% 
  rmv_stbl()
`data_2022-08-13-1007_walkup_aci_clean.xlsx` <- `data_2022-08-13-1007_walkup_aci_clean.xlsx`%>% 
  rename("Tree_Identifier" = "Tree Identifier") %>% 
  rename("Data_QC" = "...248") %>% 
  rmv_stbl()
`data_2022-08-14-0906_walkup_aci_clean.xlsx` <- `data_2022-08-14-0906_walkup_aci_clean.xlsx`%>% 
  rename("Tree_Identifier" = "Tree Identifier") %>% 
  rename("Data_QC" = "...248") %>% 
  rmv_stbl()
`data_2022-08-08-0854_walkup_aci_clean.xlsx` <- rmv_stbl(`data_2022-08-08-0854_walkup_aci_clean.xlsx`)
`data_2022-08-09-0911_walkup_aci_clean.xlsx` <- rmv_stbl(`data_2022-08-09-0911_walkup_aci_clean.xlsx`)
`data_2022-08-10-0854_walkup_aci_clean.xlsx` <- rmv_stbl(`data_2022-08-10-0854_walkup_aci_clean.xlsx`)
`data_2022-08-16-walkup_aci_clean.xlsx` <- rmv_stbl(`data_2022-08-16-walkup_aci_clean.xlsx`)



## Fix data columns into the same kind of data
`data_2022-08-06-1020_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-06-1020_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-07-0941_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-07-0941_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-11-0904_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-11-0904_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-11-0949_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-11-0949_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-11-1050_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-11-1050_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-12-0928_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-12-0928_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-13-1007_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-13-1007_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-14-0906_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-14-0906_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-08-0854_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-08-0854_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-09-0911_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-09-0911_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-10-0854_walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-10-0854_walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-16-walkup_aci_clean.xlsx`$"Leaf_position" <- as.numeric(`data_2022-08-16-walkup_aci_clean.xlsx`$"Leaf_position")
`data_2022-08-16-walkup_aci_clean.xlsx`$"hhmmss...5" <- as.character(`data_2022-08-16-walkup_aci_clean.xlsx`$"hhmmss...5")
`data_2022-08-16-walkup_aci_clean.xlsx`$"Leaf_number" <- as.character(`data_2022-08-16-walkup_aci_clean.xlsx`$"Leaf_number")
`data_2022-08-16-walkup_aci_clean.xlsx`$"hhmmss...111" <- as.character(`data_2022-08-16-walkup_aci_clean.xlsx`$"hhmmss...111")
`data_2022-08-16-walkup_aci_clean.xlsx`$"State" <- as.character(`data_2022-08-16-walkup_aci_clean.xlsx`$"State")
`data_2022-08-16-walkup_aci_clean.xlsx`$"GPIO" <- as.character(`data_2022-08-16-walkup_aci_clean.xlsx`$"GPIO")



## Merge the data together
all_aci.df <- bind_rows(`data_2022-08-06-1020_walkup_aci_clean.xlsx`,
                        `data_2022-08-07-0941_walkup_aci_clean.xlsx`,
                        `data_2022-08-11-0904_walkup_aci_clean.xlsx`,
                        `data_2022-08-11-0949_walkup_aci_clean.xlsx`,
                        `data_2022-08-11-1050_walkup_aci_clean.xlsx`,
                        `data_2022-08-12-0928_walkup_aci_clean.xlsx`,
                        `data_2022-08-13-1007_walkup_aci_clean.xlsx`,
                        `data_2022-08-14-0906_walkup_aci_clean.xlsx`,
                        `data_2022-08-08-0854_walkup_aci_clean.xlsx`,
                        `data_2022-08-09-0911_walkup_aci_clean.xlsx`,
                        `data_2022-08-10-0854_walkup_aci_clean.xlsx`,
                        `data_2022-08-16-walkup_aci_clean.xlsx`)




## Convert important columns to numeric
complete_sp <- all_aci.df %>%
    mutate(Ci = parse_number(Ci)) %>%
    mutate(A = parse_number(A)) %>% 
    mutate(Ca = parse_number(Ca)) %>% 
    mutate(E = parse_number(E)) %>%
    mutate(Qin = parse_number(Qin)) %>%
    mutate(Tleaf = parse_number(Tleaf))
    # make a new column for the tree_id using the K67 ID
complete_sp <- mutate(complete_sp, tree_id = case_match(complete_sp$Tree_Identifier,
                                         c('Maca1', 'maca1') ~ 'k6709',
                                         '1' ~ 'k6707',
                                         'Tree3' ~ 'k6708',
                                         '4' ~ 'k6706', 
                                         'Tree5' ~ 'k6704',
                                         'Tree6' ~ 'k6705',
                                         'Mela7' ~ 'k6716',
                                         c('tree8', 'Tree8') ~ 'k6715',
                                         'tree9' ~ 'k6717',
                                         'Tree10' ~ 'k6713',
                                         'tree11' ~ 'k6702',
                                         'tree12' ~ 'k6714',
                                         'tree22' ~ 'k6718'))
unique(complete_sp$tree_id)
# Once the bad data is filtered out, there will be 13 trees


# Save Complete Data File -------------------------------------------------

## Save the assembled file as a .csv
write.csv(x = complete_sp, 
          file = here("1_Raw_data/raw_combined_aci_data.csv"), 
          row.names = FALSE)





# This file is not the final version we used because we had to go back and create a new 'unique_id' column with proper leaf numbers. This was done manually, but uses the above csv as its base. Most of the excluded points are already filtered out.



# Clean and Filter the data -----------------------------------------------

## Load the data file with the fixed names
cmplt.rm_out <- read.csv(here("3_Clean_data/clean_aci_with_uniquecode.csv"),
                         sep = ",", header = TRUE,
                         fileEncoding="latin1") %>% 
    # Filter out excluded data
    filter(Data_QC == "OK") %>% 
    # Filter out outliers (impossible values)
    filter(Ci > -5 & A < 40 & A > -1) %>% 
    # Fix curve method values and column name
    mutate(Data_point = case_match(Data_point, 
                                   'Before_DAT' ~ "DAT", 'Traditional' ~ "SS"),
           curv_meth = Data_point)

write.csv(cmplt.rm_out, 
          file = here("3_Clean_data/clean_aci_noOutliers.csv"),
          row.names = FALSE)



