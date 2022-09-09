# Purpose -----------------------------------------------------------------

## This script is for cleaning and assembling the ACi curve data collected at the
## Tapajos National Forest in Aug 2022 into one file for further analysis



# Load Packages -----------------------------------------------------------

library(tidyverse)
library(googledrive)
library(readxl)


####### Set Working Directory
setwd("C:/Users/emmel/Desktop/DAT_proj")

# Load from Google Drive (in progress) ------------------------------------
# 
# ## Code to import one file
# file_name <- read_excel("2022-08-14-0906_walkup_aci_clean.xlsx", sheet = "Measurements", skip = 14)
# View(file_name)
# 
# #_____________________________________________________________________________
# ## First, authorize R to access Google Drive
# drive_auth(email = "eb0067@mix.wvu.edu")
# drive_user()
# 
# #Check the files in your drive
# drive_find(n_max = 30)
# 
# # Make a list of all the Files you want to upload
# aci_files_cln <- drive_ls(path =
#                                # "Albert_lab_ecohydromics/Gas_exchange_Aug_2022/Clean/LI-6800",
#                             as_id("https://drive.google.com/drive/u/2/folders/1nSKLtjd69BJWFlSwEl-r6d8mqK3mF-9c"),
#                             pattern = "walkup_aci_clean")
# print(aci_files_cln)
# 
# #Make a FOR loop to download all the files
# for (i in 1:length(aci_files_cln$name)) {
#   links <- drive_link(aci_files_cln$name)
#   files <- drive_download(links[i], overwrite = TRUE)
#   data <- read_excel(files[[i]], sheet = "Measurements", skip = 14)
#   #names <- filenames[[i]]
#   assign(paste0("data_", aci_files_cln$name[i]), data)
# }






# Load Data from Directory ------------------------------------------------
aci_files_nm <- list.files(path = "C:/Users/emmel/Desktop/DAT_proj",
                           pattern = "walkup_aci_clean")
print(aci_files_nm)

for (i in 1:length(aci_files_nm)) {
  files1 <- read_xlsx(path = paste0("C:/Users/emmel/Desktop/DAT_proj/", aci_files_nm[i]),
                      sheet = "Measurements", skip = 14)
  sliced <- slice(files1, -(1))
  assign(paste0("data_", aci_files_nm[i]), sliced)                     
}





# Assemble into one data file ---------------------------------------------

## Make sure all the data have the same columns (in this case, the 6th and 7th are missing
## some columns)
Reduce(setdiff, list(colnames(`data_2022-08-06-1020_walkup_aci_clean.xlsx`), 
        colnames(`data_2022-08-07-0941_walkup_aci_clean.xlsx`),
        colnames(`data_2022-08-11-0904_walkup_aci_clean.xlsx`),
        colnames(`data_2022-08-11-0949_walkup_aci_clean.xlsx`),
        colnames(`data_2022-08-11-1050_walkup_aci_clean.xlsx`),
        colnames(`data_2022-08-12-0928_walkup_aci_clean.xlsx`),
        colnames(`data_2022-08-13-1007_walkup_aci_clean.xlsx`),
        colnames(`data_2022-08-14-0906_walkup_aci_clean.xlsx`),
        colnames(`data_2022-08-08-0854_walkup_aci_clean.xlsx`),
        colnames(`data_2022-08-09-0911_walkup_aci_clean.xlsx`),
        colnames(`data_2022-08-10-0854_walkup_aci_clean.xlsx`),
        colnames(`data_2022-08-16-walkup_aci_clean.xlsx`))) #finds the unique columns


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



## Clean the data by filtering out excluded points and After DAT points
all_aci_cln <- all_aci.df %>% 
  filter(Data_QC == "OK") %>% 
  filter(Data_point != "After_DAT")
View(all_aci_cln)


## Convert important columns to numeric
all_aci_cln_num <- all_aci_cln %>%
  mutate(Ci = parse_number(Ci)) %>%
  mutate(A = parse_number(A)) %>% 
  mutate(Ca = parse_number(Ca)) %>% 
  mutate(E = parse_number(E)) %>%
  mutate(Qin = parse_number(Qin)) %>%
  mutate(Tleaf = parse_number(Tleaf))




# Adding in Species -------------------------------------------------------

### adding a column for a four-letter species code and a column for species name
complete_sp <- all_aci_cln_num %>% 
  mutate(fourlettercode = Tree_Identifier, SciName = Tree_Identifier)

complete_sp$fourlettercode <- recode(complete_sp$fourlettercode,
                                     'Maca1' = 'MAEL',
                                     'Tree3' = 'CHTU',
                                     'tree8' = 'COSP',
                                     'Tree8' = 'COSP',
                                     'tree9' = 'APCO',
                                     'Tree10' = 'VICA',
                                     'tree11' = 'COST',
                                     'tree12' = 'UNKN',
                                     'tree22' = 'ABMA',
                                     '1' = 'TACH',
                                     '4' = 'ESSP',
                                     'Tree5' = 'ABMA',
                                     'Tree6' = 'PRAP',
                                     'Mela7' = 'MELA',
                                     'maca1' = 'MAEL')
complete_sp$SciName <- recode(complete_sp$SciName,
                              'Maca1' = 'Manilkara elata',
                              'Tree3' = 'Chimaris turbinata',
                              'tree8' = 'Coussarea sp',
                              'Tree8' = 'Coussarea sp',
                              'tree9' = 'Aparisthmium cordatum',
                              'Tree10' = 'Vismia cayennensis',
                              'tree11' = 'Couratari stellata',
                              'tree12' = 'Unknown sp',
                              'tree22' = 'Abarema mataybifolia',
                              '1' = 'Tachigali chrysophylla',
                              '4' = 'Eschweilera sp.',
                              'Tree5' = 'Abarema mataybifolia',
                              'Tree6' = 'Protium apiculatum',
                              'Mela7' = 'Unknown sp.',
                              'maca1' = 'Manilkara elata')
unique(complete_sp$fourlettercode)



# Save Complete Data File -------------------------------------------------
## Save the assembled file as a .csv
write.csv(x = complete_sp, file = "clean_aci_data_one_file.csv")




###which(duplicated(names(`data_2022-08-16-walkup_aci_clean.xlsx`)))



## Lessons learned:
### Add a Data_QC column as a user constant and just populate it afterward- this avoids
### Us accidentally putting the column name in the wrong row

### Absolutely DO NOT mess with user constants/stability criteria in between files
### At the very least, just double-check that capitalization rules are the same. We 
### should standardize this just in case we lose our configuration

### It's better to just remove the stability criteria columns rather than trying to deal