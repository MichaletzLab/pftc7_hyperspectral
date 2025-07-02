### Cleaning hyperspectral measurements 
library(data.table)
library(ggplot2)
library(spectrolab)
library(stringr)
library(tidyverse) 
library(MetBrewer)
library(readxl)
library(rcartocolor)
library(ggpubr)
library(pals)

## find the files to read
allfiles = list.files("./",
                      full.names = T,
                      recursive = T,
                      pattern = 'sig')

allfiles = allfiles[72:2729] #delete practice/test scans pror to measurements commencing

# loop through, read file, correct sensor overlap, return as data table
allspectra = rbindlist(lapply(allfiles, function(x){
  y = read_spectra(x)
  y2 = match_sensors(y, c(990,1900))
  return(as.data.table(y2))
}))


#replace periods with underscore for consistent delimiter in sample name
allspectra$sample_name <- gsub(".", "_", allspectra$sample_name, fixed = TRUE)
allspectra$sample_name <- gsub("-", "_", allspectra$sample_name, fixed = TRUE)

#pull date, barcodes and scan numbers from sample name
allspectra[, c("date","barcode", "scan_num") := tstrsplit(sample_name, "_", keep = c(1,3,4))]

#remove white references, test scans
allspectra$barcode <- toupper(allspectra$barcode)
# remove white references
allspectra <- allspectra %>% filter(barcode != "WR")
# remove practice test scans
allspectra <- allspectra %>% filter(barcode != "TEST")
allspectra <- allspectra %>% filter(barcode != "SIG")
allspectra <- allspectra %>% filter(!grepl('BARCODE', barcode))

# read in metadata xlsx
metadata <- read_excel("./metadata/metadata_updated.xlsx")
View(metadata)
metadata <- data.table(metadata)
# rename columns
lookup <- c(barcode = "sample_id", scan_number = "scan number", to_delete='delete?', to_rename='rename?')
metadata <- rename(metadata, all_of(lookup))
metadata$to_delete <- toupper(metadata$to_delete)
# delete rows with no barcodes
metadata <- metadata %>% filter(!is.na(barcode))
# remove scans we decided to include in the analysis
metadata <- metadata %>% filter(to_delete != "N")
# remove rows where no scan is specified
metadata <- metadata %>% filter(!is.na(scan_number))
# make the scan number an integer
allspectra$scan_num <- lapply(allspectra$scan_num, strtoi)
# remove test rows where no scan number was added
metadata <- metadata %>% filter(!is.na(scan_number)) 
# give each scan number their own row
metadata[, c("scan_0","scan_1", "scan_2", "scan_3") := tstrsplit(metadata$scan_number, ",")]
metadata$scan_number <- NULL
metadata <- pivot_longer(metadata, c('scan_0', 'scan_1', 'scan_2', 'scan_3'), names_to=NULL, values_to='scan_number', values_drop_na=TRUE)
metadata$scan_number <- lapply(metadata$scan_number, strtoi)
# filter to the barcodes to update
rename <- metadata %>% filter(!is.na(to_rename))
# update barcode for swapped scans
allspectra[(scan_num %in% rename$scan_number) & (barcode %in% rename$barcode), barcode := rename$to_rename]


# only delete scans from metadata that are for sure Y to delete
metadata <- metadata %>% filter(to_delete == 'Y')

# finally, remove rows from allspectra that are marked to delete
# cast types
allspectra$scan_num <- as.numeric(allspectra$scan_num)
metadata$scan_number <- as.numeric(metadata$scan_number)
metadata$date <- NULL
allspectra$date <- NULL

# join metadata and spectra dataframes
allspectra <- merge(allspectra, metadata, by.x=c('barcode', 'scan_num'), by.y=c('barcode', 'scan_number'), all.x=T)
# filter out barcodes + scans that were present in metadata
allspectra <- allspectra %>% filter(is.na(allspectra$to_delete))


# read in trait dfs
traits <- read.csv("PFTC7_SA_clean_traits_2023.csv")
traits_metadata<- traits[,c(1,2,4,6,7,8,9)]
traits_metadata<- unique(traits_metadata)

allspectra[ , 4:1001] <- lapply(allspectra[ , 4:1001], as.numeric)
allspectra <- allspectra[ , c(1,4:1001)]
allspectra_merging <- allspectra %>%
  group_by(barcode) %>%  
  summarise(across(2:998, \(x) mean(x, na.rm = TRUE)))

mergedspectra <- merge(allspectra_merging, traits_metadata, by.x='barcode', by.y="id")
mergedspectra <- mergedspectra[,c(1,999:1004,2:998)]

#capitalize first letter of genus for consistency
mergedspectra <- mergedspectra %>%
  mutate(
    species = str_to_lower(species),  # make sure everything is lowercase first
    genus = str_to_title(word(species, 1)),  # capitalize first letter of genus
    epithet = word(species, 2),             # keep species epithet lowercase
    species = paste(genus, epithet)         # recombine
  ) %>%
  select(-genus, -epithet)

mergedspectra <- mergedspectra %>%
  mutate(site_id = dense_rank(elevation_m_asl))%>%
  rename(id = barcode)%>%
  select(id,date, aspect, site_id, elevation_m_asl, plot_id,plant_id, everything())

write.csv(mergedspectra, "./pftc7_spectra.csv")

