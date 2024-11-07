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

allfiles = allfiles[72:2729] #delete practice/test scans

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
metadata <- read_excel("./metadata_updated.xlsx")
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
traits_metadata<- traits[,c(1,4,6,7,8,9)]
traits_metadata<- unique(traits_metadata)

allspectra[ , 4:1001] <- lapply(allspectra[ , 4:1001], as.numeric)
allspectra <- allspectra[ , c(1,4:1001)]
allspectra_merging <- allspectra %>%
  group_by(barcode) %>%  
  summarise(across(2:998, \(x) mean(x, na.rm = TRUE)))

mergedspectra <- merge(allspectra_merging, traits_metadata, by.x='barcode', by.y="id")
mergedspectra <- mergedspectra[,c(1,999:1003,2:998)]

write.csv(mergedspectra, "./pftc7_spectra.csv")


###### Plot data ########

spectra <- mergedspectra %>%
  pivot_longer(
    cols = 7:ncol(mergedspectra),           # Select spectral columns (columns 6 to end)
    names_to = "lambda",          # New column for wavelengths
    values_to = "reflectance"  # New column for spectral values
  )


#### Plot all species
all_per_species <- spectra %>% group_by(species, lambda) %>% 
  summarize(mean_r = mean(reflectance),
            sd_r = sd(reflectance),
            n = n(),
            se_r = sd_r / sqrt(n)) %>% ungroup()
all_per_species <- data.table(all_per_species)
all_per_species$species <- gsub("_", " ", all_per_species$species)
all_per_species$species <- gsub("cf", "", all_per_species$species)
all_per_species$species <- gsub("c f", "", all_per_species$species)

all_per_species[, c("genus") := tstrsplit(species, " ", keep = c(1))]
all_per_species$lambda <- as.numeric(all_per_species$lambda)
all_per_species <- all_per_species[!(genus == "thin")] #remove unknown species
all_per_species <- all_per_species[!(genus == "purple")] #remove unknown species
all_per_species <- all_per_species[!(genus == "")] #remove unknown species



all_per_genus <- all_per_species %>% group_by(genus, lambda) %>% 
  summarize(mean_r = mean(mean_r),
            sd_r = sd(mean_r),
            n = n(),
            se_r = sd_r / sqrt(n)) %>% ungroup()
all_per_genus <- all_per_genus %>% filter(lambda > 400)
all_per_genus <- all_per_genus %>% filter(lambda < 2500)

#plot all species
allspeciesplot <- ggplot(data = all_per_genus,
                         aes(x = lambda, y = mean_r, color = genus)) +
  #  geom_ribbon(aes(ymin = 100*(mean_r-2*se_r), ymax = 100*(mean_r+2*se_r), fill = genus), alpha = 0.3, color = NA) +
  ylab("% Reflectance") + 
  xlab("Wavelength (nm)") +
  geom_line(aes(y = 100*mean_r), linewidth=1, alpha=.8) +
  scale_colour_manual(values=unname(polychrome()))+
  theme(panel.background = element_blank(),  # Removes the background color
        panel.grid.major = element_blank(),  # Removes major grid lines
        panel.grid.minor = element_blank(),  # Removes minor grid lines
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14),  # Change axis labels text size
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = "none")
allspeciesplot

table(all_per_species$species)
ggsave("per_all_species.pdf", allspeciesplot, units = "mm", width = 300, height = 150)

#plot focal species
focals <- c( 'helichrysum pallidum', "senecio glaberrimus", "helichrysum ecklonis")
sspectra <- spectra %>% filter(species %in% focals)

per_species_focal <- sspectra %>% group_by(species, lambda) %>% 
  summarize(mean_r = mean(reflectance),
            sd_r = sd(reflectance),
            n = n(),
            se_r = sd_r / sqrt(n)) %>% ungroup()

per_species_focal$lambda <- as.numeric(per_species_focal$lambda)
per_species_focal <- per_species_focal %>% filter(lambda > 400)
per_species_focal <- per_species_focal %>% filter(lambda < 2500)

speciesplot <- ggplot(data = per_species_focal,
                      aes(x = lambda, color = species, fill = species)) +
  geom_ribbon(aes(ymin = 100*(mean_r-2*se_r), ymax = 100*(mean_r+2*se_r), fill = species), alpha = 0.3, color = NA) +
  ylab("% Reflectance") + 
  xlab("Wavelength (nm)") +
  geom_line(aes(y = 100*mean_r),linewidth=1,) + 
  scale_colour_manual(values=c("#1446a0","#db3069","#f5d547"))+
  scale_fill_manual(values=c("#1446a0","#db3069","#f5d547"))+
  theme(panel.background = element_blank(),  # Removes the background color
        panel.grid.major = element_blank(),  # Removes major grid lines
        panel.grid.minor = element_blank(),  # Removes minor grid lines
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14),  # Change axis labels text size
        axis.text = element_text(size = 12),
        axis.text.y = element_blank(),
        legend.position = "none")
speciesplot
ggsave("per_focal_species.pdf", speciesplot, units = "mm", width = 200, height = 100)

comp <- ggarrange(allspeciesplot, speciesplot,
                  labels = c("A", "B"),
                  ncol = 1, nrow = 2)
comp
ggsave("focal_all.pdf", comp, units = "mm", width = 175, height = 240)
