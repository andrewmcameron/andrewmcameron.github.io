library(tidyverse)
library(kableExtra)



##  Read in NEON data (forests)

# Create vector with NEON folder names. 
folders <- list.files(paste0(getwd(), "/NEON"))

# Index of forest data folders within NEON folder
folder_indices <- 
  which(grepl("SERC", folders) | 
          grepl("SCBI", folders) |
          grepl("BLAN", folders))


# Create function to more quickly join lat/long to count data for each study area-year
neon_join <- function(x, y) {
  result <- x %>%
    left_join(y %>% 
                dplyr::select(plotID, pointID, decimalLatitude, decimalLongitude, 
                              nlcdClass, elevation, observedHabitat), 
              by = c("plotID", "pointID")
    )
  
  return(result)
}


# ------#Use `neon_join()` to join each study area-year's two csv files----------
joined_dfs <- list()   # Empty list to store the 9 joined dfs

for (i in seq_along(folder_indices)) {
  
  files <- list.files(paste0(getwd(), 
                             "/NEON/", 
                             paste(folders[folder_indices[i]]) )  )
  
  counts <-  read_csv(
    paste0(getwd(), "/NEON/", 
           paste(folders[folder_indices[i]]), "/",
           paste(files[which(grepl("countdata", files))] ) )  )
  
  points <-  read_csv(
    paste0(getwd(), "/NEON/", 
           paste(folders[folder_indices[i]]), "/",
           paste(files[which(grepl("perpoint", files))]) )  )  
  
  
  joined_df <- neon_join(counts, points)
  joined_df <- distinct(joined_df)  # to remove any unwanted results of many-to-many relationships
  
  joined_dfs[[i]] <- joined_df
  
}


# Extract individual dfs from list, then create single data frame for each study site containing all three years' data
bland2019 <- joined_dfs[[1]]
bland2020 <- joined_dfs[[2]]
bland2021 <- joined_dfs[[3]]
smith2019 <- joined_dfs[[4]]
smith2020 <- joined_dfs[[5]]
smith2021<- joined_dfs[[6]]
serc2019 <- joined_dfs[[7]]
serc2020 <- joined_dfs[[8]]
serc2021 <- joined_dfs[[9]]

bland_all <- rbind(bland2019, bland2020, bland2021)
serc_all <- rbind(serc2019, serc2020, serc2021)
smith_all <- rbind(smith2019, smith2020, smith2021)

# Create a 'plot-year' field to more readily sort and group observations
bland_all$plotYear <- paste0(bland_all$plotID, "-", year(bland_all$startDate))
serc_all$plotYear <- paste0(serc_all$plotID, "-", year(serc_all$startDate))
smith_all$plotYear <- paste0(smith_all$plotID, "-", year(smith_all$startDate))

# De-clutter environment a bit
objects <- ls(envir = globalenv())
rm(list = objects[grepl("2019|202|2021|joined|file", objects)])

# Vector of site dfs for use in for loops
all_sites <- c("bland_all", "smith_all", "serc_all")


# -- Same as previous chunk, for GRASSLANDS
folders_g <- list.files(paste0(getwd(), "/NEON"))

folder_indices <- 
  which(grepl("KONA", folders) | 
          grepl("KONZ", folders) |
          grepl("UKFS", folders))


# empty list to store joined dfs
joined_dfs_g <- list()

for (i in seq_along(folder_indices)) {
  
  files <- list.files(paste0(getwd(), 
                             "/NEON/", 
                             paste(folders[folder_indices[i]]) )  )
  
  counts <-  read_csv(
    paste0(getwd(), "/NEON/", 
           paste(folders[folder_indices[i]]), "/",
           paste(files[which(grepl("countdata", files))] ) )  )
  
  points <-  read_csv(
    paste0(getwd(), "/NEON/", 
           paste(folders[folder_indices[i]]), "/",
           paste(files[which(grepl("perpoint", files))]) )  )  
  
  
  joined_df <- neon_join(counts, points)
  joined_df <- distinct(joined_df)  # to remove any unwanted results of many-to-many relationships
  
  joined_dfs_g[[i]] <- joined_df
  
}


# Extract individual dfs from list, then create single data frame for each study site containing all three years of data
kona2019 <- joined_dfs_g[[1]]
kona2020 <- joined_dfs_g[[2]]
kona2021 <- joined_dfs_g[[3]]
konz2019 <- joined_dfs_g[[4]]
konz2020 <- joined_dfs_g[[5]]
konz2021 <- joined_dfs_g[[6]]
kans2019 <- joined_dfs_g[[7]]
kans2020 <- joined_dfs_g[[8]]
kans2021 <- joined_dfs_g[[9]]

kona_all <- rbind(kona2019, kona2020, kona2021)
konz_all <- rbind(konz2019, konz2020, konz2021)
kans_all <- rbind(kans2019, kans2020, kans2021)

# Create a 'plot-year' field to more readily sort and group observations
kona_all$plotYear <- paste0(kona_all$plotID, "-", year(kona_all$startDate))
konz_all$plotYear <- paste0(konz_all$plotID, "-", year(konz_all$startDate))
kans_all$plotYear <- paste0(kans_all$plotID, "-", year(kans_all$startDate))

# De-clutter environment a bit
objects <- ls(envir = globalenv())
rm(list = objects[grepl("2019|202|2021|joined|file", objects)])

# Vector of site dfs for use in for loops
all_sites_g <- c("kona_all", "kans_all", "konz_all")




  
  ### Richness & Diversity by Plot-Year
  
### Species richness 
richness <- list()

for (i in seq_along(all_sites)) {
  rich <- get(all_sites[i]) %>%
    group_by(plotYear) %>%
    filter(taxonRank == "species") %>%      
    summarize(richness = length(unique(taxonID)),
              elevation = max(elevation))
  
  richness[[i]] <- rich
}

# Create dfs for eventual use in analysis/modeling 
bland_final <- richness[[1]]
smith_final <- richness[[2]]
serc_final <- richness[[3]]

### Shannon Diversity Index
# Create abundance matrices and calculate Shannon vals
abun_mats <- list()

for (i in seq_along(all_sites)){
  matrix <- 
    get(all_sites[i]) %>%
    dplyr::select(plotYear, taxonID, clusterSize) %>%
    group_by(plotYear, taxonID) %>%
    summarize(abundance = sum(clusterSize)) %>%
    pivot_wider(names_from = taxonID, values_from = abundance) %>%
    column_to_rownames(var = "plotYear")
  
  # Replace NAs with 0
  df <- as.data.frame(matrix)
  df[is.na(df)] <- 0
  matrix <-as.matrix(df)
  
  abun_mats[[i]] <- matrix
}


shan_list <- lapply(abun_mats, vegan::diversity)
shan_list <- lapply(shan_list, as.data.frame)
shan_list <- lapply(shan_list, function(df) {
  rownames_to_column(df, var = "rowname")
})

# Join shan vals to 'final' dfs
bland_final <- bland_final %>%
  left_join(shan_list[[1]], 
            join_by(plotYear == rowname))

smith_final <- smith_final %>%
  left_join(shan_list[[2]], 
            join_by(plotYear == rowname))

serc_final <- serc_final %>%
  left_join(shan_list[[3]], 
            join_by(plotYear == rowname))


# Rename col with shan vals
colnames(bland_final)[4] <- "shan"
colnames(smith_final)[4] <- "shan"
colnames(serc_final)[4] <- "shan"


# Join lat/long data for mapping and spatial analysis
bland_final <- bland_final %>%
  left_join(
    bland_all %>% 
      dplyr::select(plotYear, decimalLatitude, decimalLongitude) %>%
      distinct(),
    by = "plotYear")

smith_final <- smith_final %>%
  left_join(
    smith_all %>% 
      dplyr::select(plotYear, decimalLatitude, decimalLongitude) %>%
      distinct(),
    by = "plotYear") 

serc_final <- serc_final %>%
  left_join(
    serc_all %>% 
      dplyr::select(plotYear, decimalLatitude, decimalLongitude) %>%
      distinct(),
    by = "plotYear") 



### Species richness 
richness_g <- list()

for (i in seq_along(all_sites)) {
  rich <- get(all_sites_g[i]) %>%
    group_by(plotYear) %>%
    filter(taxonRank == "species") %>%      
    summarize(richness = length(unique(taxonID)),
              elevation = max(elevation))
  
  richness_g[[i]] <- rich
}

# Create dfs for eventual use in analysis/modeling 
kona_final <- richness_g[[1]]
kans_final <- richness_g[[2]]
konz_final <- richness_g[[3]]

### Shannon Diversity Index
# Create abundance matrices and calculate SDI vals
abun_mats_g <- list()

for (i in seq_along(all_sites_g)){
  matrix <- 
    get(all_sites_g[i]) %>%
    dplyr::select(plotYear, taxonID, clusterSize) %>%
    group_by(plotYear, taxonID) %>%
    summarize(abundance = sum(clusterSize)) %>%
    pivot_wider(names_from = taxonID, values_from = abundance) %>%
    column_to_rownames(var = "plotYear")
  
  # Replace NAs with 0
  df <- as.data.frame(matrix)
  df[is.na(df)] <- 0
  matrix <-as.matrix(df)
  
  abun_mats_g[[i]] <- matrix
}

shan_list_g <- lapply(abun_mats_g, vegan::diversity)
shan_list_g <- lapply(shan_list_g, as.data.frame)
shan_list_g <- lapply(shan_list_g, function(df) {
  rownames_to_column(df, var = "rowname")
})

# Join shan vals to 'final' dfs
kona_final <- kona_final %>%
  left_join(shan_list_g[[1]], 
            join_by(plotYear == rowname))

kans_final <- kans_final %>%
  left_join(shan_list_g[[2]], 
            join_by(plotYear == rowname))

konz_final <- konz_final %>%
  left_join(shan_list_g[[3]], 
            join_by(plotYear == rowname))


# Rename col with shan vals
colnames(kona_final)[4] <- "shan"
colnames(kans_final)[4] <- "shan"
colnames(konz_final)[4] <- "shan"


# Join lat/long data for mapping and spatial analysis
kona_final <- kona_final %>%
  left_join(
    kona_all %>% 
      dplyr::select(plotYear, decimalLatitude, decimalLongitude) %>%
      distinct(),
    by = "plotYear")

kans_final <- kans_final %>%
  left_join(
    kans_all %>% 
      dplyr::select(plotYear, decimalLatitude, decimalLongitude) %>%
      distinct(),
    by = "plotYear") 

konz_final <- konz_final %>%
  left_join(
    konz_all %>% 
      dplyr::select(plotYear, decimalLatitude, decimalLongitude) %>%
      distinct(),
    by = "plotYear") 

