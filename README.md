# Pyrenean Institute of Ecology: Ecoacoustic Protocols  
![Title image](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/d1275870-6d41-46f6-8c4b-a1dd37745752)
### By [Jack A. Greenhalgh](https://www.jack-greenhalgh.com/), September 2023.     

### Contents 

- [Part 1: Ecoacoustic equipment and survey design](#ecoacoustic-equipment-and-survey-design)
- [Part 2: Acoustic indices](#acoustic-indices)
- [Part 3: Tapestry plots: Visualising long-term passive acoustic monitoring data](#tapestry-plots-visualising-long-term-passive-acoustic-monitoring-data)
- [Part 4: False-colour spectrograms: Visualising 24 hrs of passive acoustic monitoring data](#false-colour-spectrograms-plotting-24-hrs-of-acoustic-data-with-noise-reduction)

### Ecoacoustic equipment and survey design

#### Acoustic recorders

Here is a list of the most widely used acoustic recorders for long-term soundscape monitoring in a variety of environments. 

| Acoustic recorder  | Realm(s)           | Cost ($) | Link                                                                         |
|--------------------|--------------------|----------|------------------------------------------------------------------------------|
| AudioMoth          | Terrestrial        | 89.99    | https://www.openacousticdevices.info/audiomoth                               |
| HydroMoth          | Marine, Freshwater | 135.00   | https://groupgets.com/manufacturers/open-acoustic-devices/products/hydromoth |
| Song Meter Mini    | Terrestrial        | 576.00   | https://www.veldshop.nl/en/song-meter-mini.html                              |
| Song Meter SM4 BAT | Terrestrial        | 999.00   | https://www.wildlifeacoustics.com/products/song-meter-sm4bat                 |

### Types of acoustic survey design

Deciding the best way to deploy your acoustic recorders in the environment that you plan to survey is very important and worthy of careful consideration. The optimal way to deploy your recorders will depend on several factors, including, how many recorders you have, how frequently they need to be serviced, and the research questions that you are interested in answering. 

#### Line transect 
Line transects involve collecting data along predetermined linear paths (transects) within the study area. This is especially useful for studying gradients or patterns over space.

![survey design](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/2a7607ff-7031-4545-bad9-66df46b4a7d1)

#### Random sampling 

Simple random sampling:
Involves randomly selecting sample sites or individuals from the entire study area. It's useful for ensuring each unit has an equal chance of being sampled.

Stratified random sampling: 
Divides the study area into subgroups (strata) based on certain characteristics (e.g., habitat type) and then randomly samples within each stratum. This ensures representation of different ecosystem components.

 #### Recording schedule and parameters   

Something else to carefully consider is your recording schedule and parameters. Acoustic recorders can be pre-programmed to record continuously, or selectively (e.g., 1 min in every 10). There are pros and cons to each approach, but the benefits of selective recording schedules, such as longer intervals between servicing the recorders, are clear for long-term soundscape monitoring. 

![recording schedule](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/9df09336-14ef-427b-a98f-cdb8d88220c0)
AudioMoth configuration app showing: (a) recording parameters, and b) recording schedule.

You can also modify the sample rate (kHz), gain (dB), and other more detailed parameters such as trigger type and filtering. These settings will largely be dictated by your research question. For example, high sample rates (>192 kHz) and selective trigger types are ideal for bat surveys, whereas mid-range sample rates (~48 kHz) and no triggers are ideal for soundscape monitoring.        

###  Acoustic indices 

#### What are acoustic indcies?

Acoustic indcies are mathametical functions that consider spectral and temporal information obtained from audio recordings. They can be used to  analyse large audio datasets and monitor biodiversity without the need to determine species identity. 

Many different types of acoustic indcies exist, however, here is a list of the most commonly used indices. 

| Acoustic index                                | Description                                                                                                                                                                                                                                                                                                               | Developed for                                                                                      | Reference                                                 |
|-----------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------|-----------------------------------------------------------|
| Spectral entropy (Hf)                         | “Measures the evenness of the amplitude envelope over the time units… by dividing the Shannon index by its maximum” (Sueur et al., 2008)                                                                                                                                                                                  | Terrestrial fauna (in a coastal forest, Tanzania)                                                  | Sueur et al. (2008)                                       |
| Temporal entropy (Ht)                         | “A mean spectrum s(f) is first computed using a Short Time Fourier Transform (STFT) based on a nonoverlapping sliding function window of sample width τ. This mean spectrum s(f) is similarly transformed into a probability mass function S(f) of length N used to compute the spectral entropy Ht” (Sueur et al., 2008) | Terrestrial fauna (in a coastal forest, Tanzania)                                                  | Sueur et al. (2008)                                       |
| Acoustic entropy (H)                          | A function of Ht and Hf                                                                                                                                                                                                                                                                                                   | Terrestrial fauna (in a coastal forest, Tanzania)                                                  | Sueur et al. (2008)                                       |
| Acoustic richness (AR)                        | A ranked index based on the temporal entropy and amplitude of a signal                                                                                                                                                                                                                                                    | Birds (in a temperate woodland, France)                                                            | Depraetere et al. (2012)                                  |
| Acoustic evenness index (AEI)                 | “Calculated by dividing the spectrogram into bins (default 10) and taking the proportion of the signals in each bin above a threshold (default −50 dBFS). The AEI is the result of the Gini index applied to these bins” (Villanueva-Rivera, Pijanowski & Villanueva-Rivera, 2018)                                        | Birds and terrestrial biota (in forest, agricultural land and urban areas, Indiana, United States) | Villanueva-Rivera et al. (2018)                           |
| Acoustic complexity index (ACI)               | “Calculated on the basis of a matrix of the intensities extrapolated from the spectrogram (divided into temporal steps and frequency bins), the ACI calculates the absolute difference between two adjacent values of intensity in a single frequency bin” (Pieretti et al., 2011)                                        | Birds (in temperate woodland, Italy)                                                               | Pieretti et al. (2011)                                    |
| Acoustic diversity index (ADI)                | “Calculated by dividing the spectrogram into bins (default 10) and taking the proportion of the signals in each bin above a threshold (default −50 dBFS). The ADI is the result of the Shannon index applied to these bins” (Villanueva-Rivera et al., 2018)                                                              | Birds and terrestrial biota (in forest, agricultural land and urban areas, Indiana, United States) | Villanueva-Rivera, Pijanowski, Doucette, and Pekin (2011) |
| Bioacoustic index (BI)                        | Calculated as the “area under each curve included all frequency bands associated with the dB value that was greater than the minimum dB value for each curve. The area values are thus a function of both the sound level and the number of frequency bands used by the avifauna” (Boelman, Asner, Hart, & Martin, 2007)  | Birds and terrestrial biota (in forest, savanna, woodland and shrubland, Hawaii, United States)    | Boelman et al. (2007)                                     |
| Normalized difference soundscape index (NDSI) | Seeks to “estimate the level of anthropogenic disturbance on the soundscape by computing the ratio of human-generated (anthrophony) to biological (biophony) acoustic components found in field collected sound samples” (Kasten, Gage, Fox, & Joo, 2012)                                                                 | Birds and terrestrial biota (on an island in Twin Lakes, MI, United States)                        | Kasten et al. (2012)                                      |

Reference: Greenhalgh, J. A., Genner, M. J., Jones, G., & Desjonquères, C. (2020). The role of freshwater bioacoustics in ecological research. Wiley Interdisciplinary Reviews: Water, 7(3), e1416. https://doi.org/10.1002/wat2.1416  

#### Calculating acoustic indices in R Studio                                              

Multiple acoustic indices can be calculated in bulk using the'soundecology' package in R Studio (Villanueva-Rivera & Pijanowski, 2018)

```
By Jack A. Greenhalgh, 19th May, 2025.
Department of Biology, McGill University, 1205 Dr Penfield Ave, Montreal, Quebec, H3A 1B1, Canada.

#### Required Libraries ####
library(seewave)
library(tuneR)
library(soundecology)
library(tools)
library(parallel)
library(doParallel)
library(foreach)

#### Set Directory ####
audio_dir <- file.path("INSERT YOUR DIRECTORY PATH HERE")

# List all .wav files in the target directory
files <- list.files(path = audio_dir, pattern = "(?i)\\.wav$", full.names = TRUE)

#### Part 1: Soundecology indices (Parallelized) ####

# Define indices to compute
soundecology_indices <- c("acoustic_complexity", "acoustic_diversity", "acoustic_evenness",
                          "bioacoustic_index", "H", "ndsi")

# Set up parallel backend
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Compute all indices in parallel
foreach(index = soundecology_indices, .packages = "soundecology") %dopar% {
  result_file <- paste0(index, ".csv")
  multiple_sounds(audio_dir, resultfile = result_file, soundindex = index)
}

# List generated result files
resultfiles <- list.files(path = audio_dir, pattern = "\\.csv$", full.names = FALSE)

# Create an empty list to store the data frames
data_list <- list()

# Loop over each result file and read into the list
for (file in resultfiles) {
  file_path <- file.path(audio_dir, file)
  
  if (file.exists(file_path)) {
    df <- read.csv(file_path)
    if (!"FILENAME" %in% names(df)) next  # Skip unrelated CSVs
    df$FILENAME <- tolower(df$FILENAME)  # Normalize for merging
    index_name <- file_path_sans_ext(file)
    
    if ("LEFT_CHANNEL" %in% names(df)) {
      colnames(df)[colnames(df) == "LEFT_CHANNEL"] <- index_name
    }
    
    data_list[[file]] <- df
    cat("\nData for file:", file, "\n")
    print(head(df))
  }
}

# Remove duplicates by FILENAME
data_list_cleaned <- lapply(data_list, function(df) df[!duplicated(df$FILENAME), ])

# Keep only FILENAME and the index column
data_list_subset <- lapply(names(data_list_cleaned), function(file) {
  df <- data_list_cleaned[[file]]
  index_name <- file_path_sans_ext(file)
  if (index_name %in% names(df)) {
    df_subset <- df[, c("FILENAME", index_name)]
    return(df_subset)
  }
})
names(data_list_subset) <- file_path_sans_ext(resultfiles)

# Merge all soundecology index data
merged_soundecology_indices <- Reduce(function(x, y) merge(x, y, by = "FILENAME", all = TRUE), data_list_subset)
print(head(merged_soundecology_indices))

# Save merged Soundecology data
write.csv(merged_soundecology_indices, file.path(audio_dir, "merged_indices.csv"), row.names = FALSE)

#### Part 2: Seewave indices (Parallelized) ####

# Stop if no audio files found
if (length(files) == 0) stop("❌ No .wav files found in the directory.")

# Parallel processing using foreach
results_list <- foreach(file = files, .combine = rbind, .packages = c("seewave", "tuneR")) %dopar% {
  cat("Processing:", basename(file), "\n")
  
  tryCatch({
    wav <- readWave(file)
    
    if (wav@stereo) wav <- mono(wav, "left")
    if (wav@samp.rate != 44100) wav <- downsample(wav, 44100)
    
    spec_stats <- tryCatch(specprop(meanspec(wav, plot = FALSE)), error = function(e) NULL)
    
    data.frame(
      FILENAME = tolower(basename(file)),
      mean_freq = if (!is.null(spec_stats)) spec_stats$mean else NA,
      sd_freq = if (!is.null(spec_stats)) spec_stats$sd else NA,
      mode_freq = if (!is.null(spec_stats)) spec_stats$mode else NA,
      median_freq = if (!is.null(spec_stats)) spec_stats$median else NA,
      Q25 = if (!is.null(spec_stats)) spec_stats$Q25 else NA,
      Q75 = if (!is.null(spec_stats)) spec_stats$Q75 else NA,
      IQR = if (!is.null(spec_stats)) spec_stats$IQR else NA,
      skewness = if (!is.null(spec_stats)) spec_stats$skewness else NA,
      kurtosis = if (!is.null(spec_stats)) spec_stats$kurtosis else NA
    )
  }, error = function(e) {
    cat("⚠️ Error processing", basename(file), ":", e$message, "\n")
    NULL
  })
}

# Stop parallel backend
stopCluster(cl)

# Save Seewave indices
if (!is.null(results_list) && nrow(results_list) > 0) {
  write.csv(results_list, file.path(audio_dir, "acoustic_indices_summary.csv"), row.names = FALSE)
  cat("✅ Done! Results saved to 'acoustic_indices_summary.csv'\n")
} else {
  cat("⚠️ No valid audio files processed.\n")
}

#### Final Merge: Seewave + Soundecology ####

head(results_list)
head(merged_soundecology_indices)

# Ensure FILENAME columns are character, trimmed, and lowercase
results_list$FILENAME <- tolower(trimws(as.character(results_list$FILENAME)))
merged_soundecology_indices$FILENAME <- tolower(trimws(as.character(merged_soundecology_indices$FILENAME)))

# Merge
merged_all_indices <- merge(results_list, merged_soundecology_indices, by = "FILENAME")

# Save to file
write.csv(merged_all_indices, file.path(audio_dir, "all_acoustic_indices_merged.csv"), row.names = FALSE)

```
It is also possible to specify the min and max frequency, cluster size, and FFT window size. 

For example:
```
	acoustic_complexity(soundfile, min_freq = 2, max_freq = 10, j = 5, fft_w = 512)
```

| Argument  | Description                                                                                             |
|-----------|---------------------------------------------------------------------------------------------------------|
| soundfile | an object of class Wave loaded with the function readWave of the tuneR package.                         |
| min_freq  | miminum frequency to use when calculating the value, in Hertz. The default is NA.                       |
| max_freq  | maximum frequency to use when calculating the value, in Hertz. The default is the maximum for the file. |
| j         | the cluster size, in seconds.                                                                           |
| fft_w     | FFT window to use.                                                                                      |

Referecne: Villanueva-Rivera, L. J., & Pijanowski, B. C. (2018). Soundecology: Soundscape ecology. R package version 1.3.3. https://CRAN.R-project.org/package=soundecology

The package 'seewave' can also be used to calculate additional acoustic indcies, see [this link](https://cran.r-project.org/web/packages/seewave/seewave.pdf) for more information. However, 'seewave' is not optimally designed for the analysis of large datasets and is therefore not described in detail here. 

#### Calculating acoustic indices in Python                                         

Download Anaconda via: https://www.anaconda.com/download and then run the code below in Spyder. 

```
By Jack A. Greenhalgh, 10th June, 2025.
Department of Biology, McGill University, 1205 Dr Penfield Ave, Montreal, Quebec, H3A 1B1, Canada.

import os
import pandas as pd
from maad import sound
import maad.features.alpha_indices as ai
from tqdm import tqdm  # Optional: shows progress bar

# Main directory containing multiple folders with audio files
main_directory = r"C:\Users\jgreenhalgh\Downloads\Gault\Simulated 1st deployment"

# Loop through all subfolders in the main directory
for foldername in os.listdir(main_directory):
    folder_path = os.path.join(main_directory, foldername)
    
    # Proceed only if it's a directory
    if os.path.isdir(folder_path):
        print(f"\nProcessing folder: {foldername}")
        results = []
        
        # List all .wav files in this subfolder
        wav_files = [f for f in os.listdir(folder_path) if f.lower().endswith('.wav')]
        
        for filename in tqdm(wav_files, desc=f"Files in {foldername}"):
            filepath = os.path.join(folder_path, filename)
            try:
                # Load the audio
                s, fs = sound.load(filepath)

                # Generate spectrogram
                Sxx_power, tn, fn, _ = sound.spectrogram(s, fs)

                # Compute spectral alpha indices
                spectral = ai.all_spectral_alpha_indices(Sxx_power, tn, fn)
                spectral_dict = spectral[0].iloc[0].to_dict()

                # Compute temporal alpha indices
                temporal = ai.all_temporal_alpha_indices(s, fs)
                temporal_dict = temporal.iloc[0].to_dict()

                # Merge both sets of indices and tag filename with folder
                combined = {**spectral_dict, **temporal_dict}
                combined['filename'] = f"{foldername}_{filename}"  # Add folder prefix to filename

                results.append(combined)

            except Exception as e:
                print(f"Error processing {filename} in {foldername}: {e}")

        # If results were gathered, save to CSV named after the folder
        if results:
            df = pd.DataFrame(results)
            cols = ['filename'] + [c for c in df.columns if c != 'filename']
            df = df[cols].sort_values(by='filename').reset_index(drop=True)

            output_csv = os.path.join(folder_path, f"{foldername}_alpha_acoustic_indices_results.csv")
            df.to_csv(output_csv, index=False)

            print(f"Results saved for folder '{foldername}' at:\n{output_csv}")
        else:
            print(f"No audio files processed in folder '{foldername}'.")
```

#### Calculating acoustic indices using Kaleidoscope Pro

The analysis of large acoustic datasets in best handled by a program called Kaleidoscope Pro, see https://www.wildlifeacoustics.com/products/kaleidoscope/sound-level-analysis for more information.

![Kaleidoscope Pro](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/17b22f34-81c9-4899-ac74-7bfab369d5ee)

All examples of code henceforth will run using output files generated by Kaleidoscope Pro.  

### Tapestry plots: Visualising long-term passive acoustic monitoring data 

The visualisation of long-term passive acoustic monitoring data consisting of multiple acoustic indices can be challenging. In this section, I will explain the process of creating tapestry plots for the automated and intuitive  visualisation of long-term passive acoustic monitoring data.

![Tapestry plot workflow](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/be6e5223-0dd8-4993-9063-30a2f53c765c)
(Above) Tapestry plot workflow summary. 

First, the Kaleidoscope Pro output files are imported into R and subset according to recorder. Next, highly correlated variables are removed and a global PCA is calculated. The first three PCAs are the extracted and sclaed between 0-1 to form a red, green, and blue channel. Finally, these scaled PCAs are converted to a unique HEX colour code and plotted against date and time of day to create the tapestry plots. 

![J016 tapestry plot GRB (213)](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/ba3385bf-6de4-43fa-a258-68beed61a965)
(Above) Habitats characterised by dense tree cover at lower elevations, such as beech dominated deciduous woodland (1,250 m) and Scots pine forests on the valley slopes (1,250 m – 1,800 m), showed evidence of a significant amount of bird song represented as dark green, with clear demarcations for the dawn and dusk choruses. Bush cricket stridulation was also detected in the evenings (20:00 – 00:00) between July and September represented by orange.  


#### Tapestry plots: Step 1: Cleaning, subsetting, and running a global PCA. 

```
# Remove all objects in the global environment
set.seed(123)
rm(list = ls())
setwd("C:/Users/Administrador/Downloads/R/Ordesa/Data analysis Oct")

# Could it be that datetimes with multiple copies of the same are confusing the order. 
# Try ordering with date times and FOLDER.  

library(readxl)
library(dplyr)

# Read the Excel files
ordesa_1 <- read_excel("acousticindex_Ordesa1_alldata.xlsx")
ordesa_2 <- read_excel("acousticindex_Ordesa2_alldata.xlsx")
ordesa_3 <- read_excel("acousticindex_Ordesa1_newdata.xlsx")
ordesa_4 <- read_excel("acousticindex_Ordesa2_newdata.xlsx")

# Remove the 'ADI' and 'AEI' columns
ordesa_3 <- ordesa_3 %>%
  select(-ADI, -AEI)

# Remove the 'ADI' and 'AEI' columns
ordesa_4 <- ordesa_4 %>%
  select(-ADI, -AEI)

All_Data <- rbind(ordesa_1, ordesa_2, ordesa_3, ordesa_4)

# Display the combined cleaned data
head(All_Data)

# Omit NAs
All_Data <- na.omit(All_Data)

write.csv(All_Data, "All_Data.csv")

####Sub-setting data frames####

# Load all data as a shortcut 

library(readxl)
library(dplyr)

All_Data <- read.csv("All_Data.csv")

# Sub-setting to keep only rows with DURATION >= 60
All_Data <- All_Data[All_Data$DURATION >= 60, ]

All_Data <- na.omit(All_Data)

# Convert the IN.FILE to a date-time object
All_Data$Datetime <- as.POSIXct(gsub("_", " ", gsub(".WAV", "", All_Data$IN.FILE)), 
                                format = "%Y%m%d %H%M%S")

# Sort the data frame by the new Datetime column
All_Data <- All_Data[order(All_Data$Datetime), ]

head(All_Data)

min(All_Data$DATE)
max(All_Data$DATE)

# Load necessary library

# View the updated data frame
head(All_Data)

All_Data <- na.omit(All_Data)

# Remove highly correlated variables
All_Data <- All_Data %>%
  select(-X, -EVN, -SFM, -Q75, -SH, -MEAN, -IQR, -MEDIAN, -SEM, -SKEW, -CENT)

# Retrieving Global PCA subset

# Define the date range
start_date <- as.POSIXct("2023-08-12", format = "%Y-%m-%d", tz = "UTC")
end_date <- as.POSIXct("2024-08-12", format = "%Y-%m-%d", tz = "UTC")

# Filter the data based on the Datetime column
Global_PCA_Subset <- All_Data[All_Data$Datetime >= start_date & All_Data$Datetime <= end_date, ]
Global_PCA_Subset <- na.omit(Global_PCA_Subset)

# Add a 'Type' column with 'Baseline'
Global_PCA_Subset$Type <- "Baseline"

# View the filtered data frame
head(Global_PCA_Subset)

# Retrieving new data subset

# Define the start and end dates for exclusion
start_date <- as.POSIXct("2023-08-12", format = "%Y-%m-%d", tz = "UTC")
end_date <- as.POSIXct("2024-08-12", format = "%Y-%m-%d", tz = "UTC")

# Clean column names just to ensure there are no leading/trailing spaces
colnames(All_Data) <- trimws(colnames(All_Data))

# Ensure Datetime is in POSIXct format (it looks like it already is, but just for safety)
All_Data <- All_Data %>%
  mutate(Datetime = as.POSIXct(Datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC"))

# Check if Datetime is properly created
print(head(All_Data$Datetime))  # View the first few Datetime values

# Check if the Datetime column exists
if ("Datetime" %in% colnames(All_Data)) {
  # Filter out entries within the date range
  New_Data_Subset <- All_Data %>%
    filter(!is.na(Datetime) & !(Datetime >= start_date & Datetime <= end_date))
  
  # Add a 'Type' column with 'New data'
  New_Data_Subset$Type <- "New data"
  
  # View the first few rows of the new subset to confirm
  print(head(New_Data_Subset))
} else {
  stop("Datetime column is not found after mutation.")
}

# View the filtered new data subset
head(New_Data_Subset)

New_Data_Subset <- New_Data_Subset[order(New_Data_Subset$IN.FILE), ]
new_data_datetime_df <- as.data.frame(New_Data_Subset$Datetime)

Global_PCA_Subset <- Global_PCA_Subset[order(Global_PCA_Subset$IN.FILE), ]
global_datetime_df <- as.data.frame(Global_PCA_Subset$Datetime)

Global_PCA_Subset <- Global_PCA_Subset %>%
  select(-FOLDER, -IN.FILE, -CHANNEL, -OFFSET, -DURATION, -DATE, -TIME, -HOUR)

New_Data_Subset <- New_Data_Subset %>%
  select(-FOLDER, -IN.FILE, -CHANNEL, -OFFSET, -DURATION, -DATE, -TIME, -HOUR)

Global_PCA_Subset <- Global_PCA_Subset %>%
  select(-Datetime, -Type)

New_Data_Subset <- New_Data_Subset %>%
  select(-Datetime, -Type)

####Scaling data####

library(scales)

head(Global_PCA_Subset)

# Perform min-max scaling for all columns except the 'Type' column
Global_PCA_Subset_Scaled <- Global_PCA_Subset %>%
  mutate(across(everything(), ~ rescale(.x, to = c(0, 1))))

New_Data_Subset_Scaled <- New_Data_Subset %>%
  mutate(across(everything(), ~ rescale(.x, to = c(0, 1))))

# View the first few rows of both data frames
head(Global_PCA_Subset_Scaled)
head(New_Data_Subset_Scaled)

#PCA.

#### Applying global PCA to new data ####

# Step 1: Perform PCA on Global_PCA_Baseline
set.seed(123)

# Re-center the data (use the mean of the baseline data for both data sets)
pca_baseline <- prcomp(Global_PCA_Subset_Scaled, center = TRUE, scale. = FALSE)

# Extract the first three components for the baseline data
pca_components_baseline <- pca_baseline$x[, 1:3]  # First three principal components

summary(pca_baseline)
print(pca_baseline)

# Step 2: Project All_New_Data onto the PCA space
# Make sure to use the same center (mean) and scale as the baseline PCA
all_new_data_matrix <- as.matrix(New_Data_Subset_Scaled)

# Center the new data using the baseline PCA mean (since PCA already centers baseline data)
new_data_centered <- scale(all_new_data_matrix, center = pca_baseline$center, scale = F)

# Project the new data onto the PCA space using the baseline rotation matrix
pca_new_data <- new_data_centered %*% pca_baseline$rotation[, 1:3]  # First three PCA loadings

# Combine Results (Optional)
# Combine PCA components into data frames
pca_baseline_df <- as.data.frame(pca_components_baseline)
pca_new_data_df <- as.data.frame(pca_new_data)

# Add Extracted_Date and Extracted_Time to baseline_df
Global_PCA_Combined_df <- cbind(pca_baseline_df, 
                                Datetime = global_datetime_df)
Global_PCA_Combined_df$Type <- "Baseline"
colnames(Global_PCA_Combined_df)[colnames(Global_PCA_Combined_df) == "Global_PCA_Datetimes"] <- "Datetime"

# Add Extracted_Date and Extracted_Time to new_data_df
New_Data_combined_df <- cbind(pca_new_data_df, 
                              Datetime = new_data_datetime_df)
New_Data_combined_df$Type <- "New data"
colnames(New_Data_combined_df)[colnames(New_Data_combined_df) == "New_Data_Datetimes"] <- "Datetime"

head(Global_PCA_Combined_df)
head(New_Data_combined_df)

# Rename 'New_Data_Subset$Datetime' column to 'Datetime'
colnames(New_Data_combined_df)[colnames(New_Data_combined_df) == "New_Data_Subset$Datetime"] <- "Datetime"
colnames(Global_PCA_Combined_df)[colnames(Global_PCA_Combined_df) == "Global_PCA_Subset$Datetime"] <- "Datetime"

merged_data <- rbind(Global_PCA_Combined_df, New_Data_combined_df)

# Order the data frames by the Datetime column
merged_data <- merged_data[order(merged_data$Datetime), ]

# Calculate min and max values for PC1, PC2, and PC3, ignoring NAs
min_PC1 <- min(merged_data$PC1, na.rm = TRUE)
max_PC1 <- max(merged_data$PC1, na.rm = TRUE)

min_PC2 <- min(merged_data$PC2, na.rm = TRUE)
max_PC2 <- max(merged_data$PC2, na.rm = TRUE)

min_PC3 <- min(merged_data$PC3, na.rm = TRUE)
max_PC3 <- max(merged_data$PC3, na.rm = TRUE)

cat("PC1: Min =", min_PC1, "Max =", max_PC1, "\n")
cat("PC2: Min =", min_PC2, "Max =", max_PC2, "\n")
cat("PC3: Min =", min_PC3, "Max =", max_PC3, "\n")

# Normalize the PCs to the range [0, 1]
merged_data <- merged_data %>%
  mutate(
    norm_PC1 = (PC1 - min(PC1, na.rm = TRUE)) / (max(PC1, na.rm = TRUE) - min(PC1, na.rm = TRUE)),
    norm_PC2 = (PC2 - min(PC2, na.rm = TRUE)) / (max(PC2, na.rm = TRUE) - min(PC2, na.rm = TRUE)),
    norm_PC3 = (PC3 - min(PC3, na.rm = TRUE)) / (max(PC3, na.rm = TRUE) - min(PC3, na.rm = TRUE))
  )

# Display the first few rows of the updated cleaned data frame
head(merged_data)

# Combine the two data frames by columns
merged_data <- cbind(merged_data, FOLDER = All_Data$FOLDER)
merged_data <- cbind(merged_data, IN.FILE = All_Data$IN.FILE)

# Convert the IN.FILE to a date-time object
merged_data$Datetime <- as.POSIXct(gsub("_", " ", gsub(".WAV", "", merged_data$IN.FILE)), 
                                format = "%Y%m%d %H%M%S")

# Sort the data frame by the new Datetime column
merged_data <- merged_data[order(merged_data$Datetime), ]

# View the first few rows of the updated merged_data
head(merged_data)

# Identify rows where Datetime is NA
na_rows <- is.na(merged_data$Datetime)

# Extract the date and time from the IN.FILE column for rows with NA in Datetime
# Assuming the format is YYYYMMDD_HHMMSS in IN.FILE
if (any(na_rows)) {
  # Extract the datetime information from the IN.FILE column (assuming the pattern YYYYMMDD_HHMMSS)
  extracted_datetime <- strptime(str_extract(merged_data$IN.FILE[na_rows], "\\d{8}_\\d{6}"),
                                 format = "%Y%m%d_%H%M%S",
                                 tz = "UTC")
  
  # Replace the NA values in Datetime with the extracted datetime values
  merged_data$Datetime[na_rows] <- extracted_datetime
  
  # Check if the replacement was successful
  if (any(is.na(merged_data$Datetime))) {
    cat("Warning: There are still NA values in the 'Datetime' column after filling.\n")
  } else {
    cat("All missing 'Datetime' values have been filled.\n")
  }
}

# Ensure that the data is ordered by Datetime after filling
merged_data <- merged_data[order(merged_data$Datetime), ]

# View the updated data
head(merged_data)
```

#### Subsetting by AudioMoth #### 

At this point it is important to understand the data structure that this code was written for. 

Our survey area was divided into 8 zones, containing 23 sites, 5 habitats, and 36 recorders.  
| Zone        | Site | Habitat   | AudioMoth |
|-------------|------|-----------|-----------|
| Z01_PARADOR | S001 | MATORRAL  | J008      |
| Z01_PARADOR | S002 | PINO SILV | J003      |
| Z01_PARADOR | S003 | PINO SILV | J016      |
| Z01_PARADOR | S016 | MATORRAL  | J009      |

As such, the folders containing the audio files that were downloaded from the field were named using the following structure: 
Z02-PRADERASUR_S004-HAYABE_J001_220712-220901, which forms the FOLDER column in the final dataset. 

The next section of code finds unique folder names in the FOLDER column to subset AudioMoths within the kaleidoscope output files. 

Adapt the code accordingly to your data structure to find unique folder names. However, so long as the data from each AudioMoth are grouped by a unique name within a column called FOLDER, the following code should automatically subset them for you.
```

# Find unique values in the 'FOLDER' column
unique_folders <- unique(merged_data$FOLDER)

# Create a list to store data frames for each subset
subset_data_frames <- lapply(unique_folders, function(folder) {
  merged_data[merged_data$FOLDER == folder, ]
})
names(subset_data_frames) <- unique_folders

# Function to filter data frames based on folder names containing specific patterns (AudioMoth)
filter_data_frames <- function(list_of_data_frames, pattern) {
  filtered_data_frames <- lapply(names(list_of_data_frames), function(name) {
    df <- list_of_data_frames[[name]]
    # Check if the folder name contains the specified pattern
    if (grepl(pattern, name)) {
      return(df)  # Return the data frame if the pattern matches
    } else {
      return(NULL)  # Return NULL if there is no match
    }
  })
  
  # Set names for filtered data frames and filter out NULL values
  names(filtered_data_frames) <- names(list_of_data_frames)
  filtered_data_frames <- Filter(Negate(is.null), filtered_data_frames)
  
  return(filtered_data_frames)  # Return the filtered data frames
}

# Generate patterns for AudioMoth numbers (1:36)
patterns <- paste0("J", formatC(1:36, width = 3, format = "d", flag = "0"))

# Create a list to group data frames based on Audiomoth number patterns
grouped_data_frames <- lapply(patterns, function(pattern) {
  filter_data_frames(subset_data_frames, pattern)
})
names(grouped_data_frames) <- patterns

# Create a list to store concatenated data frames for each pattern
single_data_frames <- lapply(grouped_data_frames, function(data_frames_list) {
  # Remove NULL values and concatenate data frames
  data_frames_list <- Filter(Negate(is.null), data_frames_list)
  if (length(data_frames_list) > 0) {
    return(do.call(rbind, data_frames_list))
  } else {
    return(NULL)  # Return NULL if there are no matching data frames
  }
})

# Remove NULL entries from single_data_frames
single_data_frames <- Filter(Negate(is.null), single_data_frames)

# Output the final concatenated data frames grouped by AudioMoth patterns
single_data_frames

# Assuming single_data_frames contains the grouped data frames

# Directory to save CSV files
output_directory <- "C:/Users/Administrador/Downloads/R/Ordesa/Data analysis Oct"

# Create the directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# Loop through each data frame in single_data_frames and save as CSV
for (pattern in names(single_data_frames)) {
  # Get the corresponding data frame
  df <- single_data_frames[[pattern]]
  
  # Check if the data frame is not NULL and has data
  if (!is.null(df) && nrow(df) > 0) {
    # Construct the file name by replacing characters as needed
    filename <- paste0("AudioMoth_", pattern, ".csv")
    
    # Define the full file path
    file_path <- file.path(output_directory, filename)
    
    # Write the data frame to a CSV file
    write.csv(df, file = file_path, row.names = FALSE)
    
    # Optionally print a message to confirm the export
    cat("Exported:", file_path, "\n")
  } else {
    cat("No data for:", pattern, "\n")  # Indicate if there's no data to export
  }
}

```

#### Tapestry plots: Step 2: Automatically accounting for missing data and plotting.  
```
# Set seed, locale, and working directory
set.seed(123)
rm(list = ls())
Sys.setlocale("LC_TIME", "English")
setwd("C:/Users/Administrador/Downloads/R/Ordesa/Data analysis Oct")

# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(scales)
library(dplyr)

# Get the list of all AudioMoth files that match the pattern J001 to J036
files <- list.files(pattern = "AudioMoth_J0(0[1-9]|[1-2][0-9]|3[0-6]).csv")

# Define the function to fill missing dates
fill_missing_dates <- function(data, expected_dates) {
  all_dates <- data.frame(Datetime = expected_dates)
  merged_data <- merge(all_dates, data, by = "Datetime", all.x = TRUE)
  merged_data[is.na(merged_data)] <- 0
  return(merged_data)
}

# Loop through each file in the list
for (file in files) {
  # Read the data
  data <- read.csv(file)
  
  # Ensure Datetime is in the correct format
  data$Datetime <- as.character(data$Datetime)
  
  # Remove rows with NA in the Datetime column
  data <- na.omit(data)
  
  # Identify rows that contain only the date (no time)
  only_date_rows <- grepl("^\\d{4}-\\d{2}-\\d{2}$", data$Datetime)
  data$Datetime[only_date_rows] <- paste0(data$Datetime[only_date_rows], " 00:00:00")
  
  # Convert Datetime back to POSIXct
  data$Datetime <- as.POSIXct(data$Datetime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # Sort data by Datetime
  data <- data[order(data$Datetime), , drop = FALSE]
  
  # Calculate the expected sequence of dates with a 10-minute interval
  expected_dates <- seq(from = min(data$Datetime), to = max(data$Datetime), by = "10 min")
  
  # Identify missing dates
  missing_dates <- setdiff(expected_dates, data$Datetime)
  
  # Create a data frame for missing dates with proper POSIXct format
  missing_data <- data.frame(
    Datetime = as.POSIXct(missing_dates, origin = "1970-01-01", tz = "UTC"),
    norm_PC1 = 0,
    norm_PC2 = 0,
    norm_PC3 = 0,
    Type = NA,
    HEX_Codes = NA,
    FOLDER = NA
  )
  
  # Combine original data with missing data
  combined_data <- bind_rows(data, missing_data)
  combined_data <- combined_data[order(combined_data$Datetime), ]
  
  # Generate HEX color codes based on norm_PC1 (R), norm_PC2 (B), and norm_PC3 (G)
  combined_data$HEX_Codes <- rgb(
    red = combined_data$norm_PC2,
    green = combined_data$norm_PC1,
    blue = combined_data$norm_PC3,
    maxColorValue = 1
  )
  
  # Split the datetime column into date and time
  combined_data <- combined_data %>%
    mutate(
      date = as.Date(Datetime),
      time = format(Datetime, "%H:%M:%S")
    ) %>%
    select(-Datetime)
  
  # Convert date and time to proper data types
  combined_data$date <- as.Date(combined_data$date, format = "%d/%m/%Y")
  combined_data$time <- as.POSIXct(combined_data$time, format = "%H:%M")
  
  # Plot
  p <- ggplot(combined_data, aes(x = time, y = date)) +
    geom_tile(aes(fill = `HEX_Codes`)) +
    labs(x = "Time of day", y = "Date") +
    scale_x_datetime(date_breaks = "2 hours", date_labels = "%H:%M") +
    scale_y_date(date_breaks = "1 month", date_labels = "%b %Y") +
    theme_bw() +
    scale_fill_identity()
  
  # Generate file name for saving the plot (replace .csv with .png)
  plot_filename <- gsub(".csv", ".png", file)
  
  # Save the plot
  ggsave(filename = plot_filename, plot = p, width = 6, height = 6)
}

```

#### False-colour spectrograms: Plotting 24 hrs of acoustic data with noise reduction. 

![Ordesa Summer 24 hr spec](https://github.com/user-attachments/assets/48f0378d-b104-4dca-ae11-180b91eb9925)

```
library(tuneR)
library(seewave)
library(signal)
library(entropy)

# Define the directory containing the .wav files
directory <- "C:/Users/Administrador/Downloads/Ordesa sound files/Daily Spectrograms/J018/Feb 2024"

# Print directory path to confirm
print(paste("Checking directory:", directory))

# Check if the directory exists
if (!dir.exists(directory)) {
  stop("The specified directory does not exist.")
}

# List all .WAV files in the directory (considering case sensitivity)
file_list <- list.files(directory, pattern = "\\.WAV$", full.names = TRUE)

# Print file list to debug
print("List of .WAV files found:")
print(file_list)

# Check if there are any .WAV files
if (length(file_list) == 0) {
  stop("No .WAV files found in the directory.")
}

# Initialize lists to store results
entropy_list <- list()
aci_list <- list()
background_noise_list <- list()

# Function definitions for calculating metrics
calculate_entropy <- function(freq_bin) {
  total_energy <- sum(freq_bin)
  if (total_energy == 0) {
    return(0)  # Avoid division by zero
  }
  prob_distribution <- freq_bin / total_energy
  shannon_entropy <- -sum(prob_distribution * log(prob_distribution + 1e-10))
  return(shannon_entropy)
}

calculate_aci <- function(freq_bin) {
  local_maxima <- sum(diff(sign(diff(freq_bin))) == -2)
  aci_value <- local_maxima / length(freq_bin)
  return(aci_value)
}

calculate_background_noise <- function(freq_bin) {
  median_noise <- median(freq_bin)
  return(median_noise)
}

# Loop through each .wav file
for (audio_file in file_list) {
  # Load the audio file
  wave <- readWave(audio_file)
  
  # Extract the audio signal and the sampling rate
  y <- wave@left
  sr <- wave@samp.rate
  
  # Parameters for STFT
  n_fft <- 16384  # Number of FFT points
  hop_length <- n_fft / 2  # Hop length (50% overlap)
  
  # Perform Short-Time Fourier Transform (STFT) using specgram
  specgram_result <- specgram(y, n = n_fft, Fs = sr, overlap = n_fft - hop_length)
  
  # Get the magnitude spectrogram
  S <- abs(specgram_result$S)
  
  # Normalize the spectrogram
  S <- S / max(S)
  
  # Convert to dB scale
  S_db <- 10 * log10(S + 1e-10)
  
  # Calculate metrics
  entropy_values <- apply(S_db, 1, calculate_entropy)
  aci_values <- apply(S_db, 1, calculate_aci)
  background_noise_values <- apply(S_db, 1, calculate_background_noise)
  
  # Compute frequency bins
  freq_bins <- seq(0, sr / 2, length.out = n_fft / 2 + 1)[1:nrow(S_db)]
  
  # Ensure the length of freq_bins matches the length of metric values
  if (length(freq_bins) != length(entropy_values) ||
      length(freq_bins) != length(aci_values) ||
      length(freq_bins) != length(background_noise_values)) {
    warning(paste("Length mismatch in file:", basename(audio_file)))
    next
  }
  
  # Create data frames for each metric
  entropy_df <- data.frame(Frequency = freq_bins, Entropy = entropy_values, File = basename(audio_file))
  aci_df <- data.frame(Frequency = freq_bins, ACI = aci_values, File = basename(audio_file))
  background_noise_df <- data.frame(Frequency = freq_bins, Background_Noise = background_noise_values, File = basename(audio_file))
  
  # Append results to lists
  entropy_list[[basename(audio_file)]] <- entropy_df
  aci_list[[basename(audio_file)]] <- aci_df
  background_noise_list[[basename(audio_file)]] <- background_noise_df
}

# Combine all results into single data frames
entropy_all <- do.call(rbind, entropy_list)
aci_all <- do.call(rbind, aci_list)
background_noise_all <- do.call(rbind, background_noise_list)

# Print resulting data frames
print(head(entropy_all))
print(head(aci_all))
print(head(background_noise_all))

# Load necessary library
library(dplyr)

# Extract time and create the 'Time' column for entropy_all, and remove the 'File' column
entropy_all <- entropy_all %>%
  mutate(
    # Extracting the time part from the File column
    TempTime = substr(File, 10, 15),
    # Converting extracted time to hh:mm format
    Time = paste0(substr(TempTime, 1, 2), ":", substr(TempTime, 3, 4))
  ) %>%
  select(-TempTime, -File)  # Remove the intermediate 'TempTime' and 'File' columns

# View the updated data frame
head(entropy_all)


# Extract time and create the 'Time' column for aci_all, and remove the 'File' column
aci_all <- aci_all %>%
  mutate(
    # Extracting the time part from the File column
    TempTime = substr(File, 10, 15),
    # Converting extracted time to hh:mm format
    Time = paste0(substr(TempTime, 1, 2), ":", substr(TempTime, 3, 4))
  ) %>%
  select(-TempTime, -File)  # Remove the intermediate 'TempTime' and 'File' columns

# View the updated data frame
head(aci_all)

# Extract time and create the 'Time' column for background_noise_all, and remove the 'File' column
background_noise_all <- background_noise_all %>%
  mutate(
    # Extracting the time part from the File column
    TempTime = substr(File, 10, 15),
    # Converting extracted time to hh:mm format
    Time = paste0(substr(TempTime, 1, 2), ":", substr(TempTime, 3, 4))
  ) %>%
  select(-TempTime, -File)  # Remove the intermediate 'TempTime' and 'File' columns

# View the updated data frame
head(background_noise_all)

#### Scale acoustic indices ####

# Entropy

# Find min and max of Entropy
min_entropy <- min(entropy_all$Entropy)
max_entropy <- max(entropy_all$Entropy)

# Scale the Entropy column
entropy_all$Scaled_Entropy <- (entropy_all$Entropy - min_entropy) / (max_entropy - min_entropy)

# Print the result
print(data)

# Acoustic complexity

# Find min and max of aci
min_aci <- min(aci_all$ACI)
max_aci <- max(aci_all$ACI)

# Scale the aci column
aci_all$Scaled_ACI <- (aci_all$ACI - min_aci) / (max_aci - min_aci)

# Background noise

# Find min and max of background_noise
min_noise <- min(background_noise_all$Background_Noise)
max_noise <- max(background_noise_all$Background_Noise)

# Scale the background_noise column
background_noise_all$Scaled_Background_Noise <- (background_noise_all$Background_Noise - min_noise) / (max_noise - min_noise)

#### Mapping to RBG ####

# Assuming that the scaled columns are already present in the respective data frames
# Combine the relevant columns into a single data frame

combined_data <- data.frame(
  Scaled_Entropy = entropy_all$Scaled_Entropy,
  Scaled_ACI = aci_all$Scaled_ACI,
  Scaled_Background_Noise = background_noise_all$Scaled_Background_Noise
)

# Initialize RGB_Data with the same number of rows
RGB_Data <- data.frame(matrix(ncol = 1, nrow = nrow(combined_data)))
colnames(RGB_Data) <- "Color"

# Assuming combined_data is your data frame with the scaled variables
library(grDevices) # For the rgb function

# Function to generate RGB values based on different variable combinations
generate_rgb_combinations <- function(data, var1, var2, var3) {
  rgb(data[[var1]], data[[var2]], data[[var3]])
}

# All possible combinations of the variables
combinations <- list(
  c("Scaled_Entropy", "Scaled_ACI", "Scaled_Background_Noise"),
  c("Scaled_Entropy", "Scaled_Background_Noise", "Scaled_ACI"),
  c("Scaled_ACI", "Scaled_Entropy", "Scaled_Background_Noise"),
  c("Scaled_ACI", "Scaled_Background_Noise", "Scaled_Entropy"),
  c("Scaled_Background_Noise", "Scaled_Entropy", "Scaled_ACI"),
  c("Scaled_Background_Noise", "Scaled_ACI", "Scaled_Entropy")
)

# Create an empty list to store the results
RGB_Combinations <- list()

# Loop through the combinations and generate RGB colors
for (i in 1:length(combinations)) {
  combo <- combinations[[i]]
  RGB_Combinations[[i]] <- generate_rgb_combinations(combined_data, combo[1], combo[2], combo[3])
}

# Example: Storing the results in a data frame
RGB_Data <- data.frame(
  Color1 = RGB_Combinations[[1]],
  Color2 = RGB_Combinations[[2]],
  Color3 = RGB_Combinations[[3]],
  Color4 = RGB_Combinations[[4]],
  Color5 = RGB_Combinations[[5]],
  Color6 = RGB_Combinations[[6]]
)

# Print the RGB data frame to check the results
print(RGB_Data)

# Print the result
head(RGB_Data)

head(entropy_all)

##### Fit HEX codes to freq and time #####

# Extract the required columns from entropy_all and RGB_Data
frequency <- entropy_all$Frequency
time <- entropy_all$Time
color <- RGB_Data$Color1

# Combine these columns into a new data frame
combined_df <- data.frame(Frequency = frequency, Time = time, Color = color)

# Print the result
print(head(combined_df))

#### Plot the spectrogram ####

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Ensure 'Time' is in a time format (for ggplot2 to handle it correctly)
combined_df <- combined_df %>%
  mutate(Time = as.POSIXct(Time, format = "%H:%M", tz = "UTC"))

#### Removing background noise ####

# Step 1: Count the occurrences of each HEX code
hex_table <- table(unlist(RGB_Data))

# Convert the table to a data frame for easier manipulation
hex_df <- as.data.frame(hex_table)

# Rename columns for clarity
colnames(hex_df) <- c("HEX", "Count")

# Step 2: Determine the threshold for the dominant 10%
threshold_count <- quantile(hex_df$Count, 0.9)  # Get the count for the 90th percentile

# Identify dominant HEX codes (those that have counts greater than or equal to the threshold)
dominant_hexes <- hex_df$HEX[hex_df$Count >= threshold_count]

# Print dominant HEX codes for reference
print("Dominant HEX Codes:")
print(dominant_hexes)

# Step 3: Replace dominant HEX codes with black in RGB_Data
# Create a function to replace dominant HEX codes
replace_dominant_hex <- function(color) {
  if (color %in% dominant_hexes) {
    return("#000000")  # Replace with black
  } else {
    return(color)  # Keep the original color
  }
}

# Apply the function to replace colors in the first column of RGB_Data
RGB_Data$Adjusted_Color <- sapply(RGB_Data$Color1, replace_dominant_hex)

# Create a new combined_df with the adjusted colors
combined_df_adjusted <- data.frame(
  Frequency = frequency,
  Time = time,
  Color = RGB_Data$Adjusted_Color
)

# Ensure 'Time' is in POSIXct format
combined_df_adjusted <- combined_df_adjusted %>%
  mutate(Time = as.POSIXct(Time, format = "%H:%M", tz = "UTC"))

# Step 4: Plot the adjusted spectrogram
p <- ggplot(combined_df_adjusted, aes(x = Time, y = Frequency)) +
  geom_tile(aes(fill = Color)) +
  scale_fill_identity() +
  scale_x_datetime(date_breaks = "2 hours", date_labels = "%H:%M") +
  labs(x = "Time", y = "Frequency (Hz)") +
  theme_bw()

ggsave("J018 Feb 2024.pdf", plot = p, width = 8, height = 6, dpi = 300)
```

### References

Villanueva-Rivera, L. J., & Pijanowski, B. C. (2018). Soundecology: Soundscape ecology. R package version 1.3.3. https://CRAN.R-project.org/package=soundecology
