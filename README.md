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
#By Jack A. Greenhalgh, 19th May, 2025.
#Department of Biology, McGill University, 1205 Dr Penfield Ave, Montreal, Quebec, H3A 1B1, Canada.

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
#By Jack A. Greenhalgh. October, 2025.
#Department of Biology, McGill University, 1205 Dr Penfield Ave, Montreal, Quebec, H3A 1B1, Canada.

import os
import pandas as pd
from maad import sound
import maad.features.alpha_indices as ai
from tqdm import tqdm

# =============================
# Detect environment and choose executor
# =============================
def in_spyder():
    """Detect if running inside Spyder"""
    return 'SPYDER_ARGS' in os.environ or 'SPYDER_PID' in os.environ

if in_spyder():
    from concurrent.futures import ThreadPoolExecutor as Executor
    use_threads = True
else:
    from concurrent.futures import ProcessPoolExecutor as Executor
    use_threads = False

from concurrent.futures import as_completed

# =============================
# Configurable Parameters
# =============================
NFFT = 1024
noverlap = 512
window = 'hann'
temporal_threshold_db = -75

# Frequency range for spectral indices
# These parameters are passed directly to the spectral index functions
spectral_fmin = 200
spectral_fmax = 700

# Frequency range for optional spectrogram subset
# Only used if subset_spectrogram = True
spectrogram_fmin = 200
spectrogram_fmax = 700

# Option to manually subset spectrogram before index calculation
subset_spectrogram = False  # Set True to crop Sxx manually

# Main directory
main_directory = r"INSERT YOUR DIRECTORY HERE"

# =============================
# File Processing Function
# =============================
def process_file(foldername, folder_path, filename):
    """Compute spectral and temporal alpha indices for a single .wav file"""
    filepath = os.path.join(folder_path, filename)
    try:
        s, fs = sound.load(filepath)
        if s is None or len(s) == 0:
            raise ValueError("Empty or invalid audio file")

        # Compute full spectrogram
        Sxx_power, tn, fn, _ = sound.spectrogram(
            s, fs, window=window, nperseg=NFFT, noverlap=noverlap
        )

        # Optional manual cropping of spectrogram
        if subset_spectrogram:
            freq_mask = (fn >= spectrogram_fmin) & (fn <= spectrogram_fmax)
            if freq_mask.sum() == 0:
                raise ValueError(f"No frequencies found in range {spectrogram_fmin}-{spectrogram_fmax} Hz")
            Sxx_power = Sxx_power[freq_mask, :]
            fn = fn[freq_mask]

        # Compute spectral indices using spectral_fmin/spectral_fmax as function parameters
        spectral = ai.all_spectral_alpha_indices(
            Sxx_power, tn, fn, fmin=spectral_fmin, fmax=spectral_fmax
        )
        spectral_dict = spectral[0].iloc[0].to_dict()

        # Compute temporal indices
        temporal = ai.all_temporal_alpha_indices(s, fs, threshold=temporal_threshold_db)
        temporal_dict = temporal.iloc[0].to_dict()

        # Combine results
        combined = {**spectral_dict, **temporal_dict}
        combined['filename'] = f"{foldername}_{filename}"
        return combined

    except Exception as e:
        print(f"⚠ Error processing {filename} in {foldername}: {e}")
        return None

# =============================
# Main Script
# =============================
def main():
    foldernames = [f for f in os.listdir(main_directory)
                   if os.path.isdir(os.path.join(main_directory, f))]

    for foldername in foldernames:
        folder_path = os.path.join(main_directory, foldername)
        print(f"\n🎧 Processing folder: {foldername}")
        results = []

        wav_files = [f for f in os.listdir(folder_path) if f.lower().endswith('.wav')]
        if not wav_files:
            print(f"⚠ No .wav files found in {folder_path}")
            continue

        executor_type = "Threads" if use_threads else "Processes"
        print(f"Using {executor_type} for parallel processing")

        with Executor() as executor:
            futures = {executor.submit(process_file, foldername, folder_path, f): f for f in wav_files}

            for future in tqdm(as_completed(futures), total=len(futures), desc=f"Files in {foldername}"):
                result = future.result()
                if result:
                    results.append(result)

        if results:
            df = pd.DataFrame(results)

            # Reorder columns
            cols = ['filename'] + [c for c in df.columns if c != 'filename']
            df = df[cols].sort_values(by='filename').reset_index(drop=True)

            # Ensure numeric columns
            df = df.apply(pd.to_numeric, errors='ignore')

            # Save CSV
            output_csv = os.path.join(folder_path, f"{foldername}_fish_alpha_acoustic_indices_results.csv")
            df.to_csv(output_csv, index=False, encoding='utf-8-sig')
            print(f"✔ Results saved for folder '{foldername}' at:\n{output_csv}")
        else:
            print(f"⚠ No audio files processed successfully in folder '{foldername}'.")

# =============================
# Acoustic indices summary table
# =============================
def print_acoustic_indices_summary():
    """
    Print a summary table of the acoustic indices used,
    their key parameters, and a brief description.
    """
    data = [
        ["VARf", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Variance of the frequency bins in the spectrogram"],
        ["KURTf", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Kurtosis of the frequency distribution of the spectrogram"],
        ["NBPEAKS", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Number of spectral peaks in the frequency range"],
        ["BGNf", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Background noise estimate of the frequency spectrum"],
        ["EAS", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Acoustic entropy across frequency bins"],
        ["ECV", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Coefficient of variation of the energy across frequency bins"],
        ["EPS", f"threshold={temporal_threshold_db} dB", "Entropy of the temporal amplitude signal"],
        ["EPS_KURT", f"threshold={temporal_threshold_db} dB", "Kurtosis of temporal entropy"],
        ["ACI", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Acoustic Complexity Index, measuring amplitude variation across time and frequency"],
        ["rBA", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Relative Bioacoustic Index, normalized energy in the frequency band"],
        ["BI", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Bioacoustic Index, measures total energy in the band"],
        ["ADI", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Acoustic Diversity Index, reflects the number of frequency bins with significant activity"],
        ["EVNspMean", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Mean of the Event-based Normalized Spectrogram"],
        ["TFSD", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Temporal Frequency Spectrum Density"],
        ["RAOQ", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Rao's Quadratic Entropy, diversity index in frequency domain"],
        ["AGI", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Acoustic Grouping Index, measures clustering of acoustic events"],
        ["aROI", f"fmin={spectral_fmin}, fmax={spectral_fmax}", "Acoustic Region of Interest index, energy in a specific band"],
        ["MEANt", f"threshold={temporal_threshold_db} dB", "Mean amplitude of temporal signal above threshold"],
        ["SKEWt", f"threshold={temporal_threshold_db} dB", "Skewness of temporal amplitude distribution"],
        ["KURTt", f"threshold={temporal_threshold_db} dB", "Kurtosis of temporal amplitude distribution"],
        ["Ht", f"threshold={temporal_threshold_db} dB", "Shannon entropy of the temporal signal"],
        ["EVNtMean", f"threshold={temporal_threshold_db} dB", "Mean of temporal event-based normalized signal"],
        ["EVNtCount", f"threshold={temporal_threshold_db} dB", "Count of temporal events above threshold"]
    ]

    df_summary = pd.DataFrame(data, columns=["Index", "Parameters", "Description"])
    print("\n🎶 Acoustic Indices Summary Table")
    print(df_summary.to_string(index=False))

    # Optional: save summary table to CSV in main directory
    summary_csv = os.path.join(main_directory, "acoustic_indices_summary.csv")
    df_summary.to_csv(summary_csv, index=False, encoding='utf-8-sig')
    print(f"\n✔ Acoustic indices summary saved as CSV at:\n{summary_csv}")

# =============================
# Entry point
# =============================
if __name__ == "__main__":
    main()
    print_acoustic_indices_summary()
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


#### Subsetting and cleaning acoustic indices data

```
#By Jack A. Greenhalgh & José Joaquín Lahoz Monfort. Feb, 2025.
#Instituto Pirenaico de Ecología (CSIC), Avenida Nuestra Señora de la Victoria, 22700, Jaca, Huesca, España.

# load all libraries needed
library(dplyr)
library(readxl)
library(tuneR)
library(seewave)
library(signal)
library(entropy)
library(ggplot2)
library(viridis)
library(dplyr)
library(caret)
library(corrplot)
library(tidyverse)
library(scales)
library(purrr)
library(gtools)
library(lubridate)

##### Part 1) Read all data and merge in single curated dataset ---------------------

### [Run only once initially] Create raw dataset:
# Read raw .csv files 
HardDrive_1 <- read_excel("INSERT FILE PATH.xlsx")
HardDrive_2 <- read_excel("INSERT FILE PATH.xlsx")

# Merge them in a single large dataframe (accounting for possible differences in columns)
All_Data <- rbind(HardDrive_1, HardDrive_2)

# Remove entries lasting < 60 seconds
All_Data <- All_Data[All_Data$DURATION >= 60, ]

# Remove entries with NAs
All_Data <- na.omit(All_Data)

# Convert the IN FILE to a date-time object
attach(All_Data)
All_Data$Datetime <- as.POSIXct(gsub("_", " ", gsub(".WAV", "", All_Data$'IN FILE')), 
                                format = "%Y%m%d %H%M%S", tz = "UTC")
head(All_Data)
class(All_Data$Datetime)

# Check for NAs in the Datetime variable
sum(is.na(All_Data$Datetime))

# Sort the data frame by the new Datetime column
All_Data <- All_Data[order(All_Data$Datetime), ]

## Calculate and plot correlation matrix
correlation_data <- All_Data[,ADD FIRST AI COLUMN:ADD LAST AI COLUMN]  # Extract ac index columns

# Display the first few rows of the scaled data frame
head(correlation_data)

# Calculate correlation matrix
correlation_matrix <- cor(correlation_data)
correlation_matrix

## Remove acoustic indices with high correlation

threshold <- 0.8  # Set a threshold 

# Plot correlation matrix
corrplot(correlation_matrix, method = "square")

# Find highly correlated variables
highly_correlated <- findCorrelation(correlation_matrix, cutoff = threshold)

# Get names of highly correlated variables
highly_correlated_names <- colnames(correlation_matrix)[highly_correlated]
highly_correlated_names

# Remove highly correlated variables
All_Data <- All_Data %>%
  select(-VAR1, -VAR2, -VAR3, ETC..)

# Save final dataframe as csv
write.csv(All_Data, file="INSERT YOUR FILE NAME HERE.csv")

```

#### Subset dataset by recording device

```
#By Jack A. Greenhalgh & José Joaquín Lahoz Monfort. Feb, 2025.
#Instituto Pirenaico de Ecología (CSIC), Avenida Nuestra Señora de la Victoria, 22700, Jaca, Huesca, España.

# Find unique values in the 'FOLDER' column
unique_folders <- unique(All_Data$FOLDER)

# Create a list to store data frames for each subset
subset_data_frames <- lapply(unique_folders, function(folder) {
  All_Data[All_Data$FOLDER == folder, ]
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
# NOTE: CHOOSE HERE WHAT AUDIOMOTHS WILL BE USED
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
output_directory <- "./1_individual_AM_dataframes"

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
#### Add missing dates to the timeseries 

```
#By Jack A. Greenhalgh & José Joaquín Lahoz Monfort. Feb, 2025.
#Instituto Pirenaico de Ecología (CSIC), Avenida Nuestra Señora de la Victoria, 22700, Jaca, Huesca, España.

# Step 1: Read all CSV files in the directory matching the pattern
files <- list.files(path = "./1_individual_AM_dataframes", pattern = "AudioMoth_J0(0[1-9]|[1-2][0-9]|3[0-6]).csv")

# Step 2: Define fixed date and time ranges for the processing
date_range <- seq(as.Date(start_date), as.Date(end_date), by = "day")

time_range <- format(seq(from = as.POSIXct("00:00:00", format = "%H:%M:%S", tz = "UTC"), 
                         to = as.POSIXct("23:50:00", format = "%H:%M:%S", tz = "UTC"), 
                         by = "10 min"), "%H:%M:%S")

# Step 3: Loop through each file and process it
for (file in files) {
  # Read the CSV file
  data <- read.csv(paste0("./1_individual_AM_dataframes/",file))
  
  # Step 4: Dealing with midnight values that appear as N/As
  data$TIME <- ifelse(
    grepl("000000\\.WAV$", data$IN.FILE), "00:00:00",  # Explicitly set midnight time
    substr(data$IN.FILE, 10, 15)  # Extract HHMMSS from file name for non-midnight cases
  )
  
  # Fix the format of TIME to be consistent
  data$TIME <- gsub("^(\\d{2})(\\d{2})(\\d{2})$", "\\1:\\2:\\3", data$TIME)
  
  # Step 5: Create the Datetime column combining DATE and TIME
  data$Datetime <- as.POSIXct(
    paste(data$DATE, data$TIME),
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
  )
  
  # Step 7: Ensure that start_time and end_time are finite
  start_time <- min(data$Datetime, na.rm = TRUE)
  end_time <- max(data$Datetime, na.rm = TRUE)
  
  # Check for valid start_time and end_time
  if (is.finite(start_time) & is.finite(end_time)) {
    # Generate expected timestamps every 10 minutes
    expected_timestamps <- seq(from = start_time, to = end_time, by = "10 min")
  } else {
    # If invalid, print error message and skip the file
    cat("Invalid Datetime range in file:", file, "\n")
    next
  }
  
  # Step 8: Identify missing timestamps by comparing expected sequence with actual Datetime
  missing_timestamps <- setdiff(expected_timestamps, data$Datetime)
  
  # Convert missing timestamps to POSIXct format
  missing_timestamps <- as.POSIXct(missing_timestamps, origin = "1970-01-01", tz = "UTC")
  
  # Step 9: Create a data frame for missing timestamps with placeholder values (0)
  NA_value <- 0  # Define what value NA will take
  
  missing_data <- data.frame(
    Datetime = missing_timestamps,
    FOLDER = NA_value,  
    IN.FILE = NA_value,  
    CHANNEL = NA_value,  
    OFFSET = NA_value,  
    DURATION = NA_value,  
    DATE = as.Date(missing_timestamps),
    TIME = format(missing_timestamps, "%H:%M:%S"),
    HOUR = NA_value,  
    MEAN = NA_value,  
    SD = NA_value,  
    MODE = NA_value,  
    Q25 = NA_value,  
    KURT = NA_value,  
    NDSI = NA_value,  
    ACI = NA_value,  
    BI = NA_value,  
    BGN = NA_value,  
    SNR = NA_value,  
    ACT = NA_value,  
    LFC = NA_value,  
    MFC = NA_value,  
    HFC = NA_value,  
    CENT = NA_value,
    PC1 = NA_value,
    PC2 = NA_value,
    PC3 = NA_value,
    PC4 = NA_value,
    PC5 = NA_value,
    PC6 = NA_value,
    PC7 = NA_value,
    PC8 = NA_value,
    PC9 = NA_value,
    PC10 = NA_value,
    PC11 = NA_value,
    PC12 = NA_value,
    PC13 = NA_value,
    PC14 = NA_value,
    norm_PC1 = NA_value,
    norm_PC2 = NA_value,
    norm_PC3 = NA_value,
    Identifier = sub(".csv", "", file)  # Dynamically set Identifier to file name
  )
  
  # Retain only the specified columns in the 'data' data frame
  data <- data %>%
    select(SD, MODE, Q25, KURT, NDSI, ACI, BI, BGN, SNR, ACT, LFC, MFC, HFC, CENT, Datetime,
           PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,norm_PC1,norm_PC2,norm_PC3)
 
  missing_data <- missing_data %>%
    select(SD, MODE, Q25, KURT, NDSI, ACI, BI, BGN, SNR, ACT, LFC, MFC, HFC, CENT, Datetime,
           PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,norm_PC1,norm_PC2,norm_PC3)
  
  # Step 9: Scale the specified columns between 0 and 1 (Min-Max Scaling)
  cols_to_scale <- c("SD","MODE", "Q25", "KURT", "NDSI", "ACI", "BI", "BGN", "SNR", "ACT", "LFC", "MFC", "HFC", "CENT")
  
  # Apply scaling to both 'data' and 'missing_data' data frames
  data[cols_to_scale] <- lapply(data[cols_to_scale], function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  missing_data[cols_to_scale] <- lapply(missing_data[cols_to_scale], function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  
  # Replace all NAs in 'data' with NA_value
  data[is.na(data)] <- NA_value
  
  # Replace all NAs in 'missing_data' with NA_value
  missing_data[is.na(missing_data)] <- NA_value
  
  # Step 10: Combine the missing data with the original dataset
  data_filled <- rbind(data, missing_data)
  
  # Step 11: Sort the data by Datetime after combining
  data_filled <- data_filled[order(data_filled$Datetime), ]
  
  # Step 12: Ensure the date range is consistent for all processing
  data_filled <- data_filled %>%
    mutate(date = as.Date(Datetime)) %>%
    dplyr::filter(date %in% date_range)
  
  # Step 15: Export the processed `data_filled` data frame as a CSV file
  output_file <- paste0("./1b_individual_AM_dataframes_full/","Processed_", sub(".csv", "", file), ".csv")
  write.csv(data_filled, output_file, row.names = FALSE)
  
}
```

#### False-colour spectrograms: Plotting 24 hrs of acoustic data with noise reduction. 

![Ordesa Summer 24 hr spec](https://github.com/user-attachments/assets/48f0378d-b104-4dca-ae11-180b91eb9925)

```
#By Jack A. Greenhalgh. Feb, 2025.
#Instituto Pirenaico de Ecología (CSIC), Avenida Nuestra Señora de la Victoria, 22700, Jaca, Huesca, España.

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

#### Tapestry plots: Using one acoustic index 

```
#By Jack A. Greenhalgh1 & José Joaquín Lahoz Monfort2. Feb, 2025.
#1Department of Biology, McGill University, 1205 Dr Penfield Ave, Montreal, Quebec, H3A 1B1, Canada.
#2Instituto Pirenaico de Ecología (CSIC), Avenida Nuestra Señora de la Victoria, 22700, Jaca, Huesca, España.

# # Load required libraries
# library(ggplot2)
# library(dplyr)

# Step 1: Read all CSV files in the directory matching the pattern
files <- list.files(path = "INSERT YOUR SUBSETTED AND CLEANED DATA HERE")
print(files)

# Step 2: Define fixed date and time ranges for the plot
date_range <- seq(as.Date(start_date), as.Date(end_date), by = "day")
time_range <- format(seq(from = as.POSIXct("00:00:00", format = "%H:%M:%S", tz = "UTC"), 
                         to = as.POSIXct("23:50:00", format = "%H:%M:%S", tz = "UTC"), 
                         by = "10 min"), "%H:%M:%S")

# Step 3: Loop through each file and process it
for (file in files) {
  # Read the CSV file
  data <- read.csv(paste0("./1b_individual_AM_dataframes_full/",file))
  
  # Fix missing midnight values in the Datetime column
  data$Datetime <- as.character(data$Datetime)  # Convert to character for manipulation
  missing_midnight <- grepl("^\\d{4}-\\d{2}-\\d{2}$", data$Datetime)  # Identify rows missing time
  data$Datetime[missing_midnight] <- paste0(data$Datetime[missing_midnight], " 00:00:00")  # Add "00:00:00"
  
  # Convert back to POSIXct
  data$Datetime <- as.POSIXct(data$Datetime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # Step 4: Split Datetime into Date and Time for plotting
  data <- data %>% 
    mutate(
      date = as.Date(Datetime),  # Extract the date part
      time = format(Datetime, "%H:%M:%S")  # Extract the time part as character
    )
  
  # Convert `time` back to POSIXct for plotting
  data$time <- as.POSIXct(data$time, format = "%H:%M:%S", tz = "UTC")
  
  # LOOP FOR ALL ACOUSTIC INDICES:
  for (index in 1:13) {
      
      acousticindex <- names(data)[index]
    
      # Step 5: Generate HEX codes for each acoustic index
      data$HEX_codes <- rgb(
        red = data[,index],
        green = data[,index],
        blue = data[,index],
        maxColorValue = 1
      )
      
      # Step 6: Generate the plot with HEX color codes and standardized date/time ranges
      p <- ggplot(data, aes(x = time, y = date)) + 
        geom_tile(aes(fill = HEX_codes)) + 
        labs(x = "Time of day", y = "Date", title = paste(acousticindex," for", file)) + 
        scale_x_datetime(
          date_breaks = "2 hours", 
          date_labels = "%H:%M", 
          limits = c(as.POSIXct("00:00:00", format = "%H:%M:%S", tz = "UTC"), 
                     as.POSIXct("23:50:00", format = "%H:%M:%S", tz = "UTC"))
        ) + 
        scale_y_date(
          date_breaks = "1 month", 
          date_labels = "%b %Y", 
          limits = c(as.Date(start_date), as.Date(end_date))
        ) + 
        theme_bw() + 
        scale_fill_identity()
      
      # Step 7: Save the plot with a dynamic file name based on the current file
      ggsave(paste0(acousticindex,"_", sub(".csv", "", file), ".png"), 
             path = "./2_single_indiv_tapestries", plot = p, height = 6, width = 6)

  } # end of inner loop (acoustic indices)
  
} # end of outer loop (individual device files)
```

#### Tapestry plots: Using three acoustic indices

```
#By Jack A. Greenhalgh1 & José Joaquín Lahoz Monfort2. Feb, 2025.
#1Department of Biology, McGill University, 1205 Dr Penfield Ave, Montreal, Quebec, H3A 1B1, Canada.
#2Instituto Pirenaico de Ecología (CSIC), Avenida Nuestra Señora de la Victoria, 22700, Jaca, Huesca, España.

# library(ggplot2)
# library(dplyr)

# Step 1: Read all CSV files in the directory matching the pattern
files <- list.files(path = "INSERT YOUR SUBSETTED AND CLEANED DATA HERE")
print(files)

# Step 2: Define fixed date and time ranges for the plot
date_range <- seq(as.Date(start_date), as.Date(end_date), by = "day")
time_range <- format(seq(from = as.POSIXct("00:00:00", format = "%H:%M:%S", tz = "UTC"), 
                         to = as.POSIXct("23:50:00", format = "%H:%M:%S", tz = "UTC"), 
                         by = "10 min"), "%H:%M:%S")

# Step 3: Loop through each file and process it
for (file in files) {
  # Read the CSV file
  data <- read.csv(paste0("./1b_individual_AM_dataframes_full/",file))
  
  # Fix missing midnight values in the Datetime column
  data$Datetime <- as.character(data$Datetime)  # Convert to character for manipulation
  missing_midnight <- grepl("^\\d{4}-\\d{2}-\\d{2}$", data$Datetime)  # Identify rows missing time
  data$Datetime[missing_midnight] <- paste0(data$Datetime[missing_midnight], " 00:00:00")  # Add "00:00:00"
  
  # Convert back to POSIXct
  data$Datetime <- as.POSIXct(data$Datetime, format="%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # Split Datetime into Date and Time for plotting
  data <- data %>% 
    mutate(
      date = as.Date(Datetime),  # Extract the date part
      time = format(Datetime, "%H:%M:%S")  # Extract the time part as character
    )
  
  # Convert `time` back to POSIXct for plotting
  data$time <- as.POSIXct(data$time, format = "%H:%M:%S", tz = "UTC")
  
  # Part 3: Generate HEX codes using a column - CHOOSE INDICES HERE
  data$HEX_codes <- rgb(
    red = data$MODE,
    green = data$BI,
    blue = data$MFC,
    maxColorValue = 1
  )
  combo_name <- "MODE-BI-MFC"  # Change name accordingly
  
  # Step 5: Generate the plot with HEX color codes and standardized date/time ranges
  p <- ggplot(data, aes(x = time, y = date)) + 
    geom_tile(aes(fill = HEX_codes)) + 
    labs(x = "Time of day", y = "Date", title = paste(combo_name, " for", file)) + 
    scale_x_datetime(
      date_breaks = "2 hours", 
      date_labels = "%H:%M", 
      limits = c(as.POSIXct("00:00:00", format = "%H:%M:%S", tz = "UTC"), 
                 as.POSIXct("23:50:00", format = "%H:%M:%S", tz = "UTC"))
    ) + 
    scale_y_date(
      date_breaks = "1 month", 
      date_labels = "%b %Y", 
      limits = c(as.Date(start_date), as.Date(end_date))
    ) + 
    theme_bw() + 
    scale_fill_identity()
  
  # Step 6: Save the plot with a dynamic file name based on the current file
  ggsave(paste0(combo_name,"_", sub(".csv", "", file), ".png"), 
         path = "./3_combined_tapestries", plot = p, height = 6, width = 6)
}
```

### Seasonal drift code

```
#By Jack A. Greenhalgh1, Dec, 2025.
#1Department of Biology, McGill University, 1205 Dr Penfield Ave, Montreal, Quebec, H3A 1B1, Canada.
#### Part 1. Loading, cleaning, and scaling data ####

# Load packages
library(corrplot)
library(caret)
library(stringr)
library(lubridate)
library(dplyr)
library(hms)
library(ggplot2)

# Set working directory
setwd("C:/Users/Administrador/OneDrive - McGill University/Gault Data/Raw acoustic indices (terrestrial)")

# Load data
J001_May_July_2025 <- read.csv("J001_Maple_Beech_alpha_acoustic_indices_results_May_July_2025.csv")
J001_July_August_2025 <- read.csv("J001_Maple_Beech_alpha_acoustic_indices_results_July_August_2025.csv")
J002_May_July_2025 <- read.csv("J002_Oak_alpha_acoustic_indices_results_May_July_2025.csv")
J002_July_August_2025 <- read.csv("J002_Oak_alpha_acoustic_indices_results_July_August_2025.csv")
J003_May_July_2025 <- read.csv("J003_Lake_Shore_alpha_acoustic_indices_results_May_July_2025.csv")
J003_July_August_2025 <- read.csv("J003_Lake_Shore_alpha_acoustic_indices_results_July_August_2025.csv")
M006_May_July_2025 <- read.csv("M006_Beaver_Pond_alpha_acoustic_indices_results_May_July_2025.csv")
M006_July_August_2025 <- read.csv("M006_Beaver_Pond_alpha_acoustic_indices_results_July_August_2025.csv")
M007_May_July_2025 <- read.csv("M007_Wetland_alpha_acoustic_indices_results_May_July_2025.csv")  
M007_July_August_2025 <- read.csv("M007_Wetland_alpha_acoustic_indices_results_July_August_2025.csv")  
M008_May_July_2025 <- read.csv("M008_Oak_alpha_acoustic_indices_results_May_July_2025.csv")
M008_July_August_2025 <- read.csv("M008_Oak_alpha_acoustic_indices_results_July_August_2025.csv")
M009_May_July_2025 <- read.csv("M009_Maple_Beech_alpha_acoustic_indices_results_May_July_2025.csv")
M009_July_August_2025 <- read.csv("M009_Maple_Beech_alpha_acoustic_indices_results_July_August_2025.csv")
M010_May_July_2025 <- read.csv("M010_Beaver_Pond_alpha_acoustic_indices_results_May_July_2025.csv")  
M010_July_August_2025 <- read.csv("M010_Beaver_Pond_alpha_acoustic_indices_results_July_August_2025.csv")  

# Merge all data sets into one
merged_data <- rbind(
  J001_May_July_2025, J001_July_August_2025,
  J002_May_July_2025, J002_July_August_2025,
  J003_May_July_2025, J003_July_August_2025,
  M006_May_July_2025, M006_July_August_2025,
  M007_May_July_2025, M007_July_August_2025,
  M008_May_July_2025, M008_July_August_2025,
  M009_May_July_2025, M009_July_August_2025,
  M010_May_July_2025, M010_July_August_2025
)
head(merged_data)

# Extract filename
Filename <- merged_data$filename

# Subset numeric columns from 2 to 61
numeric_data <- merged_data[, 2:61]

# Compute correlation matrix
cor_matrix <- cor(numeric_data, use = "complete.obs")

# Plot correlation matrix
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.cex = 0.7, tl.col = "black", addCoef.col = "black", number.cex = 0.5)
print(cor_matrix)

# Find indices of highly correlated variables (threshold > 0.8)
high_corr_indices <- findCorrelation(cor_matrix, cutoff = 0.8, names = TRUE)

# Remove them from the dataset
filtered_data <- numeric_data[, !colnames(numeric_data) %in% high_corr_indices]

head(filtered_data)

# Apply z-transformation
z_scaled_data <- as.data.frame(scale(filtered_data))

# View the result
head(z_scaled_data)
summary(z_scaled_data)

# Add filename back
z_scaled_data <- cbind(Filename, z_scaled_data)
head(z_scaled_data)

# Extract site information and store in a new column called site
z_scaled_data$Site <- str_extract(z_scaled_data$Filename, "^.*?_.*?(?=_)")
head(z_scaled_data[c("Filename", "Site")])

#Extract datetime
z_scaled_data <- z_scaled_data %>%
  mutate(
    # Extract datetime string (e.g., "20250513_120000") using regex
    datetime_str = str_extract(Filename, "\\d{8}_\\d{6}"),
    
    # Parse into POSIXct format (YYYYMMDD_HHMMSS)
    Datetime = as.POSIXct(datetime_str, format = "%Y%m%d_%H%M%S")
  )

# Round timestamps to the nearest 10
z_scaled_data <- z_scaled_data %>%
  mutate(
    Datetime = round_date(Datetime, unit = "10 minutes")
  )

head(z_scaled_data)

# ===============================
# 🔹 Extract Site and Date from filename
# ===============================
z_scaled_data <- z_scaled_data %>%
  mutate(
    # Extract site code from the start of the filename (everything before first underscore)
    SiteCode = str_extract(Filename, "^[A-Z0-9]+"),
    
    # Extract date string from filename (YYYYMMDD)
    Date_str = str_extract(Filename, "\\d{8}"),
    
    # Convert date string to Date object
    Date = as.Date(Date_str, format = "%Y%m%d")
  )

# ===============================
# 🔹 Subset to unique days × site
# ===============================
unique_days_site <- z_scaled_data %>%
  group_by(SiteCode, Date) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup()

# ===============================
# 🔹 Check result
# ===============================
head(unique_days_site)

# ===============================
# 🔹 Select numeric features for PCA
# ===============================
numeric_cols <- unique_days_site %>%
  select(where(is.numeric)) %>%
  colnames()

# ===============================
# 🔹 Remove highly correlated features (optional)
# ===============================
numeric_data <- unique_days_site %>% select(all_of(numeric_cols))

# Remove constant or all-NA columns
numeric_data_clean <- numeric_data %>%
  select(where(~ !all(is.na(.)) & var(., na.rm = TRUE) > 0))

# Remove highly correlated variables
if(ncol(numeric_data_clean) > 1){
  corr_matrix <- cor(numeric_data_clean, use = "pairwise.complete.obs")
  high_corr <- findCorrelation(corr_matrix, cutoff = 0.9)
  numeric_data_uncor <- numeric_data_clean[, -high_corr]
} else {
  numeric_data_uncor <- numeric_data_clean
}

# ===============================
# 🔹 Run PCA on all numeric features
# ===============================
pca_res <- prcomp(numeric_data_uncor, scale. = TRUE)

# ===============================
# 🔹 Project each site × day into PCA space
# ===============================
pca_scores <- as.data.frame(pca_res$x[, 1:2]) # first two PCs
pca_scores$Site <- unique_days_site$SiteCode
pca_scores$Date <- unique_days_site$Date

# Rename columns for clarity
centroids <- pca_scores %>%
  rename(PC1 = PC1, PC2 = PC2)

# ===============================
# 🔹 Check results
# ===============================
head(centroids)

library(ggplot2)
library(grid)

# ===============================
# 🔹 Prepare arrow segments per site
# ===============================
arrow_segments <- centroids %>%
  arrange(Site, Date) %>%
  group_by(Site) %>%
  mutate(PC1_end = lead(PC1),
         PC2_end = lead(PC2)) %>%
  filter(!is.na(PC1_end) & !is.na(PC2_end)) %>%
  ungroup()

library(dplyr)
library(lubridate)

# Compute weekly centroids for all sites
centroids_weekly <- centroids %>%
  mutate(Week = floor_date(Date, unit = "week")) %>%  # Convert Date to week
  group_by(Site, Week) %>%                            # Group by Site and Week
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    Month = first(Month)  # Keep month for coloring
  ) %>%
  arrange(Site, Week) %>%  # Ensure ordering within each site
  ungroup() %>%
  group_by(Site) %>%        # Add a time index per site
  mutate(TimeIndex = row_number()) %>%
  ungroup()

# Get the unique sites
sites <- unique(centroids_weekly$Site)

# Loop through each site
for (site_id in sites) {
  
  # Filter data for the current site
  site_data <- centroids_weekly %>% filter(Site == site_id)
  
  # Prepare arrow segments
  arrow_segments <- site_data %>%
    mutate(
      PC1_end = lead(PC1),
      PC2_end = lead(PC2)
    ) %>%
    filter(!is.na(PC1_end) & !is.na(PC2_end))
  
  # Create the plot
  p <- ggplot(site_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Month), size = 7, alpha = 0.7) +
    geom_segment(
      data = arrow_segments,
      aes(x = PC1, y = PC2, xend = PC1_end, yend = PC2_end, color = Month),
      arrow = arrow(length = unit(0.3, "cm")),
      linetype = "solid"
    ) +
    geom_text(aes(label = TimeIndex), color = "white", size = 3, vjust = 0.5, hjust = 0.5) +
    theme_bw() +
    labs(
      title = paste("Site:", site_id),
      x = "PC1",
      y = "PC2",
      color = "Month"
    ) +
    scale_color_manual(
      values = colorRampPalette(c("darkred", "orange"))(length(unique(site_data$Month)))
    ) +
    coord_cartesian(xlim = pc1_range, ylim = pc2_range) +
    theme(legend.position = "bottom")
  
  # Save the plot as a PNG file named after the site
  ggsave(
    filename = paste0("plot_", site_id, ".png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}
```

### References

Villanueva-Rivera, L. J., & Pijanowski, B. C. (2018). Soundecology: Soundscape ecology. R package version 1.3.3. https://CRAN.R-project.org/package=soundecology
