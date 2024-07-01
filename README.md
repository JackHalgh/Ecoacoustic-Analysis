# Pyrenean Institute of Ecology: Ecoacoustic Protocols  
![Title image](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/d1275870-6d41-46f6-8c4b-a1dd37745752)
### By [Jack A. Greenhalgh](https://www.jack-greenhalgh.com/), September 2023.     

### Contents 

- [Part 1: Ecoacoustic equipment and survey design](#ecoacoustic-equipment-and-survey-design)
- [Part 2: Acoustic indices](#acoustic-indices)
- [Part 3: Data handling and manipulation](#data-handling-and-manipulation)
- [Part 4: Data visulisation](#data-visulisation)
- [Part 5: Dealing with spatial replication](#dealinig-with-spatial-replication)

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
install.packages("soundecology")
library(soundecology)

multiple_sounds("name_of_your_respository", resultfile = "name_of_your_output.csv",
                soundindex = "acoustic_complexity")

multiple_sounds("name_of_your_respository", resultfile = "name_of_your_output.csv",
                soundindex = "acoustic_diversity")

multiple_sounds("name_of_your_respository", resultfile = "name_of_your_output.csv",
                soundindex = "acoustic_evenness")

multiple_sounds("name_of_your_respository", resultfile = "name_of_your_output.csv",
                soundindex = "bioacoustic_index")

multiple_sounds("name_of_your_respository", resultfile = "name_of_your_output.csv",
                soundindex = "H")

multiple_sounds("name_of_your_respository", resultfile = "name_of_your_output.csv",
                soundindex = "ndsi")
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


#### Tapestry plots: Step 1. 

```
# Remove all objects in the global environment
set.seed(123)
rm(list = ls())
setwd("C:/Users/Administrador/Downloads/R/Ordesa/Data")

library(readxl)

# Read the Excel files
disk1 <- read_excel("acousticindex1.xlsx")
disk2 <- read_excel("acousticindex2.xlsx")

# Append disk2 to the bottom of disk1
All_Data <- rbind(disk1, disk2)

# Remove AEI and ADI 
All_Data <- subset(All_Data, select = -c(SKEW, KURT, CENT))

# Sub-setting to keep only rows with DURATION >= 60
All_Data_filtered <- All_Data[All_Data$DURATION >= 60, ]

# Sub-setting to keep only columns with acoustic indices
All_Data_trimmed <- All_Data_filtered[, -(1:8)]

# Define a function to scale data between 0 and 1
scale_0_1 <- function(x) {
  return ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# Convert relevant columns to numeric
numeric_columns <- sapply(All_Data_trimmed, is.character)
All_Data_trimmed[numeric_columns] <- lapply(All_Data_trimmed[numeric_columns], as.numeric)

# Scale each numeric column between 0 and 1
All_Data_scaled <- as.data.frame(lapply(All_Data_trimmed, function(x) 
  if(is.numeric(x)) scale_0_1(x) else x))

# Display the first few rows of the scaled data frame
head(All_Data_scaled)

All_Data_scaled <- read.csv("All_Data_Scaled.csv")

# Identify and remove highly correlated variables

#review correlation matrix and remove highly correlated variables 

library(corrplot)

# Calculate correlation matrix
correlation_matrix <- cor(All_Data_scaled)

correlation_matrix

# Plot correlation matrix
corrplot(correlation_matrix, method = "square")

# Remove highly correlated variables
library(caret)

# Set a threshold 
threshold <- 0.8 

# Find highly correlated variables
highly_correlated <- findCorrelation(correlation_matrix, cutoff = threshold)

# Get names of highly correlated variables
highly_correlated_names <- colnames(correlation_matrix)[highly_correlated]

# Print the names of highly correlated variables
print(highly_correlated_names)

# Find indices of highly correlated variables in the dataframe
highly_correlated_indices <- match(highly_correlated_names, colnames(All_Data_scaled))

# Remove highly correlated variables from the entire dataframe
All_Data_Cleaned <- All_Data_scaled[, -highly_correlated_indices]

# Print the cleaned dataframe
head(All_Data_Cleaned)

#### Step 2. Perform Principal Component Analysis ####

indices_pca <- prcomp(All_Data_Cleaned, scale. = F)
indices_pca$PC1 <- indices_pca$x[,1]
indices_pca$PC2 <- indices_pca$x[,2]
indices_pca$PC3 <- indices_pca$x[,3]
indices_pca$PC4 <- indices_pca$x[,4]
indices_pca$PC5 <- indices_pca$x[,5]
indices_pca$PC6 <- indices_pca$x[,6]
indices_pca$PC7 <- indices_pca$x[,7]

pca_coef <- cbind(indices_pca$PC1, indices_pca$PC2,
                  indices_pca$PC3, indices_pca$PC4,
                  indices_pca$PC5, indices_pca$PC6,
                  indices_pca$PC7)
rm(indices_pca)

coef_min_max <- pca_coef[,1:3]

# Scale the PCA coefficients between 0 and 1 so they can be 
# mapped to red, green and blue channels.

# Define the normalise function
normalise <- function(x, min_val, max_val) {
  (x - min_val) / (max_val - min_val)
}

coef_min_max_norm <- coef_min_max
min.values <- NULL
max.values <- NULL
for (i in 1:3) {
  min <- unname(quantile(pca_coef[,i], probs = 0.0, na.rm = TRUE))
  max <- unname(quantile(pca_coef[,i], probs = 1.0, na.rm = TRUE))
  min.values <- c(min.values, min)
  max.values <- c(max.values, max)
  coef_min_max_norm[,i]  <- normalise(coef_min_max[,i], min, max)
}

head(coef_min_max_norm)

# Combine PCAs to data frame 
colnames(coef_min_max_norm) <- c("PC1", "PC2", "PC3")
All_Data_Cleaned <- cbind(All_Data_Cleaned, coef_min_max_norm)

#### Step 4. Generate HEX colour codes from RGB #### 

library(ggplot2)
library(reshape2)

Hex_Colour_Codes <- rgb(All_Data_Cleaned$PC1, All_Data_Cleaned$PC2, All_Data_Cleaned$PC3)
head(Hex_Colour_Codes)

# Combine HEX colour codes, recorder and datetime to the data frame 
All_Data_HEX  <- cbind(All_Data_Cleaned,Hex_Colour_Codes)
All_Data_Final  <- cbind(All_Data_HEX,All_Data_filtered$`IN FILE`, All_Data_filtered$FOLDER, All_Data_filtered$DATE)

head(All_Data_Final)

write.csv(All_Data_Final, "All Data Final.csv")

#### Step 5. Sub-setting by recorder #### 

# Find unique values in the 'FOLDER' column
unique_folders <- unique(All_Data_Final$`All_Data_filtered$FOLDER`)

# Create a list to store data frames for each subset
subset_data_frames <- list()

# Loop through each unique folder name
for (folder in unique_folders) {
  # Subset the original data frame based on the current folder
  subset_data <- All_Data_Final[All_Data_Final$`All_Data_filtered$FOLDER` == folder, ]
  
  # Store the subset data frame in the list with folder name as key
  subset_data_frames[[folder]] <- subset_data
}

# Function to filter data frames based on folder names containing specific patterns
filter_data_frames <- function(list_of_data_frames, pattern) {
  filtered_data_frames <- list()
  for (folder_name in names(list_of_data_frames)) {
    if (grepl(pattern, folder_name)) {
      filtered_data_frames[[folder_name]] <- list_of_data_frames[[folder_name]]
    }
  }
  return(filtered_data_frames)
}

# Group based on Audiomoth number (adjust 1:36 based on number of recorders)
grouped_data_frames <- list()
patterns <- paste0("J", formatC(1:36, width = 3, format = "d", flag = "0"))

for (pattern in patterns) {
  grouped_data_frames[[pattern]] <- filter_data_frames(subset_data_frames, pattern)
}

# Create a list to store single data frames for each pattern
single_data_frames <- list()

# Loop through each pattern
for (pattern in names(grouped_data_frames)) {
  # Get the list of data frames for the current pattern
  data_frames_list <- grouped_data_frames[[pattern]]
  
  # Concatenate all data frames within the pattern
  if (length(data_frames_list) > 0) {
    combined_df <- do.call(rbind, data_frames_list)
    
    # Store the combined data frame in the new list
    single_data_frames[[pattern]] <- combined_df
  }
}

# Export all single data frames as .csv files
for (pattern in names(single_data_frames)) {
  # Define the file name for the current pattern
  file_name <- paste0(pattern, "_data.csv")
  
  # Write the data frame to a .csv file
  write.csv(single_data_frames[[pattern]], file = file_name, row.names = FALSE)
}
```

#### Tapestry plots: Step 2. 
```
# Set seed, locale, and working directory
set.seed(123)
rm(list = ls())
Sys.setlocale("LC_TIME", "English")
setwd("C:/Users/Administrador/Downloads/R/Ordesa/Data")

# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(scales)

# Get list of all J00x files in the directory matching J001 to J036
files <- list.files(pattern = "J0(0[1-9]|[1-2][0-9]|3[0-6])_data_wo_CENT.csv")

# Function to handle missing dates
fill_missing_dates <- function(data, expected_dates) {
  all_dates <- data.frame(All_Data_filtered.DATE = expected_dates)
  merged_data <- merge(all_dates, data, by = "All_Data_filtered.DATE", all.x = TRUE)
  merged_data[is.na(merged_data)] <- 0
  return(merged_data)
}

# Loop through each file
for (file in files) {
  # Load data
  data <- read.csv(file)
  
  # Convert DATE column to Date format
  data$All_Data_filtered.DATE <- as.Date(data$All_Data_filtered.DATE, format = "%Y-%m-%d")
  
  # Sort the dates in ascending order
  data <- data[order(data$All_Data_filtered.DATE), , drop = FALSE]
  
  # Calculate the expected sequence of dates
  expected_dates <- seq(min(data$All_Data_filtered.DATE), max(data$All_Data_filtered.DATE), by = "day")
  
  # Identify missing dates
  missing_dates <- setdiff(expected_dates, data$All_Data_filtered.DATE)
  missing_dates_as_dates <- as.Date(missing_dates, origin = "1970-01-01")
  formatted_missing_dates <- format(missing_dates_as_dates, "%Y/%m/%d")
  
  # Check for missing dates and handle accordingly
  if (length(formatted_missing_dates) == 0) {
    print(paste("***No missing dates in time series for", file, "*** No missing dates."))
  } else {
    print(paste("***Missing dates in time series for", file, "*** Filling missing dates with zero values."))
    data <- fill_missing_dates(data, expected_dates)
  }
  
  # Convert IN.FILE to datetime format
  data$datetime <- as.POSIXct(strptime(substr(data$All_Data_filtered..IN.FILE., 1, 14), format = "%Y%m%d_%H%M%S"), tz = "UTC")
  
  # Ensure no missing datetime values
  data <- data[!is.na(data$datetime),]
  
  # Sort the dataframe by datetime
  data <- data[order(data$datetime),]
  
  # Calculate time differences between consecutive rows
  time_diff <- c(0, diff(data$datetime))
  
  # Identify where there are gaps greater than 10 minutes
  gap_indices <- which(time_diff > (10 * 60))
  
  # Create a new dataframe to store the gaps
  gaps_df <- data.frame(All_Data_filtered..IN.FILE. = character(),
                        PC1 = numeric(),
                        PC2 = numeric(),
                        PC3 = numeric(),
                        datetime = as.POSIXct(character(), tz = "UTC"),
                        stringsAsFactors = FALSE)
  
  # Iterate over the dataframe and fill gaps
  for (i in 1:(length(data$datetime) - 1)) {
    current_time <- data$datetime[i]
    next_time <- data$datetime[i + 1]
    time_diff <- as.numeric(difftime(next_time, current_time, units = "mins"))
    
    if (!is.na(time_diff) && time_diff > 10) {
      num_intervals <- ceiling(time_diff / 10)
      for (j in 1:(num_intervals - 1)) {
        new_row <- data.frame(IN.FILE = paste("Gap_", i, "_", j, sep = ""),
                              PC1 = 0,
                              PC2 = 0,
                              PC3 = 0,
                              datetime = current_time + (j * (10 * 60)))
        gaps_df <- rbind(gaps_df, new_row)
      }
    }
  }
  
  # Sub-set main data frame to include only file name, PCAs, and datetime
  subset <- data[, c("All_Data_filtered..IN.FILE.", "PC1", "PC2", "PC3", "datetime")]
  colnames(subset)[colnames(subset) == "All_Data_filtered..IN.FILE."] <- "IN.FILE"
  
  # Change column names in gaps_df to match
  colnames(gaps_df)[colnames(gaps_df) == "X1"] <- "PC1"
  colnames(gaps_df)[colnames(gaps_df) == "X2"] <- "PC2"
  colnames(gaps_df)[colnames(gaps_df) == "X3"] <- "PC3"
  
  # Merge gaps_df with the original dataframe
  rgb_data <- rbind(subset, gaps_df)
  
  # Sort the dataframe by datetime
  imputed_zeros_RGB <- rgb_data[order(rgb_data$datetime),]
  
  # Remove duplicated rows based on the datetime column
  imputed_zeros_RGB <- imputed_zeros_RGB[!duplicated(imputed_zeros_RGB$datetime),]
  
  # Generate HEX colour codes from scaled PCAs
  combined_colors <- rgb(imputed_zeros_RGB$PC3, imputed_zeros_RGB$PC2, imputed_zeros_RGB$PC1)
  colour_codes <- as.data.frame(combined_colors)
  combined_df <- cbind(imputed_zeros_RGB$datetime, colour_codes)
  colnames(combined_df)[colnames(combined_df) == "combined_colors"] <- "HEX colour code"
  colnames(combined_df)[colnames(combined_df) == "imputed_zeros_RGB$datetime"] <- "datetime"
  
  # Ensure datetime is properly parsed as POSIXct
  combined_df$datetime <- as.POSIXct(combined_df$datetime, format = "%Y-%m-%d %H:%M:%S")
  
  # Split the datetime column into date and time
  combined_df <- combined_df %>%
    mutate(
      date = as.Date(datetime),
      time = format(datetime, "%H:%M:%S")
    ) %>%
    select(-datetime)
  
  # Convert date and time to proper data types
  combined_df$date <- as.Date(combined_df$date, format = "%d/%m/%Y")
  combined_df$time <- as.POSIXct(combined_df$time, format = "%H:%M")
  
  # Plot using ggplot2
  p <- ggplot(combined_df, aes(x = time, y = date)) +
    geom_tile(aes(fill = `HEX colour code`)) +
    labs(x = "Time of day", y = "Date") +
    scale_x_datetime(date_breaks = "2 hours", date_labels = "%H:%M") +
    scale_y_date(date_breaks = "1 month", date_labels = "%b %Y") +
    theme_bw() +
    scale_fill_identity()
  
  # Save the plot
  plot_name <- paste0(sub("_data_wo_CENT.csv", "", file), " tapestry plot BGR (321).jpeg")
  ggsave(plot_name, plot = p, width = 8, height = 8, dpi = 300)
}
```

#### Calculate hourly means from values of acoustic indices

1) Using base R:

```
Data$Hour <- as.factor(Data$Hour)

mean_values <- aggregate(. ~ Hour, data = Data, mean)

print(mean_values)

write.csv(mean_values, "Hourly_means.csv")
```

2) Using for loops:

```
SiteCode_Index <- read.table("Site_code_acoustic_index.txt", header = T, sep = "\t")
attach(SiteCode_Index)

Index_Means <- as.numeric()

for (i in 1: put number of hours sample here) {
  SiteCode_Index_Subset <- SiteCode_Index[SiteCode_Index$Hour==i,]
  Index_Means[i] <- mean(SiteCode_Index_Subset$Index)
}

SiteCode_HourlyMeans <- data.frame("Hour"= c(1:put number of hours sample here), "Index value"= Index_Means)
attach(SiteCode_HourlyMeans)

write.csv(SiteCode_HourlyMeans, "SiteCode_HourlyMeans.csv")
```


### Data visulisation 

#### Visualising daily and seasonal variation 

There are many ways to visualise soundscape variation over time in R Studio. Here, we are going to take a look at some of the most popular methods. 

#### 1. Generalized additive models

Here, we will use a generalized additive model to predict daily acoustic variation of a pond soundscape from hourly means of the Bioacoustic Index. You can follow along with this example by using the HourlyMeans dataset in this repository. 

```
#Load hourly means data
HourlyMeans <- read.table("HourlyMeans.txt", sep = "\t", header=T)
head(HourlyMeans)
  Hour   BioMean
1    0 33.169258
2    1 28.401034
3    2 23.632810
4    3 18.864585
5    4 14.096361
6    5  9.328136

attach(HourlyMeans)

#Check variables are numeric
HourlyMeans$Hour <- as.numeric(HourlyMeans$Hour)
class(BioMean)
class(Hour)

#Run GAM. You may have to change the k (knots) value to avoid over-fitting the model
library(mgcv)
library(nlme)
gam_1 <- gam(BioMean ~ s(Hour, k=10), data = HourlyMeans, 
                  method = "REML")

#Plot GAM
plot(gam_1, residuals = T, rug = F, pch = 1, cex = 1,
     shade = T, shade.col = "lightblue", seWithMean = TRUE)

#Print results 
summary(gam_1)

#Check model
gam.check(gam_1, old.style=F)
qq.gam(gam_1,rep=100, pch = 1, cex = 1)
hist(gam_1$residuals)

#Export model predictions 
model_p <- predict(gam_1, type = "response")
write.csv(model_p, "gam_1.csv")

#Import formatted model predictions
gam_1_response <- read.table("gam_1.txt", header = T, sep = "\t")

#Plot model predictions using ggplot2
library(ggplot2)
jpeg("Daily GAM.jpeg", width = 7, height = 7, units = 'in', res = 300)
ggplot() +
  geom_line(data = gam_1_response, aes(x=Hour, y=Fit), 
            color="darkblue") +  
  theme_classic() + ylab("Bioacoustic index") + xlab("Hour of day")
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18,20,22,24))
dev.off()

```
![Daily GAM](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/3e77e87a-0526-43ad-abb5-573ba6aa0338)

#### 2. Rose plots

Next, we'll look at the use of rose plots to visualise seasonal soundscape data. 

```
head(Data)
  Season     Month      Value
1 Winter Janurary  6.896953
2 Winter Janurary  2.858804
3 Winter Janurary  7.197239
4 Winter Janurary  4.475334
5 Winter Janurary  11.957710
6 Winter Janurary  12.339196

jpeg("Rose Plot.jpeg", width = 7, height = 7, units = 'in', res = 300)
ggplot(data=Data,aes(x=Month,y=Season,fill=Value))+ 
  geom_tile(colour="black",size=0.1)+ 
  scale_fill_gradientn(name="Bioacoustic index", colours=c("lightgreen","darkgreen"))+
  coord_polar()+xlab("")+ylab("") + theme_minimal() 
dev.off()
```
![Annual rose plot](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/68d2ceb1-f558-4ff3-bf7f-81a383b039d9)

Of course, either method can be used to visualise daily, weekly, monthly, or seaosnal data. These are just examples and you should find the best method for your data. 

### Automated monitoring of bird populations using BirdNET 

BirdNET is a free online resource developed by the Cornell Lab of Ornithology that offers automated analysis of multiple sound files to identify brid calla and can be downloaded here: https://birdnet.cornell.edu/ 

During the analysis, BirdNET produces a csv result file for every minuite of recorded audio. Therefore, it is necessary to combine the output files for further analysis of the data. 

The following code combines multiple csv output files into a single csv file. 

```
# Set the directory where your BirdNET output CSV files are located
setwd("DEFINE PATH")

# List all the CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")

# Initialize an empty dataframe to store the combined data
combined_data <- data.frame()

# Loop through each CSV file and read it into a temporary dataframe
for (file in csv_files) {
  birdNET_data <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  
  # Combine the data from the temporary dataframe with the combined_data
  combined_data <- rbind(combined_data, birdNET_data)
}

# Reset the working directory to its original value (optional)
setwd("DEFINE PATH")

# Write the combined data to a single CSV file
write.csv(combined_data, file = "Combined_Results.csv", row.names = FALSE)
```

Next, we can subset the collated data to investigate all of the detections of a single species. 

```
library(lessR)

####AudioMoth####

combined_data <- read.table("AudioMoth_Combined_Results.txt", header = T, sep = "\t")

#Subset data frame by species and calculate median and IQRs

Barn_Owl <- combined_data[.(Common.name=="Barn Owl"), .(Scientific.name:Confidence)]
Barn_Owl_IQR <- quantile(Barn_Owl$Confidence)

Black_bellied_Plover <- combined_data[.(Common.name=="Black-bellied Plover"), .(Scientific.name:Confidence)]
Black_bellied_Plover_IQR <- quantile(Black_bellied_Plover$Confidence)

....
```


### Dealinig with spatial replication

### References

Villanueva-Rivera, L. J., & Pijanowski, B. C. (2018). Soundecology: Soundscape ecology. R package version 1.3.3. https://CRAN.R-project.org/package=soundecology
