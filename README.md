# Pyrenean Institute of Ecology: Ecoacoustic Protocols  
![Title image](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/d1275870-6d41-46f6-8c4b-a1dd37745752)
### By [Jack A. Greenhalgh](https://www.jack-greenhalgh.com/), September 2023.     

### Contents 

- [Part 1: Ecoacoustic Survey Design](#ecoacoustic-equipment-and-survey-design)
- [Part 2: Calculating acoustic indices](#calculating-acoustic-indices)
- [Part 3: Data handling and manipulation](#data-handling-and-manipulation)
- [Part 4: Data visulisation](#data-visulisation)
- [Part 5: Dealing with spatial replication](#dealinig-with-spatial-replication)

### Ecoacoustic equipment and survey design

#### Acoustic recorders

Here is a list of the most widely used acoustic recorders for long-term soundscape monitoring in a variety of environments. 

| Acoustic recorder | Realm(s)           | Cost ($) | Link                                                                         |
|-------------------|--------------------|----------|------------------------------------------------------------------------------|
| AudioMoth         | Terrestrial        | 89.99    | https://www.openacousticdevices.info/audiomoth                               |
| HydroMoth         | Marine, Freshwater | 135.00   | https://groupgets.com/manufacturers/open-acoustic-devices/products/hydromoth |
| Song Meter Mini   | Terrestrial        | 576.00   | https://www.veldshop.nl/en/song-meter-mini.html                              |

![IMG-1223](https://github.com/JackHalgh/Ecoacoustic-Analysis/assets/74665965/63e1240d-2659-4dc5-bd51-a619a99f1dc0)
AudioMoth

#### Transect 

#### Grid 

### Calculating acoustic Indices 

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

### Data handling and manipulation

#### Add a column for a catagorical variable to large datasets using dplyr

In this case, we are going to add a new column to our dataset called 'Season'. Adding this variable to a dataset with >500,000 rows using Excel is challenging, but it's simple in R Studio. 

```
install.packages("dplyr")
install.packages("cluster")
library(dplyr)
library(cluster)

# Load data 
Data <- read.table("your_dataframe.txt", header = TRUE, sep = "\t")

# Define your groups and row numbers
group_info <- data.frame(Season = c("Spring", "Summer", "Autumn", "Winter"),
                         StartRow = c(1, 121819, 230832, 372640),
                         EndRow = c(121818, 230831, 372639, 465395))

# Ensure that group_info is sorted by StartRow
group_info <- group_info[order(group_info$StartRow), ]

# Use mutate to add the new column to the dataframe 
Data <- Data %>%
  mutate(Season = factor(findInterval(rownames(Data), group_info$StartRow), 
                        labels = group_info$Group))
```		


### Data visulisation 

#### Visualising daily and seasonal variation 

There are many ways to visualise soundscape variation over time in R Studio. Here, we are going to take a look at some of the most popular methods. 


#### 1. Rose plots

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


### Dealinig with spatial replication

### References

Villanueva-Rivera, L. J., & Pijanowski, B. C. (2018). Soundecology: Soundscape ecology. R package version 1.3.3. https://CRAN.R-project.org/package=soundecology
