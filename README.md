# Pyrenean Institute of Ecology: Ecoacoustic Protocols  
### By Jack A. Greenhalgh, September 2023

### Contents 

1. Part 1: Survey design
2. Part 2: Calculating acoustic indices
3. Part 3: Data handling and manipulation
4. Part 4: Data visulisation 
5. Part 5: Dealing with spatial replication

### Part 1: Ecoacoustic Survey Design

### Part 2: Calculating Acoustic Indices 

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

### References

Villanueva-Rivera, L. J., & Pijanowski, B. C. (2018). Soundecology: Soundscape ecology. R package version 1.3.3. https://CRAN.R-project.org/package=soundecology
