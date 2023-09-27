library(cluster)
library(dplyr)
library(ggplot2)

####Subsample####

#Load data
Data <- read.table("SeasonalDataGroups10.txt", header=T, sep="\t")
attach(Data)
head(Data)

#Ensure the order factor levels for month and season 
Data$Month <- factor(Data$Month, levels=c("January", "Feburary", "March", "April", "May", 
                                          "June","July", "August", "September", "October", 
                                          "November", "December"))
Data$Season <- factor(Data$Season, levels=c("Winter", "Autumn", "Summer", "Spring"))

#display factor levels for region
levels(Data$Month)

#Plot Rose Plot
jpeg("Rose Plot.jpeg", width = 7, height = 7, units = 'in', res = 300)
ggplot(data=Data,aes(x=Month,y=Season,fill=Value))+ 
  geom_tile(colour="black",size=0.1)+ 
  scale_fill_gradientn(name="Bioacoustic index", colours=c("lightgreen","darkgreen"))+
  coord_polar()+xlab("")+ylab("") + theme_minimal() 
dev.off()


####Full data set####

# Read your data from the file
BI <- read.table("Bioacoustic Index all months.txt", header = TRUE, sep = "\t")

# Define your groups and row numbers (-1 row number from Excel)
group_info <- data.frame(Group = c("Spring", "Summer", "Autumn", "Winter"),
                         StartRow = c(1, 121819, 230832, 372640),
                         EndRow = c(121818, 230831, 372639, 465395))

# Ensure that group_info is sorted by StartRow
group_info <- group_info[order(group_info$StartRow), ]

# Use mutate to add the new column
BI <- BI %>%
  mutate(Group = factor(findInterval(rownames(BI), group_info$StartRow), 
                        labels = group_info$Group))

#Ensure the order factor levels for month and season for plotting 
BI$Month <- factor(BI$Month, levels=c("January", "February", "March", "April", "May", 
                                          "June","July", "August", "September", "October", 
                                          "November", "December"))
BI$Group <- factor(BI$Group, levels=c("Winter", "Autumn", "Summer", "Spring"))

#display factor levels for region
levels(BI$Month)

#Plot Rose Plot
jpeg("Rose Plot.jpeg", width = 7, height = 7, units = 'in', res = 300)
ggplot(data=Data,aes(x=Month,y=Season,fill=Value))+ 
  geom_tile(colour="black",size=0.1)+ 
  scale_fill_gradientn(name="Bioacoustic index", colours=c("lightgreen","darkgreen"))+
  coord_polar()+xlab("")+ylab("") + theme_minimal() 
dev.off()

