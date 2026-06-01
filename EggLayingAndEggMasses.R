#clear workspace
rm(list=ls())

#set working directory
setwd("") #where csv of data is

#load libraries
library(dplyr)

#read in data
eggs <- read.csv("SlugEggsPositioningCorals.csv")

egg_mass <- read.csv("EggMassMeasurements.csv")

#create a function for standard error
se <- function(x) sqrt(var(x)/length(x))


########################## Egg-Laying Analysis ########################## 

#subset to include all nudibranchs that are Medium or Large
ML <- eggs[c(which(eggs$Size=="Large"|eggs$Size=="Medium"),6:8),] 
#this subset includes the "Unknown" nudibranchs from 3-22, 1 of which we know is Medium (7) and 1 of which is Large (8)
#and the one nudibranch sizable from photos: 6 (L)

ML <- ML[-6,] #remove Individual 17 (died after 1 day, incomplete data)

#average eggs laid/day over time they were alive
ML <- ML %>%
      mutate_at(vars(EggsDay2, EggsDay3, EggsDay4, EggsDay5, EggsDay6, EggsDay7, EggsDay8, EggsDay9, EggsDay10), as.numeric) #convert all egg count columns to numeric
      ML$AvgPerDay <- rowMeans(ML[,7:11], na.rm=TRUE) 
      #currently includes the first five days only 
      #(Day6 results being the 24 hours starting on Day5)

#average means of all possible egg-laying nudibranchs (M and L)
mean(ML$AvgPerDay) #currently includes only the first 5 days in which all egg-laying nudibranchs laid at least once
      sd(ML$AvgPerDay)
      se(ML$AvgPerDay)

########################## Eggs Per Mass ########################## 
      
#calculate mean eggs per egg mass for the egg masses examined in the lab
mean(egg_mass[which(egg_mass$SizeClass=="Large"),6])
      sd(egg_mass[which(egg_mass$SizeClass=="Large"),6])
      se(egg_mass[which(egg_mass$SizeClass=="Large"),6])

mean(egg_mass[which(egg_mass$SizeClass=="Small"),6])
      sd(egg_mass[which(egg_mass$SizeClass=="Small"),6])
      se(egg_mass[which(egg_mass$SizeClass=="Small"),6])


