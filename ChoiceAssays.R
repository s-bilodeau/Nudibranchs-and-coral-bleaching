#clear workspace
rm(list=ls())

#set working directory
setwd("") #where csv of data is

#load libraries
library(reshape)
library(ggplot2)
library(rstatix)

#read in data
choice <- read.csv("AllTenelliaChoiceAssays2021.csv") #this csv contains compiled data from the timepoints in each 2021 assay where most individuals made a choice
choice2 <- read.csv("All2025ChoiceAssays.csv") #this csv contains compiled data from the timepoints in each 2025 assay where most individuals made a choice

#create separate dataframes for each separate assay
corals <- choice[c(2:5, 7:10, 14:20, 41:45), c(2,5,6,7)] 
corals2 <- choice[c(21:40), c(2,10:11)]
porites <- choice[c(1:20), c(2:4)]
bleach <- choice[c(21:37), c(2, 8:9)]
positions <- choice2[36:57, 2:7] 
cyclodius <- choice2[61:72, c(2,9:10)]
corallio <- choice2[1:15, c(2, 11:12)]
drupella <- choice2[1:15, c(2, 13:14)]
chlorodi <- choice2[c(2,4:12, 14:15), c(2, 15:17)]
jbleach <- choice2[16:35, c(2, 18:19)]
      

########################## Coral Choice Assays ########################## 
      
#Cochran's Q takes a formula of the form a ~ b | c, where a is the outcome variable name (0 or 1),
#b is the within-subjects factor variables (coral type),
#and c (factor) is the column name containing individuals/subjects identifier (UID)

#assess coral choices for the overnight timepoint; this is the assay reported in the appendix of the manuscript (Fig A1)

#reformat data as binary for test, starting with a mostly-empty dataframe to fill
nightchoice <- as.data.frame(corals$UID) #copy nudibranch identifiers
      names(nightchoice)[names(nightchoice)=="corals$UID"] <- "UID" #rename column
      nightchoice$RUS <- NA #add empty columns for each coral
      nightchoice$POR <- NA
      nightchoice$DAM <- NA
      nightchoice$PAV <- NA
      nightchoice$ACR <- NA

#loop through rows and for each column representing a coral in the new dataframe, 
#check if that coral is the one specified at the overnight timepoint (1) or not (0)
for(i in 1:nrow(corals)){
      nightchoice[i,2] <- ifelse(corals[i,4]=="RUS",1,0); 
      nightchoice[i,3] <- ifelse(corals[i,4]=="POR",1,0);
      nightchoice[i,4] <- ifelse(corals[i,4]=="DAM",1,0);
      nightchoice[i,5] <- ifelse(corals[i,4]=="PAV",1,0);
      nightchoice[i,6] <- ifelse(corals[i,4]=="ACR",1,0)
}

mnight <- melt(nightchoice) #melt data into a better format for the test

mnight$value <- as.factor(mnight$value) #ensure that our binary outcomes are coded as a factor

mnight$UID <- as.factor(mnight$UID) #ensure that our individual nudibranchs are coded as a factor

cochran_qtest(mnight, value~variable|UID) #run the test

pairwise_mcnemar_test(mnight, value~variable|UID, p.adjust.method="holm") #pairwise comparisons between groups 
#using the Holm adjustment method for multiple planned pairwise comparisons because it has higher power than Bonferroni
#and we are making a lot of comparisons here!

#Repeat for the overnight timepoint (21hrs) of coral choice 2; this is the assay reported in the main mansucript (Fig 1)
#reformat data as binary for test, starting with a mostly-empty dataframe to fill
nightchoice2 <- as.data.frame(corals2$UID) #copy nudibranch identifiers
      names(nightchoice2)[names(nightchoice2)=="corals2$UID"] <- "UID" #rename column
      nightchoice2$RUS <- NA #add empty columns for each coral
      nightchoice2$MON <- NA
      nightchoice2$VER <- NA
      nightchoice2$ACR <- NA

#loop through rows and for each column representing a coral in the new dataframe, 
#check if that coral is the one specified at the overnight timepoint (1) or not (0)
for(i in 1:nrow(corals2)){
      nightchoice2[i,2] <- ifelse(corals2[i,3]=="RUS",1,0); 
      nightchoice2[i,3] <- ifelse(corals2[i,3]=="MON",1,0);
      nightchoice2[i,4] <- ifelse(corals2[i,3]=="VER",1,0);
      nightchoice2[i,5] <- ifelse(corals2[i,3]=="ACR",1,0)
}
      
#drop NA row where the nudibranch made no choice
nightchoice2 <- na.omit(nightchoice2)

mnight2 <- melt(nightchoice2) #melt data into a better format for the test

mnight2$value <- as.factor(mnight2$value) #ensure that our binary outcomes are coded as a factor

mnight2$UID <- as.factor(mnight2$UID) #ensure that our individual nudibranchs are coded as a factor

cochran_qtest(mnight2, value~variable|UID) #run the test

pairwise_mcnemar_test(mnight2, value~variable|UID, p.adjust.method="holm") #pairwise comparisons between groups 
#using the Holm adjustment method for multiple planned pairwise comparisons because it has higher power than Bonferroni
#and we are making a lot of comparisons here!

#format data for plotting coral choices
coral_bar <- as.data.frame(matrix(nrow=5)) #create counts for bar chart
      coral_bar$name <- c("P. rus", "P. lobata", "P. damicornis", "A. pulchra", "P. cactus") #n=20
      coral_bar$value <- c(nrow(corals[which(corals$Coral_choice_night_14hrs=="RUS"),]), nrow(corals[which(corals$Coral_choice_night_14hrs=="POR"),]), 
                           nrow(corals[which(corals$Coral_choice_night_14hrs=="DAM"),]), nrow(corals[which(corals$Coral_choice_night_14hrs=="ACR"),]),
                           nrow(corals[which(corals$Coral_choice_night_14hrs=="PAV"),]))

#create a bar chart showing one bar for each choice option for corals (nighttime; Fig A1)
cplot <- ggplot(coral_bar) +
      geom_bar(aes(x=reorder(name, -value), y=value), stat="identity", alpha=0.7, fill=c("#067F9C", "#067F9C","#1DA9BF", "#1DA9BF","#1DA9BF")) + #plot bars for both treatments, reorder by counts
      labs(y="Number of Nudibranchs on Coral", x="Coral Species") + #add titles to both axes
      scale_y_continuous(limits = c(0,20), expand = c(0, 0)) + #set y-axis limits and remove space between bars and x-axis
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), #remove background and grid
            plot.margin=margin(t=20, b=25, l=10), #add padding at top and sides of plot so axis doesn't run off page
            axis.line = element_line(colour = "black"), #add axis lines
            axis.title.x=element_text(margin=margin(b=25), vjust=-8, size=25), #move x-axis title lower and increase bottom whitespace
            axis.title.y=element_text(margin=margin(r=25), size=25),
            axis.text=element_text(color="black", size=15)) + #move y-axis title farther left and increase leftmost whitespace
      annotate(geom="text", x = c("P. rus", "P. lobata", "P. damicornis", "A. pulchra", "P. cactus"), y=(coral_bar$value+2), label=c("a", "ab", "ab", "b", "b"), size=7) #add letters for significance groupings
cplot

#create a bar chart showing one bar for each choice option for the second set of corals (nighttime; Fig 1)
coral_bar2 <- as.data.frame(matrix(nrow=4)) #create counts for bar chart
      coral_bar2$name <- c("Porites", "Montipora", "Pocillopora", "Acropora") #n=19
      coral_bar2$value <- c(nrow(corals2[which(corals2$Coral_choice_2_21hrs=="RUS"),]), nrow(corals2[which(corals2$Coral_choice_2_21hrs=="MON"),]), 
                            nrow(corals2[which(corals2$Coral_choice_2_21hrs=="VER"),]), nrow(corals2[which(corals2$Coral_choice_2_21hrs=="ACR"),]))

c2plot <- ggplot(coral_bar2) +
      geom_bar(aes(x=reorder(name, -value), y=value), stat="identity", alpha=0.7, fill=c("#067F9C","#1DA9BF", "#1DA9BF","#1DA9BF")) + #plot bars for both treatments, reorder by counts
      labs(y="Number of Nudibranchs on Coral", x="Coral Genus") + #add titles to both axes
      scale_y_continuous(limits = c(0,20), expand = c(0, 0)) + #set y-axis limits and remove space between bars and x-axis
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), #remove background and grid
            plot.margin=margin(t=20, b=25, l=10), #add padding at top and sides of plot so axis doesn't run off page
            axis.line = element_line(colour = "black"), #add axis lines
            axis.title.x=element_text(margin=margin(b=25), vjust=-8, size=25), #move x-axis title lower and increase bottom whitespace
            axis.title.y=element_text(margin=margin(r=25), size=25),
            axis.text=element_text(color="black", size=15)) + #move y-axis title farther left and increase leftmost whitespace
      annotate(geom="text", x = c("Porites", "Montipora", "Pocillopora", "Acropora"), y=(coral_bar2$value+2), label=c("a", "b", "b", "b"), size=7) #add letters for significance groupings
c2plot


#use binomial tests to compare outcomes for two-choice assays (Porites, bleaching, and juvenile bleaching)
binom.test(c(nrow(porites[which(porites$Porites_choice_night_240min=="RUS"),]), nrow(porites[which(porites$Porites_choice_night_240min=="POR"),])), 
           p=0.5, alternative="two.sided", conf.level = 0.95) #adult Porites only 

binom.test(c(nrow(bleach[which(bleach$Bleached_choice_17hrs=="B"),]), nrow(bleach[which(bleach$Bleached_choice_17hrs=="H"),])), 
           p=0.5, alternative="two.sided", conf.level = 0.95) #adult bleached vs. healthy P. rus 

binom.test(c(nrow(jbleach[which(jbleach$Juvenile_bleaching_17hrs=="B"),]), nrow(jbleach[which(jbleach$Juvenile_bleaching_17hrs=="H"),])), 
           p=0.5, alternative="two.sided", conf.level = 0.95) #juvenile bleached vs. healthy P. rus 

#compare egg laying on RUS and POR (7 on RUS, 2 on POR)
binom.test(c(7,2), p=0.5, alternative="two.sided", conf.level=0.95)

#plot adult and juvenile choices from the bleaching assay
bleach_barC <- as.data.frame(matrix(nrow=2)) #create dataframe for bar chart with bleaching counts
      bleach_barC$name <- c("Bleached", "Healthy") #n=17
      bleach_barC$value <- c(nrow(bleach[which(bleach$Bleached_choice_17hrs=="B"),]), nrow(bleach[which(bleach$Bleached_choice_17hrs=="H"),]))

#create a bar chart showing adult bleaching choice by counts
blplotC <- ggplot(bleach_barC) +
      geom_bar(aes(x=reorder(name, -value), y=value), stat="identity", fill=c("#61C3D2","#067F9C"), alpha=0.7) + #plot counts, bleached first
      labs(y="Number Chosen", x="") + #add y-axis title
    scale_y_continuous(limits = c(0,12), expand = c(0, 0)) + #set y-axis limits and remove space between bars and x-axis
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), #remove background and grid
          plot.margin=margin(t=20, b=5, l=10), #add padding at top and sides of plot so axis doesn't run off page
          axis.line = element_line(colour = "black"), #add axis lines
          axis.title.x=element_text(margin=margin(b=50), vjust=-8, size=25), #move x-axis title lower and increase bottom whitespace
          axis.title.y=element_text(margin=margin(r=25), size=25),
          axis.text=element_text(color="black", size=20)) + #move y-axis title farther left and increase leftmost whitespace
    annotate(geom="text", x = c(x=1.5), y=10.5, label="p=0.63", size=8) + #add letters for significance groupings
    annotate(geom="text", x = c(x=0.9), y=11.8, label="A. Adults", size=10) #add plot label at top
blplotC

#create corresponding bar chart for juvenile bleaching choices
jbleach_bar25 <- as.data.frame(matrix(nrow=3)) #create counts for bar chart
      jbleach_bar25$name <- c("Healthy", "Bleached", "No Choice") #n=20, 4 NAs
      jbleach_bar25$value <- c(nrow(choice2[which(choice2$Juvenile_bleaching_17hrs=="H"),]), 
                               nrow(choice2[which(choice2$Juvenile_bleaching_17hrs=="B"),]), 
                               nrow(choice2[which(is.na(choice2$Juvenile_bleaching_17hrs)==TRUE),]))

#create a bar chart showing juvenile belaching choice by counts
jblplot25 <- ggplot(jbleach_bar25[1:2,]) + #exclude NA row
      geom_bar(aes(x=reorder(name, -value), y=value), stat="identity", fill=c("#067F9C", "#61C3D2"), alpha=0.7) + #plot counts, bleached first
      labs(y="Number Chosen", x="") + #add y-axis title
      scale_y_continuous(limits = c(0,12), expand = c(0, 0)) + #set y-axis limits and remove space between bars and x-axis
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), #remove background and grid
            plot.margin=margin(t=20, b=5, l=10), #add padding at top and sides of plot so axis doesn't run off page
            axis.line = element_line(colour = "black"), #add axis lines
            axis.title.x=element_text(margin=margin(b=50), vjust=-8, size=25), #move x-axis title lower and increase bottom whitespace
            axis.title.y=element_text(margin=margin(r=25), size=25),
            axis.text=element_text(color="black", size=20)) + #move y-axis title farther left and increase leftmost whitespace
      annotate(geom="text", x = c(x=1.5), y=8.5, label="p=1", size=8) + #add letters for significance groupings
      annotate(geom="text", x = c(x=0.9), y=11.8, label="B. Juveniles", size=10) #add plot label at top (use x=2 for right aligned label)
jblplot25


########################## Crab and Snail Choice Assays ########################## 

#use binomial test to compare outcomes for two-choice assays
binom.test(c(nrow(chlorodi[which(chlorodi$Chlorodiella_16hrs=="C"),]), nrow(chlorodi[which(chlorodi$Chlorodiella_16hrs=="T"),])), 
           p=0.5, alternative="two.sided", conf.level = 0.95) #Chlorodiella crab choice assay with corals

binom.test(c(nrow(cyclodius[which(cyclodius$Cyclodius_30min=="C"),]), nrow(cyclodius[which(cyclodius$Cyclodius_30min=="T"),])), 
           p=0.5, alternative="two.sided", conf.level = 0.95) #Cyclodius crab choice assay with corals

binom.test(c(nrow(corallio[which(corallio$Coralliophila_16hrs=="C"),]), nrow(corallio[which(corallio$Coralliophila_16hrs=="T"),])), 
           p=0.5, alternative="two.sided", conf.level = 0.95) #Coralliophila snail choice assay with corals

binom.test(c(nrow(drupella[which(drupella$Drupella_16hrs=="C"),]), nrow(drupella[which(drupella$Drupella_16hrs=="T"),])), 
           p=0.5, alternative="two.sided", conf.level = 0.95) #Drupella snail choice assay with corals


########################## Nudibranch Overnight Positioning ########################## 

#use Cochran's Q to compare positions (sheltered vs. exposed) at daylight and nighttime timepoints
#positions "Under" and "Bottom" have been classified as "Sheltered" while "Side", "Mid", and "Top" are "Exposed"

#reformat data as binary for test, starting with an empty dataframe to fill
simple_pos <- as.data.frame(matrix(nrow=66, ncol=0))
      simple_pos$UID <- rep(1:22,3) #add repeated individual identifiers
      simple_pos$Position <- NA #add a column for position
      simple_pos$Time <- c(rep("afternoon",22), rep("midnight",22), rep("predawn",22)) #add a column for time of day

#merge location signifiers for the sheltered category
positions[positions=="Under"] <- "Bottom"

#loop through the three timepoints and classify each nudibranch as sheltered (0) or exposed (1)
for(i in 1:22){
      simple_pos[i,2] <- ifelse(positions[i,2]=="Bottom",0,1);
      simple_pos[i+22,2] <- ifelse(positions[i,3]=="Bottom",0,1);
      simple_pos[i+44,2] <- ifelse(positions[i,4]=="Bottom",0,1)
}

#some smaller nudibranchs can become trapped in the surface tension or otherwise separated from the coral 
#for mobility reasons rather than making a choice to leave or stay away; 
#therefore, non-association with a coral may not be equivalent to non-sheltering behavior on a coral,
#so we only evaluate exposed/sheltered positions for nudibranchs that show a consistent coral association

#reference the original positions dataframe to remove individuals that did not stay associated with the coral the whole time
simple_pos <- simple_pos[-c(1:2,4,9:10,13:14,16:18,20:21,23:24,26,31:32,35:36,38:40,42:43,45:46,48,53:54,57:58,60:62,64:65),]

simple_pos$Position <- as.factor(simple_pos$Position) #ensure that our binary outcomes are coded as a factor

simple_pos$UID <- as.factor(simple_pos$UID) #ensure that our individual nudibranchs are coded as a factor

#given potential dockside activity near the water tables after dark, 
#we will compare the afternoon positioning with predawn positioning 
#(after nudibranchs have been undisturbed for multiple hours)

#compare afternoon and predawn
pos1<- simple_pos[-c(11:20),] #remove intermediate midnight timepoint
      cochran_qtest(pos1, Position~Time|UID) #run the test


