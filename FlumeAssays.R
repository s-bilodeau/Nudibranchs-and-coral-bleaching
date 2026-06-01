#clear workspace
rm(list=ls())

#set working directory
setwd("") #where csv of data is

#load libraries
library(reshape) #for melt function
library(ggpubr) #for barplot
library(stringr) #for wrapping labels in plots
library(diptest) #for analyzing bimodal results

#read in data
flume <- read.csv("TenelliaNightPoritesFlumeAssays.csv")
flume2 <- read.csv("TenelliaNightBleachingFlumeAssays.csv")


#create a function for standard error
se <- function(x) sqrt(var(x)/length(x))

########################## Porites tracking in the flume (2021) ########################## 

#create column for percent not in cue
flume$PercentContBoth <- ((60-flume$TotalCueBoth-flume$TotalMid)/60)*100 

#collapse dataframe to exclude "duplicate" rows for each trial (i.e., only one row per trial)
row_odd <- seq_len(nrow(flume)) %% 2
      flume_odd <- flume[row_odd == 1, ] 

#Chi square test (as used in Dixson et al. 2008; compares total trials with preference)
#(i.e. >50% time) and total without, assumption of 1:1 ratio

#create dataset and run Chi square test (counts of preference and not, total to 20)
ChiCounts <- c(sum(flume_odd$PercentCueBoth>50), sum(flume_odd$PercentCueBoth<=50)) 
      chisq.test(ChiCounts) 
      #default performs a goodness-of-fit test (x is a one-dimensional contingency table)
      #the null is that the population probabilities are equal (since p is not specified)

#Wilcoxon sign-rank test (allows amount of time to be taken into account with rankings)
flume_odd$PreferenceIndex <- (2*(flume_odd$TotalCueBoth)-60+flume_odd$TotalMid)/(60-flume_odd$TotalMid) #create Preference Index between -1 and 1 (0=no preference)
      #using formula (CueTime - ControlTime)/TotalTime 
      #adjusted to exclude middle times, i.e. reflective of adjusted percentages
      
    wilcox.test(flume_odd$PreferenceIndex, mu=0) #run single-sample Wilcoxon signed-rank test to compare preference index against null hypothesis of 0 (no preference)
    wilcox.test(flume_odd$PercentCueAdj, mu=50) #run single-sample Wilcoxon signed-rank test to compare un-adjusted cue percentage against null hypothesis of 50% (no preference)
    #both tests give equivalent results, confirming we did the Preference Index correctly
    #both give warnings because of ties and zero values (therefore the p-value is only approximated based on the normal distribution)

    wilcox.test(flume_odd[-c(3,5,9,18), 43], mu=0) #re-running the test with the tied values (3 100s) and the one zero value removed (therefore, n=16) still yields a significant result


########################## Bleached vs. healthy Porites in the flume (2025) ########################## 

#create column for percent not in cue
flume2$PercentContBoth <- ((60-flume2$TotalCueBoth-flume2$TotalMid)/60)*100 

#collapse dataframe to exclude "duplicate" rows for each trial (i.e., only one row per trial)  
row_odd25 <- seq_len(nrow(flume2)) %% 2
      flume_odd25 <- flume2[row_odd25 == 1, ] 

#create dataset and run Chi square test (counts of preference and not, total to 20)
ChiCounts25 <- c(sum(flume_odd25$PercentCueBoth>50), sum(flume_odd25$PercentCueBoth<=50))
      chisq.test(ChiCounts25) 
      #default performs a goodness-of-fit test (x is a one-dimensional contingency table)
      #and the null is that the population probabilities are equal (since p is not specified)

#Wilcoxon sign-rank test (allows amount of time to be taken into account with rankings)
flume_odd25$PreferenceIndex <- (2*(flume_odd25$TotalCueBoth)-60+flume_odd25$TotalMid)/(60-flume_odd25$TotalMid) #create Preference Index between -1 and 1 (0=no preference)
      #using formula (CueTime - ControlTime)/TotalTime 
      #adjusted to exclude middle times, i.e. reflective of adjusted percentages

      wilcox.test(flume_odd25$PreferenceIndex, mu=0) #run single-sample Wilcoxon signed-rank test to compare preference index against null hypothesis of 0 (no preference)
      wilcox.test(flume_odd25$PercentCueAdj, mu=50) #run single-sample Wilcoxon signed-rank test to compare adjusted cue percentage against null hypothesis of 50% (no preference)
      #both tests are equivalent
      #both give warnings because of ties and zero values (therefore the p-value is only approximated based on the normal distribution)

      wilcox.test(flume_odd25[-c(3,4,6,7,11,13,14,15,16,18,20), 40], mu=0) #re-running the test with the tied values and zeros removed (therefore, n=9)
      #this is likely not a fair test to run having artificially eliminated the strongest preferences in both directions

#run Hardigan's Dip Test to determine if the responses are bimodal (or at least not unimodal)
dip.test(flume_odd25$PreferenceIndex)


###################################### Create combined flume figure ###################################################

#combine and format data from 2021 and 2025 assays
mflumeexp <- melt(flume_odd[,40]) #melt 2021 data for combined graphing (using adjusted percentages)
  mflumeexp$variable <- c(rep("P.rus", 20))
  flumeplot <- flume_odd25[,40] #subset adjusted percentages for 2025, "Cue" = Bleached
  mflumeplot <- melt(flumeplot)
  mflumeplot$variable <- c(rep("Bleached", 20))
  combflume <- rbind(mflumeexp, mflumeplot) #combine flume data from 2021 and 2025

#add mean and se values
combflume$mean <- c(rep(mean(combflume[which(combflume$variable=="P.rus"),1]),20), rep(mean(combflume[which(combflume$variable=="Bleached"),1]),20))
combflume$se <- c(rep(se(combflume[which(combflume$variable=="P.rus"),1]),20), rep(se(combflume[which(combflume$variable=="Bleached"),1]),20))

#created (wrappable) labels for the x axis
flumelabs <- str_wrap(c("P. rus", "Bleached"), width=20) 

#create p-value, n, and axis labels using geom_text
annotation <- data.frame( 
      x = c("P.rus", "Bleached"),
      y = c(25,25),
      label = c("p=0.001", "p=0.73")
      )
annotation2 <- data.frame( 
      x = c("P.rus","Bleached"),
      y = c(21,21),
      label = c("n=20", "n=20")
      )
annotationbottom <- data.frame( 
      x = c(1,2),
      y = c(-12.5,-12.5),
      label = c("(vs. seawater control)", "(vs. healthy coral)")
      )

#combined barplot
combflumeplot <- ggbarplot(
        combflume, x = "variable", y = "value", 
        add = c("mean", "jitter"), 
        fill= "variable", palette = c("#067F9C", "#61C3D2"),
        position = position_dodge(0.8),
        alpha=0.7,
        size=0 #setting size to 0 should remove the black box around the bars, but in Windows this only works consistently when saving with another graphics program (e.g. Cairo)
      ) +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), #remove background and grid
            axis.line = element_line(colour = "black"), #add axis lines
            axis.title.x=element_text(margin=margin(b=45), vjust=-5, size=25), #move x-axis title lower and increase bottom whitespace
            axis.title.y=element_text(margin=margin(r=30), size=25),
            axis.text=element_text(color="black", size=17),
            plot.margin = unit(c(20,20,10,20), "pt"), #add space in outer margins of plot
            legend.position = "none") +
      coord_cartesian(ylim = c(0,100), #this allows the annotations for the n labels to go below the plot area
                      clip = 'off') +
      labs(y="Percent Time Spent in Cue", x="") + #add titles to both axes
      scale_x_discrete(labels=str_wrap(flumelabs, width=20)) + #replace default x-axis labels with our custom labels, using stringr function to wrap long text
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.3, colour="black", alpha=0.1, size=0.8) + #add standard error bars
      geom_hline(yintercept=50, linetype="dashed", color="dark gray", size=1) + #add a dashed line denoting 50%
      geom_text(data=annotation, aes( x=x, y=y, label=label),                  
                  color="black", 
                  size=5) +
      geom_text(data=annotation2, aes( x=x, y=y, label=label),                  
                  color="black", 
                  size=4.5) +
      geom_text(data=annotationbottom, aes( x=x, y=y, label=label),                  
                  color="black", 
                  size=5.5)
combflumeplot

