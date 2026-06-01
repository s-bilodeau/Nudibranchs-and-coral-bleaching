#clear workspace
rm(list=ls())

#set working directory
setwd("") #where csv of data is

#load libraries
library(reshape)
library(lubridate)
library(ggpubr) #for boxplot
library(rstatix) #for repeated measures ANOVA
library(stringr) #for wrapping labels in plots

#read in data
PAM <- read.csv("NudibranchFeedingEffectsPAMReadings.csv") #PAM readings from nudibranch-exposed corals
flume_corals <- read.csv("BleachedControlCoralsPAMdata.csv") #contains Control corals not exposed to nudibranchs

#create a function for standard error
se <- function(x) sqrt(var(x)/length(x))


#################################### Compare PAM readings before and after nudibranch feeding #################################### 

#remove rows with insufficient baseline and those where the coral bleached apparently unrelated to nudibranch predation
PAM <- PAM[-c(3,5,9:14),]

#create columns to average the three readings from each timepoint/condition
PAM$BaselineAvg <- (PAM$Baseline1+PAM$Baseline2+PAM$Baseline3)/3
PAM$HealthyPostAvg <- (PAM$HealthyPost1+PAM$HealthyPost2+PAM$HealthyPost3)/3
PAM$DamagedPostAvg <- (PAM$DamagedPost1+PAM$DamagedPost2+PAM$DamagedPost3)/3

#melt dataframe to allow modeling that accounts for nudibranch ID
PAMsub <- PAM[,c(1:2, 15:17)]
mPAM <- melt(PAMsub)

#compare Baseline, Healthy, and Damaged
bxp <- ggboxplot(mPAM, x="variable", y="value", add="point")
bxp #quickly visualize data

#remove second (and third) fragments from a given nudibranch, reducing n but also eliminating nudibranch duplicates
PAMsingle <- PAMsub[-c(7,10,13,14,15,16,19),]
      mPAMsingle <- melt(PAMsingle)
      mPAMsingle$value <- mPAMsingle$value/1000 #rescale Y-axis values

#test assumptions for ANOVA, starting with outliers
mPAMsingle %>%
  group_by(variable) %>%
  identify_outliers(value) #no outliers

#test for normality
shapiro.test(mPAMsingle[which(mPAMsingle$variable=="BaselineAvg"),4])
shapiro.test(mPAMsingle[which(mPAMsingle$variable=="HealthyPostAvg"),4])
shapiro.test(mPAMsingle[which(mPAMsingle$variable=="DamagedPostAvg"),4]) #all normally distributed (assumption 1 for ANOVA)

ggqqplot(mPAMsingle, "value", facet.by="variable") #all are approximately normal

#repeated measures ANOVA
a2 <- anova_test(data=mPAMsingle, dv=value, wid=Coral, within=variable)
get_anova_table(a2)

pwc2 <- mPAMsingle %>%
  pairwise_t_test(
    value ~ variable, paired = TRUE,
    p.adjust.method = "holm"
  )
pwc2 #pairwise comparisons


#get mean and se for each grouping
mean(mPAMsingle[which(mPAMsingle$variable=="BaselineAvg"),4]) 
      se(mPAMsingle[which(mPAMsingle$variable=="BaselineAvg"),4]) 
mean(mPAMsingle[which(mPAMsingle$variable=="HealthyPostAvg"),4]) 
      se(mPAMsingle[which(mPAMsingle$variable=="HealthyPostAvg"),4]) 
mean(mPAMsingle[which(mPAMsingle$variable=="DamagedPostAvg"),4]) 
      se(mPAMsingle[which(mPAMsingle$variable=="DamagedPostAvg"),4]) 

      

#################################### Compare start and end PAM values for control corals #################################### 

#isolate only control corals from larger dataset
fcontrol <- flume_corals[which(flume_corals$Treatment=="Control"),c(1:2, 9:10)]

#compare start and end PAM readings for each of the 15 corals
t.test(x=fcontrol$HealthyAvg, y=fcontrol$BleachedAvg, paired=TRUE)


#################################### Compare controls to nudibranch damage #################################### 

#use multiple (adjusted) t-tests to compare DamagedPostAvg, HealthyPostAvg (healthy tissue of damaged corals) and BleachedAvg (healthy tissue of non-bleached control corals)
aPAM <- as.data.frame(cbind(fcontrol$BleachedAvg, c(PAMsingle$HealthyPostAvg, NA, NA, NA), c(PAMsingle$DamagedPostAvg, NA, NA, NA)))

t1 <- t.test(aPAM$V1, aPAM$V2, paired=FALSE, var.equal=FALSE) #compare Control with HealthyPost
      t1
t2 <- t.test(aPAM$V1, aPAM$V3, paired=FALSE, var.equal=FALSE) #compare Control with DamagedPost
      t2 
t3 <- t.test(aPAM$V2, aPAM$V3, paired=FALSE, var.equal=FALSE) #compare HealthyPost with DamagedPost
      t3 

# Extract p-values
p_values <- c(t1$p.value, t2$p.value, t3$p.value)

# Perform Holm correction
p_adjusted <- p.adjust(p_values, method = "holm")
      p_adjusted 

#convert aPAM into a plottable format
abPAM <- matrix(nrow=15, ncol=0)
      abPAM <- as.data.frame(abPAM)
      
      abPAM$Control <- aPAM$V1/1000
      abPAM$Indirect <- aPAM$V2/1000
      abPAM$Direct <- aPAM$V3/1000
      
      mabPAM <- melt(abPAM)
      
      mabPAMlabs <- c("Control", "Indirect Feeding Effects", "Direct Feeding Effects") #create x-axis group labels
      
#create a boxplot with control corals and separated PAM readings for direct and indirect nudibranch damage
PAMboxplot <- ggboxplot(mabPAM, x="variable", y="value", color="variable", palette = c("#045568", "#117F94", "#1DA9BF"),
                        size=0.8, add = c("mean_se", "jitter")) +
      labs(y="Effective Quantum Yield (Y)", x="") + #add titles to both axes
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), #remove background and grid
        axis.line = element_line(colour = "black"), #add axis lines
        axis.title.x=element_text(margin=margin(b=45), vjust=-5, size=25), #move x-axis title lower and increase bottom whitespace
        axis.title.y=element_text(margin=margin(r=30), size=25),
        axis.text=element_text(color="black", size=18),
        plot.margin = unit(c(20,20,40,20), "pt"),
        legend.position="none") + #add space in outer margins of plot
      coord_cartesian(ylim = c(0, 0.8), #this allows the annotations for the n labels to go below the plot area
                  clip = 'off') +
      scale_x_discrete(labels=str_wrap(mabPAMlabs, width=20)) + #replace default x-axis labels with our custom labels, using stringr function to wrap long text
      annotate(geom="text", x = c("Control", "Indirect", "Direct"), y=c(0.7, 0.68, 0.44), label=c("a", "b", "c"), size=6) + #add letters for significance groupings
      annotate(geom="text", x = c("Direct"), y=c(0.7), label=c("p≤0.03"), size=7) + #add a p-value
      annotate(geom="text", x = c("Control", "Indirect", "Direct"), y=c(-0.14, -0.14, -0.14), label=c("n=15", "n=12", "n=12"), size=6) #add n for each group

PAMboxplot



#calculate mean days for nudibranch damaage to reach 70% (requires converting dates using lubridate)
PAM$BaselineDate_5PM <- ymd(PAM$BaselineDate_5PM)
      PAM$DamagedDate_5PM <- ymd(PAM$DamagedDate_5PM)
      PAM$TotalDays <- (PAM$DamagedDate_5PM - PAM$BaselineDate_5PM) #creates difftime class object (reported in days)

      mean(PAM[-c(7,10,13,14,15,16,19), 18]) #subset removes the rows with multiple corals from the same nudibranch, making this equivalent to PAMsingle

