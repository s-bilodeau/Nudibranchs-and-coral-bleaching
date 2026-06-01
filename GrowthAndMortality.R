#clear workspace
rm(list=ls())

#set working directory
setwd("") #where csv of data is

#load libraries
library(ggplot2) #for plotting
library(ggpubr) #for boxplot
library(dunn.test) #for Kruskal-Wallis post-hoc comparisons

#read in data
growth <- read.csv("TenelliaGrowthData.csv")


########################## Compare Mortality Across Treatments ########################## 

#Use a Chi square to compare mortality across groups
chisq.test(c(nrow(growth[which(growth$Treatment=="lobata_healthy"&growth$Mortality==1),]), nrow(growth[which(growth$Treatment=="pulchra_healthy"&growth$Mortality==1),]),
             nrow(growth[which(growth$Treatment=="rus_healthy"&growth$Mortality==1),]),nrow(growth[which(growth$Treatment=="rus_bleached"&growth$Mortality==1),])))
              

#compare just Porites treatments' mortality (since mortality is 20% for both, this test is not really necessary)
binom.test(c(nrow(growth[which(growth$Treatment=="lobata_healthy"&growth$Mortality==1),]), nrow(growth[which(growth$Treatment=="rus_healthy"&growth$Mortality==1),])))
           

#compare bleached P. rus and healthy A. pulchra mortality (we are teasing out which groups from the Chi Square differ)
binom.test(c(nrow(growth[which(growth$Treatment=="rus_bleached"&growth$Mortality==1),]), nrow(growth[which(growth$Treatment=="pulchra_healthy"&growth$Mortality==1),])))


#Neither ANOVA nor Kruskal-Wallis is intended for binary responses (like mortality), 
#but we can run an ANOVA and Tukey post-hoc just to see how different treatments group
#and compare that to the apparent groupings from the Chi square and binomial tests above
anova1 <- aov(Mortality~Treatment, growth) 
      summary(anova1)
      TukeyHSD(anova1)

      #check ANOVA assumptions
      hist(growth$Mortality)
      qqnorm(anova1$residuals)
      qqline(anova1$residuals)
      shapiro.test(growth$Mortality) #data is DEFINITELY not normal, 
      #but we are not actually relying on this test to identify significant differences, 
      #simply to confirm the groupings identified above


########################## Compare Growth Across Treatments ########################## 

#calculate change in length as a new column
growth$Length_Change_mm <- growth$End_Length_mm - growth$Start_Length_mm #calculate change in length

#use an ANOVA or Kruskal-Wallis to compare length change (Growth) across groups
anova2 <- aov(Length_Change_mm~Treatment, growth)
      summary(anova2)
      TukeyHSD(anova2) #this version excludes those who died (their growth is NA), so it is biased towards survivors in the bleached treatment

#create a new dataframe which counts death (NA growth) as growth "0" (includes measurements for dead individuals where possible)
growth_all <- growth       
      growth_all$Length_Change_mm[is.na(growth_all$Length_Change_mm)] <- 0 
      
#create ANOVA with new length that includes deaths
anova3 <- aov(Length_Change_mm~Treatment, growth_all)
      summary(anova3)
      TukeyHSD(anova3)
      
      #check ANOVA assumptions
      hist(growth_all$Length_Change_mm)
      qqnorm(anova3$residuals)
      qqline(anova3$residuals)
      shapiro.test(growth_all$Length_Change_mm) #data is DEFINITELY not normal; 
      #while this assumption can sometimes be violated with large data sets,
      #Kruskal-Wallis is likely more appropriate here
      
#run a Kruskal-Wallis test and a Dunn post-hoc test to get groupings for the Kruskal-Wallis results
kruskal.test(x=list(growth_all$Length_Change_mm[which(growth_all$Treatment=="lobata_healthy")], growth_all$Length_Change_mm[which(growth_all$Treatment=="rus_healthy")],
                    growth_all$Length_Change_mm[which(growth_all$Treatment=="rus_bleached")], growth_all$Length_Change_mm[which(growth_all$Treatment=="pulchra_healthy")]))

dunn.test(x=list(growth_all$Length_Change_mm[which(growth_all$Treatment=="lobata_healthy")], growth_all$Length_Change_mm[which(growth_all$Treatment=="rus_healthy")],
                 growth_all$Length_Change_mm[which(growth_all$Treatment=="rus_bleached")], growth_all$Length_Change_mm[which(growth_all$Treatment=="pulchra_healthy")]))
#lobata = 1, rus = 2, bleached rus = 3, pulchra = 4 (for group comparisons)

  
#since not all nudibranchs were the same starting length, we actually want to compare percent change in length across groups
    
#create a column for percent change
growth_all$PercentChange <- ((growth_all$End_Length_mm - growth_all$Start_Length_mm)/abs(growth_all$Start_Length_mm))*100

#replace NA values with zeros for percent change
for(i in 1:length(growth_all$PercentChange)) {
  
  if (is.na(growth_all$PercentChange[i])==TRUE) growth_all$PercentChange[i] <- 0
  else growth_all$PercentChange[i] = growth_all$PercentChange[i]
}


#run a Kruskal-Wallis test and a Dunn post-hoc test to get groupings for the Kruskal-Wallis results
kruskal.test(x=list(growth_all$PercentChange[which(growth_all$Treatment=="lobata_healthy")], growth_all$PercentChange[which(growth_all$Treatment=="rus_healthy")],
                    growth_all$PercentChange[which(growth_all$Treatment=="rus_bleached")], growth_all$PercentChange[which(growth_all$Treatment=="pulchra_healthy")]))

dunn.test(x=list(growth_all$PercentChange[which(growth_all$Treatment=="lobata_healthy")], growth_all$PercentChange[which(growth_all$Treatment=="rus_healthy")],
                 growth_all$PercentChange[which(growth_all$Treatment=="rus_bleached")], growth_all$PercentChange[which(growth_all$Treatment=="pulchra_healthy")]))
#lobata = 1, rus = 2, bleached rus = 3, pulchra = 4 (for group comparisons)


########################## Plot Mortality and Growth ##########################

#note that plot labels are not italicized in this format and may overlap depending on the dimensions of your plotting window

#compile data for a bar chart
live_bar <- as.data.frame(matrix(nrow=4)) 
      live_bar$name <- c("P. lobata", "P. rus", "P. rus bleached", "A. pulchra")
      live_bar$value <- c((sum(growth[which(growth$Treatment=="lobata_healthy"),7])/20)*100, (sum(growth[which(growth$Treatment=="rus_healthy"),7])/20)*100, 
                    (sum(growth[which(growth$Treatment=="rus_bleached"),7])/20)*100, (sum(growth[which(growth$Treatment=="pulchra_healthy"),7])/20)*100)

#create a bar chart showing one bar for percent mortality in each treatment (P. rus bleached, P. rus unbleached, P. lobata, and A. pulchra)
mplot <- ggplot(live_bar) +
  geom_bar(aes(x=reorder(name, value), y=value), stat="identity", fill=c("#067F9C", "#045568", "#1DA9BF", "#848484"), alpha=0.7) + #plot means for both treatments, reorder bars to place cue first
  labs(y="Percent Mortality", x="Treatment") + #add titles to both axes
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) + #set y-axis limits and remove space between bars and x-axis
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), #remove background and grid
        axis.line = element_line(colour = "black"), #add axis lines
        axis.title.x=element_text(margin=margin(b=25), vjust=-7, size=25), #move x-axis title lower and increase bottom whitespace
        axis.title.y=element_text(margin=margin(r=25), size=25),
        axis.text=element_text(size=20, colour="black"),
        plot.margin = unit(c(20,20,20,20), "pt")) + #add space in outer margins of plot) + #move y-axis title farther left and increase leftmost whitespace
  geom_text(aes(x = name, y = value, label = c("a","a","b","b")), vjust = -.8, size=7) + #add letters for significance groupings
  geom_text(x = 0.65, y=98, label=c("A"), size=9)#add plot label A
  
mplot


#create a box plot for percent change in length (using the growth_all dataset above)
plbox <- ggboxplot(growth_all, x="Treatment", y="PercentChange", order=c("lobata_healthy", "rus_healthy", "rus_bleached", "pulchra_healthy"), 
                  color="Treatment", 
                  palette=c( "#067F9C", "#045568", "#1DA9BF", "#848484"), size=0.8, add = c("mean_se", "jitter")) + #set color and show mean, se, and points
  labs( x = "Treatment", y = "Percent Length Change", color = "") + #add titles to both axes
  theme(axis.title.x=element_text(margin=margin(b=25), vjust=-5, size=25), #move x-axis title lower and increase bottom whitespace
        axis.title.y=element_text(margin=margin(r=25), size=25),
        legend.position = "none",
        axis.text=element_text(size=20, colour="black"),
        plot.margin = unit(c(20,20,20,20), "pt")) +
  scale_x_discrete(breaks=c("lobata_healthy", "rus_healthy", "rus_bleached", "pulchra_healthy"), labels=c("P. lobata", "P. rus", "P. rus bleached", "A. pulchra")) +
  annotate(geom="text", x = 4, y=185, label=c("p<0.001"), size=7) + #add p-value
  annotate(geom="text", x = c("lobata_healthy", "rus_healthy", "rus_bleached", "pulchra_healthy"), y=c(-30,-30,-30,-30), label=c("a", "a", "b", "c"), size=7) + #add letters for significance groupings
  annotate(geom="text", x = 4.3, y=225, label=c("B"), size=9) #add plot label B


plbox
