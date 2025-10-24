#--------------------------------------------------------------------------
# Analysis 1 from "The Confidence Database" paper by Rahnev, Desender, Lee, et al.
# 
# This analysis explores the relation between confidence on the one hand, and choice RTs and confidence RTs on the other hand.
# This is done on all datasets including these variables.
#
# To run this analysis, all the files of the Confidence Database should be placed in a folder called 'Confidence Database' located in your current WD
#
# Written by Kobe Desender. Last update: Sep 16, 2019.
#--------------------------------------------------------------------------
rm(list=ls());library(here);library(pwr)

#Set data path 
rootPath = "..." #add your current wd here
setwd(rootPath)
dataPath = paste0(rootPath, "/Confidence Database/")

#This function easily merges two dataframes
source('fastmerge.R')
#This function computes VIFs for mixed models
source('vif_mer.R')

#Get the names of the files that will be read in
data_files <- list.files(path=dataPath,pattern = "^data_*", recursive = FALSE)
readme_files <- list.files(path=dataPath,pattern = "^readme_*", recursive = FALSE)


#########################################################################################
#Analysis 1: The relation between confidence and choice/confidence response times (RTs) #
#########################################################################################
#pre-allocate empty variables 
rt_conf <- rep(NA,0);rtconf_conf <- rep(NA,0);p_outlier1 <- rep(NA,0);p_outlier2 <- rep(NA,0);which_data_included <- rep(NA,0)
counter <- 1 #participant counter
counter2 <- 1 #study counter

#Loop over all files
for(i in 1:length(data_files)){
  Data <- read.csv(paste0(dataPath,data_files[i]),fileEncoding="UTF-8-BOM")
  
  #Only include data with the variables RT_conf, RT_dec and Confidence
  if(any(names(Data)=="RT_conf")&any(names(Data)=="RT_dec"&any(names(Data)=="Confidence"))){
    
    #loop over subjects in this dataset
    subs <- unique(Data$Subj_idx)  
    for(j in 1:length(subs)){
      temp <- subset(Data,Subj_idx==subs[j])
      
      #compute standardized RTs
      temp$RT_dec_z <- scale(temp$RT_dec)
      temp$RT_conf_z <- scale(temp$RT_conf)
      
      #exclude trials with NaNs for each of these variables
      n1 <- dim(temp)[1]
      temp <- temp[complete.cases(temp$Confidence),]
      temp <- temp[complete.cases(temp$RT_dec),]
      temp <- temp[complete.cases(temp$RT_conf),]
      #exclude trials with RT==0 and slower than 5s
      temp <- subset(temp, c(RT_dec > 0.01 & RT_dec < 5))
      temp <- subset(temp, c(RT_conf > 0.01 & RT_conf < 5))
      p_outlier1[counter] <- dim(temp)[1]/n1 #how many data is retained after 100ms < RT < 5s
      
      #exclude trials exceeding + or - 3sds from the grand average
      n2 <- dim(temp)[1]
      temp <- subset(temp, c(RT_dec_z < 3 & RT_dec_z > -3))
      temp <- subset(temp, c(RT_conf_z < 3 & RT_conf_z > -3))
      p_outlier2[counter] <- dim(temp)[1]/n2 #how many data is retained after excluding 3sds
      
      #Only compute correlations for participants with variation in confidence
      if(length(unique(temp$Confidence))>1){
        #An (absolute) minimum of 3 trials is needed to compute correlations
        if(dim(temp)[1]>2){
          
          #Compute correlations Confidence/Choice RT, and confidence
          temp$Confidence <- as.numeric(temp$Confidence)
          rt_conf[counter] <- cor(temp$RT_dec,temp$Confidence)
          rtconf_conf[counter] <- cor(temp$RT_conf,temp$Confidence)
          counter <- counter+1
        }
      }
      j<-j+1
    }
    which_data_included[counter2] <- i; counter2 <- counter2+1
  }
  print(paste('Computing dataset ',i,'from ',length(data_files)))
}

#Which data are included
data_files[which_data_included]
length(which_data_included) #studies included
length(data_files) #studies in the database
length(rtconf_conf) #final N

#Average number of outliers
1-mean(p_outlier1) #faster/slower than 0.01 or 5
1-mean(p_outlier2,na.rm=T) #slower/faster than 3SD

#Do these correlatie with each other
par(mar = c(5.6, 5, 4.1, 2.1))
plot(rt_conf,rtconf_conf,col="black",bg="white",ylim=c(-1,1),xlim=c(-1,1),
     ylab=expression(~italic("r")["confidence, confidence RTs"]),xlab=expression(~italic("r")["confidence, choice RTs"]),
     cex.lab=1.5,pch=21,frame=F)
abline(h=0,lty=3,col='grey',lwd=2); abline(v=0,lty=3,col='grey',lwd=2); 
abline(lm(rtconf_conf~rt_conf),col="black",lwd=3,lty=3)
cor.test(rtconf_conf,rt_conf) #r=.203

#Summary statistics
t.test(rt_conf,conf.level=.999)
sd(rt_conf)
d <- mean(rt_conf)/sd(rt_conf) #cohen's d
pwr.t.test(9,d,.05,NULL,type='one.sample') #power analysis

t.test(rtconf_conf,conf.level=.999)
sd(rtconf_conf)
d <- mean(rtconf_conf)/sd(rtconf_conf) #cohen's d
pwr.t.test(33,d,.05,NULL,type='one.sample') #power analysis
pwr.t.test(NULL,d,.05,.95,type='one.sample') #power analysis


#Plot these correlations using a rainbow plot
library(plotly);library(reshape);
source('hackySolution.R') #
df_wide <- cbind(rt_conf,rtconf_conf)
df <- melt(df_wide)
names(df) <- c("Subjects","Condition","r")
lb <- function(x) mean(x) - sd(x);ub <- function(x) mean(x) + sd(x)
sumld<- ddply(df, ~Condition, summarise, mean = mean(r), median = median(r), lower = lb(r), upper = ub(r))

g <- ggplot(data = df, aes(y = r, x = Condition, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = r, color = Condition), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  # coord_flip() +
  theme_bw() +
  raincloud_theme
g






