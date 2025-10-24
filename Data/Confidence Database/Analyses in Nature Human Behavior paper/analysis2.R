#--------------------------------------------------------------------------
# Analysis 2 from "The Confidence Database" paper by Rahnev, Desender, Lee, et al.
# 
# This analysis explores serial dependence in confidence RTs, up to lag 7.
# This is done on all datasets including this variable.
#
# To run this analysis, all the files of the Confidence Database should be placed in a folder called 'Confidence Database' located in your current WD
#
# Written by Kobe Desender. Last update: Feb, 2020.
#--------------------------------------------------------------------------
rm(list=ls());library(here);library(pwr)

#Set data path 
rootPath = "C:/Users/kobe/Documents/projecten/The Confidence Database"
setwd(rootPath)
dataPath = paste0(rootPath, "/Confidence Database/")

#This function easily merges two dataframes
source('fastmerge.R')
#This function computes VIFs for mixed models
source('vif_mer.R')

#Get the names of the files that will be read in
data_files <- list.files(path=dataPath,pattern = "^data_*", recursive = FALSE)
readme_files <- list.files(path=dataPath,pattern = "^readme_*", recursive = FALSE)

#################################################
#Analysis 2: Serial dependence in confidence RT #
#################################################
#Preallocate one big dataframe to save all data in (caution, this takes quite a while)
DataAll <- data.frame(matrix(NA,1,0));which_data_included <- rep(NA,0)
counter <- 1 #participant counter 
counter2 <- 1 #study counter

#Loop over all files and add each dataset to the big DF
for(i in 1:length(data_files)){
  Data <- read.csv(paste0(dataPath,data_files[i]),fileEncoding="UTF-8-BOM")
  
  #Only include a study if it has a column RT_conf
  if(any(names(Data)=="RT_conf")){
    
      #loop over subjects in this dataset
      subs <- unique(Data$Subj_idx)  
      for(j in 1:length(subs)){
        temp <- subset(Data,Subj_idx==subs[j])
        temp <- temp[,c('Subj_idx','RT_conf')]
        
        #create lagged trials
        temp$RT_conf_1[2:length(temp$RT_conf)] <- temp$RT_conf[1:(length(temp$RT_conf)-1)]
        temp$RT_conf_2[3:length(temp$RT_conf)] <- temp$RT_conf[1:(length(temp$RT_conf)-2)]
        temp$RT_conf_3[4:length(temp$RT_conf)] <- temp$RT_conf[1:(length(temp$RT_conf)-3)]
        temp$RT_conf_4[5:length(temp$RT_conf)] <- temp$RT_conf[1:(length(temp$RT_conf)-4)]
        temp$RT_conf_5[6:length(temp$RT_conf)] <- temp$RT_conf[1:(length(temp$RT_conf)-5)]
        temp$RT_conf_6[7:length(temp$RT_conf)] <- temp$RT_conf[1:(length(temp$RT_conf)-6)]
        temp$RT_conf_7[8:length(temp$RT_conf)] <- temp$RT_conf[1:(length(temp$RT_conf)-7)]
        
        #Mixed models can deal with NaNs, so no need to exclude any NaNs
        
        #Only include data with (at least) more than 7 trials, otherwise it's not possible to estimate lag 7
        if(dim(temp)[1]>7){
          temp$Subj_idx <- counter #create a new participant counter
          DataAll <- fastmerge(DataAll,temp) #merge this paticipant with all data
          counter <- counter+1
        }
      }
      which_data_included[counter2] <- i;counter2 <- counter2+1
    }
  print(paste('merging dataset ',i,'from ',length(data_files)))
}
DataAll <- DataAll[2:length(DataAll$Subj_idx),]

#Which data are included
data_files[which_data_included]
length(which_data_included)
length(data_files)
length(unique(DataAll$Subj_idx))


#Exclude outliers & conf RTs < 0 or > 5
RTdeadline <- 5 #5 before
subs <- unique(DataAll$Subj_idx)
DataAll$outlier <- 0
for(i in 1:length(subs)){
  temp <- subset(DataAll,Subj_idx==i)
  for(j in 1:8){
    DataAll$outlier[DataAll$Subj_idx==i][temp[,j+1]>RTdeadline] <- 1 #Slower than 5s
    DataAll$outlier[DataAll$Subj_idx==i][temp[,j+1]<0.001] <- 2 #No RT
    DataAll$outlier[DataAll$Subj_idx==i][scale(temp[,j+1])>3] <- 3 #slower than 3sds
  }
  print(paste('running ', i, 'from ', length(subs)))
}
sum(table(DataAll$outlier)[2:3])/dim(DataAll)[1]
table(DataAll$outlier)[4]/dim(DataAll)[1]

#Select data
DataSelect <- subset(DataAll,outlier==0)

#Run mixed model
#Note that here, we treated all participants as coming from the same experiment. Ideally, one would model that participants are nested within experiment, using the following random-effect structure: (1|Experiment/Subj_idx), and also build a random-effects structure for the fixed effects, e.g., (RT_conf_1 + ...|Experiment/Subj_idx).
library(lmerTest);library(multcomp)
fit <- lmer(RT_conf ~ RT_conf_1 + RT_conf_2 + RT_conf_3 + RT_conf_4 + RT_conf_5 + RT_conf_6 + RT_conf_7 + (1|Subj_idx),data=DataSelect)
vif.mer(fit) #check VIFs for the predictors
summary(fit) #model output

# Extract the fixed effect estimates & confints, and plot these
tmp <- as.data.frame(confint(glht(fit))$confint)[2:8,]
tmp$History <- c('RTconf n-1','RTconf n-2','RTconf n-3','RTconf n-4','RTconf n-5','RTconf n-6','RTconf n-7')
ggplot(tmp, aes(x = History, y = Estimate, ymin = lwr, ymax = upr)) +
  geom_errorbar() + geom_point()  + 
  ylim(0, .15) + 
  theme(axis.text=element_text(size=12),axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16))








