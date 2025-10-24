#--------------------------------------------------------------------------
# Analysis 3 from "The Confidence Database" paper by Rahnev, Desender, Lee, et al.
# 
# This analysis explores negative metacognitive sensitivity and potential explanations for this 
# This is done on all datasets including the variables accuracy and confidence
#
# To run this analysis, all the files of the Confidence Database should be placed in a folder called 'Confidence Database' located in your current WD
#
# Written by Kobe Desender. Last update: Sep 16, 2019.
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
#Analysis 3: Negative metacognitive sensitivity #
#################################################
#pre-allocate empty variables 
acc_conf <- rep(NA,0);p_outlier <- rep(NA,0);ntrials <- rep(NA,0);mean_acc <- rep(NA,0);median_rt <- rep(NA,0);median_confrt <- rep(NA,0);pMostCommonConf <- rep(NA,0);which_data_included <- rep(NA,0);whichStudy <- rep(NA,0)
counter <- 1 #ppt counter
counter2 <- 1 #experiment counter

#Loop over studies
for(i in 1:length(data_files)){
  Data <- read.csv(paste0(dataPath,data_files[i]),fileEncoding="UTF-8-BOM")
  
  #Only include datasets containing variables confidence and accuracy
  if(any(names(Data)=="Confidence")&any(names(Data)=="Accuracy")){
    
    #to exclude studies on subjective difficulty, where it doesn't make sense to look at accuracy~confidence
    if(Data$Confidence[100]!='easy' & Data$Confidence[100]!='diff'){
      
      #loop over subjects in this dataset
      subs <- unique(Data$Subj_idx)  
      for(j in 1:length(subs)){
        temp <- subset(Data,Subj_idx==subs[j])
        
        #drop trials with missing data
        temp <- temp[complete.cases(temp$Accuracy),]
        temp <- temp[complete.cases(temp$Confidence),]
        
        #Only include data when (at least) more than 2 trials
        if(dim(temp)[1]>2){
          
          #to exclude constant confidence responses
          if(length(unique(temp$Confidence))>1){
            
            #Exclude datasets with continouus measures who don't have 0/1 accuracy
            if(length(unique(temp$Accuracy)) < 3){
              
              #some datasets have words instead of numerics
              if(temp$Accuracy[1]=="Correct"|temp$Accuracy[1]=="Error"){
                temp$acc_temp <- temp$Accuracy
                temp$Accuracy <- NA
                temp$Accuracy[temp$acc_temp=="Correct"] <- 1
                temp$Accuracy[temp$acc_temp=="Error"] <- 0
              }
              
              #scale confidence to get standardized betas and perform logistic regression
              temp$Confidence <- scale(temp$Confidence)
              acc_conf[counter] <- glm(Accuracy~Confidence,data=temp)$coefficients[2]
              
              #save some other variables
              ntrials[counter] <- length(temp$Confidence)
              mean_acc[counter] <- mean(temp$Accuracy)
              pMostCommonConf[counter] <- max(table(temp$Confidence))/dim(temp)[1]
              
              #Save median choice and confidence RTs if these exist
              if(any(names(temp)=="RT_dec")){
                median_rt[counter] <- median(temp$RT_dec)
              }else{median_rt[counter] <- NA}
              if(any(names(temp)=="RT_conf")){
                median_confrt[counter] <- median(temp$RT_conf)
              }else{median_confrt[counter] <- NA}
              whichStudy[counter] <- i
              
              counter <- counter+1
            }
          }
        }
      }
      which_data_included[counter2] <- i;counter2 <- counter2+1
    }
  }
  print(paste('merging dataset ',i,'from ',length(data_files)))
}

#Which data are included
data_files[which_data_included]
length(which_data_included)
length(data_files)
length(acc_conf)

#Summary statistics
t.test(acc_conf)
plot(acc_conf,col=whichStudy)
abline(h=0,lty=2)
length(which(acc_conf<0))/length(acc_conf)

#proportion of negative betas per study
acc_conf_all <- data.frame(cbind(acc_conf,whichStudy))  
acc_conf_all$negative_beta <- 0
acc_conf_all$negative_beta[acc_conf_all$acc_conf<0] <- 1

perStudy <- with(acc_conf_all,aggregate(negative_beta,by=list(whichStudy),mean))
plot(perStudy$x) 
abline(h=.25,lty=2) #only 1 study with >%25


#Plot this
library(plotly);library(reshape);source('hackySolution.R')
df_wide <- cbind(acc_conf)
df <- melt(df_wide)
names(df) <- c("Subjects","Condition","beta")
lb <- function(x) mean(x) - sd(x); ub <- function(x) mean(x) + sd(x)
sumld<- ddply(df, ~Condition, summarise, mean = mean(beta), median = median(beta), lower = lb(beta), upper = ub(beta))

g <- ggplot(data = df, aes(y = beta, x = Condition, fill = Condition)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = beta, color = Condition), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 3.75) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  #coord_flip() +
  theme_bw() +
  raincloud_theme
g


#################################
#Link acc_conf to other variables
#################
#number of trials
par(mar = c(5.6, 5, 4.1, 2.1))
plot(acc_conf~ntrials,col="black",bg="white",ylim=c(-.3,.4),xlim=c(0,4200),
     xlab="# trials",ylab=expression(~italic("b")["accuracy~confidence"]),
     cex.lab=1.5,pch=21,frame=F)
abline(h=0,lty=3,col='grey',lwd=2); 
abline(lm(acc_conf~ntrials),col="black",lwd=3,lty=3)
cor.test(acc_conf,ntrials)

#How many people with acc_conf<0 also have ntrials<50
length(acc_conf[acc_conf<0&ntrials<50])/length(acc_conf[acc_conf<0])


################################
# Average accuracy
par(mar = c(5.6, 5, 4.1, 2.1))
plot(acc_conf~mean_acc,col="black",bg="white",ylim=c(-.3,.4),xlim=c(0,1),
     xlab="Average accuracy (%)",ylab=expression(~italic("b")["accuracy~confidence"]),
     cex.lab=1.5,pch=21,frame=F)
abline(h=0,lty=3,col='grey',lwd=2); abline(v=.5,lty=3,col='grey',lwd=2); 
abline(lm(acc_conf~mean_acc),col="black",lwd=3,lty=3)
cor.test(acc_conf,mean_acc)

#How many people with acc_conf<0 also have accuracy < 55%
length(mean_acc[acc_conf<0&mean_acc<.55])/length(mean_acc[acc_conf<0])


################################
# median RT
plot(acc_conf[!is.na(median_rt)&median_rt<5]~median_rt[!is.na(median_rt)&median_rt<5],col="black",bg="white",ylim=c(-.3,.4),xlim=c(0,5),
     xlab="Median choice RT (s)",ylab=expression(~italic("b")["accuracy~confidence"]),
     cex.lab=1.5,pch=21,frame=F)
abline(h=0,lty=3,col='grey',lwd=2); abline(v=.5,lty=3,col='grey',lwd=2); 
abline(lm(acc_conf[!is.na(median_rt)&median_rt<5]~median_rt[!is.na(median_rt)&median_rt<5]),col="black",lwd=3,lty=3)
cor.test(acc_conf[!is.na(median_rt)&median_rt<5],median_rt[!is.na(median_rt)&median_rt<5])

#How many people with acc_conf<0 also have RT < .2
length(median_rt[!is.na(median_rt)&median_rt<5&acc_conf<0&median_rt<.2])/length(median_rt[!is.na(median_rt)&median_rt<5&acc_conf<0])

################################
# median confrt
plot(acc_conf[!is.na(median_confrt)&median_confrt<5]~median_confrt[!is.na(median_confrt)&median_confrt<5],col="black",bg="white",ylim=c(-.3,.4),xlim=c(0,5),
     xlab="Median confidence RT (s)",ylab=expression(~italic("b")["accuracy~confidence"]),
     cex.lab=1.5,pch=21,frame=F)
abline(h=0,lty=3,col='grey',lwd=2); abline(v=.5,lty=3,col='grey',lwd=2); 
abline(lm(acc_conf[!is.na(median_confrt)&median_confrt<5]~median_confrt[!is.na(median_confrt)&median_confrt<5]),col="black",lwd=3,lty=3)
cor.test(acc_conf[!is.na(median_confrt)&median_confrt<5],median_confrt[!is.na(median_confrt)&median_confrt<5])

#How many people with acc_conf<0 also have RT < .2
length(median_confrt[!is.na(median_confrt)&median_confrt<5&acc_conf<0&median_confrt<.2])/length(median_confrt[!is.na(median_confrt)&median_confrt<5&acc_conf<0])

################################
#most common confidence judgment
mean(pMostCommonConf);sd(pMostCommonConf)
par(mar = c(5.6, 5, 4.1, 2.1))
plot(acc_conf~pMostCommonConf,col="black",bg="white",ylim=c(-.3,.4),xlim=c(0,1),
     xlab="p(most common confidence)",ylab=expression(~italic("b")["accuracy~confidence"]),
     cex.lab=1.5,pch=21,frame=F)
abline(h=0,lty=3,col='grey',lwd=2); 
abline(lm(acc_conf~pMostCommonConf),col="black",lwd=3,lty=3)
cor.test(acc_conf,pMostCommonConf)

#How many people with acc_conf<0 only use one confidence rating for 95% of the time
length(pMostCommonConf[acc_conf<0&pMostCommonConf>.95])/length(pMostCommonConf[acc_conf<0])


########################
#Sum up: How many people with acc_conf<0 have acc<.55,pconf>.95 or ntrials<50
length(median_rt[acc_conf < 0]) #272 people with negative betas
length(median_rt[acc_conf<0 & c(mean_acc<.55 | pMostCommonConf>.95 | ntrials < 50)]) #81 of these have acc<.55 or pconf>.95 or ntrials<50

#How many people with acc_conf<0 have ANY of the above
#176 with negative betas
length(median_rt[c(!is.na(median_rt) & median_rt<5 & !is.na(median_confrt) & median_confrt<5 & acc_conf<0)])
#55 with one of the above
length(median_rt[c(!is.na(median_rt) & median_rt<5 & !is.na(median_confrt) & median_confrt<5 & acc_conf<0) & c(median_rt<.2 | median_confrt<.2 | mean_acc<.55 | pMostCommonConf>.95 | ntrials < 50)])


#Which people with negative betas that have acc_conf<0 but acc>.55 & pconf>.95& ntrials>50 don't have rt & confrt data
table(is.na(median_rt[acc_conf<0 & c(mean_acc>=.55 & pMostCommonConf<=.95 & ntrials >= 50)])) #49
table(is.na(median_confrt[acc_conf<0 & c(mean_acc>=.55 & pMostCommonConf<=.95 & ntrials >= 50)])) #49







