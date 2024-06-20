#Haley Greenyer 
#modified 2022-02-01
#Network Inference Algorithm Performance Comparisons 

#packages-----------------
library(tidyverse)
library(infotheo) # for MI/clr 
library(minet) #for DCA
library(ggplot2)
library(WGCNA) #for MI
library(dplyr)
library(devtools)
library(igraph) #for outputting edge lists 
library(pROC) #for calculating auroc
library(ppcor) # for partial correlation 
library(propr)
library(Rcpp)
library(mltools)
library(apcluster)
library(MLmetrics) #for auc pr

#network inference methods
library(ggm)
library(GENIE3)
library(tigress)
library(GeneNet)

#functions
source('~/Thesis/final_code/connectivity_matrix.R') #to get primary/secondary/tertiary connection reference matrices 
source('~/Thesis/final_code/get_auroc.R')
source('~/Thesis/final_code/get_pr.R')

#inference methods
source('~/Thesis/final_code/inference_a/PCLRC_CORR.R') #PCLRC
source('~/Thesis/final_code/inference_a/LOPC_R.R') #for partial correlations 
source('~/Thesis/final_code/inference_a/CLRCR.R') 
source('~/Thesis/final_code/inference_a/test_method.R') 
source('~/Thesis/final_code/inference_a/fuzzy_clust.R') 


#Load simulated data and stoichiometric matrix-----------------
setwd("~/Thesis/final_code")

#load simulated data 

#edoardo simulation 
sim_e <- read.csv('~/Thesis/edoardo_sim/Steady_state.csv', header=TRUE, fileEncoding="UTF-8-BOM") ##SS no noise, 500 samples
#load reference matrices
prim_e <-as.matrix(read.csv('~/Thesis/edoardo_sim/Primary_adj.csv', header=TRUE, fileEncoding="UTF-8-BOM"))
sec_e <-as.matrix(read.csv('~/Thesis/edoardo_sim/Secondary_adj.csv', header=TRUE, fileEncoding="UTF-8-BOM"))
ter_e <-as.matrix(read.csv('~/Thesis/edoardo_sim/Tertiary_adj.csv', header=TRUE, fileEncoding="UTF-8-BOM"))

#load stoichiometric matrix 
stoich_m <- read.csv('~/Thesis/final_code/simulations/stoich_matlab.csv', header=FALSE, fileEncoding="UTF-8-BOM") 
stoich_m <- as.data.frame(t(stoich_m))

#get statistics for algorithms to test DO NOT RUN -----------------

#CORR: regular correlation (pearson's, spearman's and kendall's) 
#bicor: bi-weight mid-correlation 
#MI: mutual information 
#clust: a method that first clusters metabolite profiles and then applies correlation
#CLRCR: Corr-CLR with resampling (and no thresholding) (p,s,k,b)
#PCLRC: probabilistic clr with correlation (p,s,k,b)
#PCORR: partial correlation (p,s,k)
#GENIE: a bioconductor method that uses decision trees 

alg_labels <- c('CORR_p','CORR_s','CORR_k','bicor','pcor','MI', 'clr', 'GENIE', 'PCLRC_p', 'PCLRC_s','PCLRC_b','PCLRC_k')

str1 <- '~/Thesis/final_code/full_stats/'
str2 <- '.csv'

#function to generate statistics n times and
avg_tests <- function(sim_data, primary_ref, a, n, s1, s2){
  #ARGS:
  # sim_data: simulated data
  # _ref: reference networks
  # a: algorithm to be used with get_auc function
  # n: number of iterations to be run


  #get file name
  file_name <- paste(s1,a,s2)

  #run once to create file
  out <- get_auc(sim_data, primary_ref, method = a)
  write.csv(out, file = file_name, row.names = FALSE)
  #repeat and add results
  for(i in 1:(n-1)){
    #read in previous results
    out <- read.csv(file_name, header = TRUE,  fileEncoding="UTF-8-BOM")
    #add new results
    out <- out + get_auc(sim_data, primary_ref, method = a)
    write.csv(out, file = file_name, row.names = FALSE)
    print(i)
  }

  #divide by number of iterations to get avg
  out <- read.csv(file_name, header = TRUE,  fileEncoding="UTF-8-BOM")
  out <- out/n
  write.csv(out, file = file_name, row.names = FALSE)
}


avg_pr <- function(sim_data, primary_ref, a, n){
  #ARGS:
  # sim_data: simulated data
  # _ref: reference networks
  # a: algorithm to be used with get_auc function
  # n: number of iterations to be run
  
  
  #get file name
  s1 <-  '~/Thesis/final_code/PR_curve/'
  s2 <- 'P.csv'
  file_name <- paste(s1,a,s2)
  
  #run once to create file
  out <- get_pr(sim_data, primary_ref, method = a)
  write.csv(out, file = file_name, row.names = FALSE)
  #repeat and add results
  for(i in 1:(n-1)){
    #read in previous results
    out <- read.csv(file_name, header = TRUE,  fileEncoding="UTF-8-BOM")
    #add new results
    out <- out + get_pr(sim_data, primary_ref, method = a)
    write.csv(out, file = file_name, row.names = FALSE)
    print(i)
  }
  
  #divide by number of iterations to get avg
  out <- read.csv(file_name, header = TRUE,  fileEncoding="UTF-8-BOM")
  out <- out/n
  write.csv(out, file = file_name, row.names = FALSE)
}


#sim_n <- scale(sim_e)

# Get statistics for every algorithm to test

for(alg in c('PCLRC_k')){
  print(alg)
  avg_pr(sim_e, prim_e, a = 'alg', n = 100)
}

avg_pr(sim_e, prim_e, a = 'mrnet', n = 100)


n <-0
test <- 0

samples=sample(dim(sim_e)[1], 100)
df_s <- sim_e[samples,]

pred <- getPCLRC(df_s, m = 'kendall')
p_v <- get_vec(prim_e)
v <- get_vec(pred)

test <- test + PRAUC(v, p_v)
n <- n+1

test/100


#analysis of statistics (to load stats) -----------------------------

#function to load in all algorithm stat files 
getFiles <- function(folder_path){
  
  #get list of files from folder 
  files <- list.files(path = folder_path, full.names = TRUE)
  
  for(f in files){
    print(f)
  }
  
  
  df_list <- list()
  
  
  #read in files 
  for(i in 1:length(files)){
    
    df <- read.csv(files[i], header=TRUE, fileEncoding="UTF-8-BOM")
    df_list[[i]] <- df
  }
  
  
  return(df_list)
  
  
}

#function to extract specific statistic for each algorithm
merge_stats <- function(files, stat, c){
  
  #ARGS:
  # files: a list of stat data frames
  # stat: the statistic to be extracted
  # c: list of column names for output
  
  n_f <- length(files)
  
  out <- as.data.frame(matrix(nrow = nrow(files[[1]]), ncol = n_f + 1))
  
  #assign sample sizes to first column
  temp <- files[[1]]
  out[,1] <- temp[,1]
  
  i <- 2
  for(f in files){
    out[,i] <- f[,stat]
    i <- i+1
  }
  colnames(out) <- c
  return(out)
}


#get all stat file names 
stats <- getFiles('~/Thesis/final_code/full_stats/')

#should reflect order in which files are read in 
cn <- c('sample_size','bicor','clr','CORR_k', 'CORR_p','CORR_s','GENIE','MI','mrnet','PCLRC_b','PCLRC_k','PCLRC_p','PCLRC_s', 'pcor')

stat_list <- colnames(stats[[1]])

#create data frame for each statistic -
#AUROCs
p_AUROC <- merge_stats(stats, 'p_AUROC',cn)

#Thesholds 
best_thresh <- merge_stats(stats, 'best_thresh',cn)

#FDRs
p_fdr <- merge_stats(stats, 'p_fdr', cn)

#PPVs
p_ppv <- merge_stats(stats, 'p_ppv', cn)

#MCC
p_mcc <- merge_stats(stats, 'p_mcc', cn)

p_alph <- merge_stats(stats, 'alpha', cn)

p_K <- merge_stats(stats, 'KS',cn)

p_Kp <- merge_stats(stats, 'KS_p',cn)

#RUNTIME
runtime <- merge_stats(stats, 'runtime', cn)


# PLOTS ---------------------------------------------------
library(reshape2)
library(RColorBrewer)

#for legend labels
CORRp <- bquote(CORR[p])
CORRs <- bquote(CORR[s])
CORRk <- bquote(CORR[k])
CORRb <- bquote(bicor)
pCORR <- bquote(pcor)
mi <- bquote(MI)
CLR <- bquote(CLR)
MRNET <- bquote(MRNET)
genie <- bquote(GENIE3)
PCLRCp <- bquote(PCLRC[p])
PCLRCs <- bquote(PCLRC[s])
PCLRCk <- bquote(PCLRC[k])
PCLRCb <- bquote(PCLRC[b])

corr_labs <- c(pCORR, CORRp, CORRs, CORRk, CORRb)
pclrc_labs <- c(PCLRCp, PCLRCs, PCLRCk, PCLRCb)
mi_labs <- c(mi, CLR, MRNET)
all_labs1 <- c(CORRb, PCLRCb, CLR, genie)
all_labs12 <- c(CORRk, PCLRCb, CLR, genie)


#groups 
corr_g <- c('sample_size','pcor','CORR_p','CORR_s', 'CORR_k', 'bicor')

pclrc_g <- c('sample_size','PCLRC_p','PCLRC_s','PCLRC_k','PCLRC_b')


all <- c('sample_size','pcor','CORR_p','CORR_s', 'CORR_k', 'bicor','PCLRC_p','PCLRC_s','PCLRC_k','PCLRC_b',
         'MI','clr','mrnet','GENIE')

#clrs <- c("#6A6A68", "#C8A867", "#007A93", "#3BBE9D", '#796A16')
#clrs2 <- c('#DD715B', '#9BD7BE','#8B5777','#589899')

clrs_c <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", '#fb9a99')
clrs_p <- c("#e31a1c", "#fdbf6f","#ffff99", "#ff7f00")
clrs_i <- c("#6a3d9a", "#cab2d6", "#b15928")
clrs_a <- c("#33a02c", "#ff7f00", "#cab2d6", "#6b2336")
clrs_a2 <- c("#fb9a99", "#ff7f00", "#cab2d6", "#6b2336")

clrs_all <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", '#fb9a99',"#e31a1c", "#fdbf6f","#ffff99", "#ff7f00","#6a3d9a", "#cab2d6", "#b15928","#6b2336")
#note that all following plots can easily be altered for primary/secondary/tertiary

# AUROC plots  --------

#CORRELATION
pAUROC_corr_df <- melt(p_AUROC[,corr_g] ,  id.vars = 'sample_size', variable.name = 'algorithm')

y1 <- as.numeric(rev(p_AUROC[1,corr_g]) +.03)
y2 <- as.numeric(rev(p_AUROC[2,corr_g]) +.03)
y3 <- as.numeric(rev(p_AUROC[3,corr_g]) +.03)
y4 <- as.numeric(rev(p_AUROC[4,corr_g]) +.03)
y5 <- as.numeric(rev(p_AUROC[5,corr_g]) +.03)
y6 <- as.numeric(rev(p_AUROC[6,corr_g]) +.03)

ggplot(pAUROC_corr_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_c, labels = corr_labs)+
  labs(x = 'sample size', y='AUROC', fill = 'Algorithm') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16)) +
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5)) +
  annotate("text",x=c(6.36, 6.18, 6, 5.81, 5.63),y=y6[1:5],label=c('A', 'B', 'B', 'BC', 'C'), size = 4)+
  annotate("text",x=c(5.36, 5.18, 5, 4.81, 4.63),y=y5[1:5],label=c('A', 'B', 'B', 'BC', 'C'), size = 4)+
  annotate("text",x=c(4.36, 4.18, 4, 3.81, 3.63),y=y4[1:5],label=c('A', 'B', 'B', 'BC', 'C'), size = 4)+
  annotate("text",x=c(3.36, 3.18, 3, 2.81, 2.63),y=y3[1:5],label=c('A', 'AB', 'AB', 'BC', 'C'), size = 4)+
  annotate("text",x=c(2.36, 2.18, 2, 1.81, 1.63),y=y2[1:5],label=c('A', 'A', 'A', 'AB', 'B'), size = 4)+
  annotate("text",x=c(1.36, 1.18, 1, .81, .63),y=y1[1:5],label=c('A', 'A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(6.36, 5.36, 4.36),y=c(y6[1],y5[1],y4[1])+.02,label=c('*'), size = 6)

#PCLRC
pAUROC_pclrc_df <- melt(p_AUROC[,pclrc_g] ,  id.vars = 'sample_size', variable.name = 'algorithm')

y1 <- as.numeric(rev(p_AUROC[1,pclrc_g]) +.03)
y2 <- as.numeric(rev(p_AUROC[2,pclrc_g]) +.03)
y3 <- as.numeric(rev(p_AUROC[3,pclrc_g]) +.03)
y4 <- as.numeric(rev(p_AUROC[4,pclrc_g]) +.03)
y5 <- as.numeric(rev(p_AUROC[5,pclrc_g]) +.03)
y6 <- as.numeric(rev(p_AUROC[6,pclrc_g]) +.03)

ggplot(pAUROC_pclrc_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_p, labels = pclrc_labs)+
  labs( x = 'sample size', y='AUROC', fill = 'Algorithm') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.33, 6.13, 5.88, 5.65),y=y6[1:4],label=c('A', 'B', 'B', 'B'), size = 4)+
  annotate("text",x=c(5.33, 5.13, 4.88, 4.65),y=y5[1:4],label=c('A', 'B', 'B', 'B'), size = 4)+
  annotate("text",x=c(4.33, 4.13, 3.88, 3.65),y=y4[1:4],label=c('A', 'B', 'B', 'B'), size = 4)+
  annotate("text",x=c(3.33, 3.13, 2.88, 2.65),y=y3[1:4],label=c('A', 'AB', 'AB', 'B'), size = 4)+
  annotate("text",x=c(2.33, 2.13, 1.88, 1.65),y=y2[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.33, 1.13, 0.88, 0.65),y=y1[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(6.33, 5.33, 4.33),y=c(y6[1],y5[1],y4[1])+.02,label=c('*'), size = 6)


#MI
pAUROC_mi_df <- melt(p_AUROC[,c('sample_size','MI','clr', 'mrnet')] ,  id.vars = 'sample_size', variable.name = 'algorithm')

m <- c('MI','clr', 'mrnet')

y1 <- as.numeric(rev(p_AUROC[1,m]) +.03)
y2 <- as.numeric(rev(p_AUROC[2,m]) +.03)
y3 <- as.numeric(rev(p_AUROC[3,m]) +.03)
y4 <- as.numeric(rev(p_AUROC[4,m]) +.03)
y5 <- as.numeric(rev(p_AUROC[5,m]) +.03)
y6 <- as.numeric(rev(p_AUROC[6,m]) +.03)


ggplot(pAUROC_mi_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + scale_fill_manual(values = clrs_i, labels = mi_labs)+
  labs( x = 'sample size', y='AUROC', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.23, 5.99, 5.77),y=y6[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(5.23, 4.99, 4.77),y=y5[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(4.23, 3.99, 3.77),y=y4[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(3.23, 2.99, 2.77),y=y3[1:3],label=c('B', 'A', 'A'), size = 4)+
  annotate("text",x=c(2.23, 1.99, 1.77),y=y2[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.23, 0.99, 0.77),y=y1[1:3],label=c('A', 'A', 'A'), size = 4)


#ALL
pAUROC_all_df <- melt(p_AUROC[,c('sample_size','bicor', 'PCLRC_b', 'clr','GENIE')] ,  id.vars = 'sample_size', variable.name = 'algorithm')

a <- c('sample_size','bicor', 'PCLRC_b', 'clr','GENIE')

y1 <- as.numeric(rev(p_AUROC[1,a]) +.03)
y2 <- as.numeric(rev(p_AUROC[2,a]) +.03)
y3 <- as.numeric(rev(p_AUROC[3,a]) +.03)
y4 <- as.numeric(rev(p_AUROC[4,a]) +.03)
y5 <- as.numeric(rev(p_AUROC[5,a]) +.03)
y6 <- as.numeric(rev(p_AUROC[6,a]) +.03)


ggplot(pAUROC_all_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_a2, labels = all_labs1)+
  labs( x = 'sample size', y='AUROC', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=16)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.33, 6.11, 5.88, 5.65),y=y6[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(5.33, 5.11, 4.88, 4.65),y=y5[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(4.33, 4.11, 3.88, 3.65),y=y4[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(3.33, 3.11, 2.88, 2.65),y=y3[1:4],label=c('AB', 'B', 'A', 'AB'), size = 4)+
  annotate("text",x=c(2.33, 2.11, 1.88, 1.65),y=y2[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.33, 1.11, 0.88, 0.65),y=y1[1:4],label=c('A', 'A', 'A', 'A'), size = 4)



# AUC PR plots  --------

pr_curve <- read.csv('~/Thesis/final_code/PR_curve/PR_scores.csv', header=TRUE, fileEncoding="UTF-8-BOM")


#CORRELATION
pr_corr_df <- melt(pr_curve[,corr_g] ,  id.vars = 'sample_size', variable.name = 'algorithm')

y1 <- as.numeric(rev(pr_curve[1,corr_g]) +.03)
y2 <- as.numeric(rev(pr_curve[2,corr_g]) +.03)
y3 <- as.numeric(rev(pr_curve[3,corr_g]) +.03)
y4 <- as.numeric(rev(pr_curve[4,corr_g]) +.03)
y5 <- as.numeric(rev(pr_curve[5,corr_g]) +.03)
y6 <- as.numeric(rev(pr_curve[6,corr_g]) +.03)

ggplot(pr_corr_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_c, labels = corr_labs)+
  labs(x = 'sample size', y='AUPR', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16)) +
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5)) +
  annotate("text",x=c(6.36, 6.18, 6, 5.81, 5.63),y=y6[1:5],label=c('A', 'B', 'B', 'C', 'D'), size = 4)+
  annotate("text",x=c(5.36, 5.18, 5, 4.81, 4.63),y=y5[1:5],label=c('A', 'B', 'B', 'C', 'C'), size = 4)+
  annotate("text",x=c(4.36, 4.18, 4, 3.81, 3.63),y=y4[1:5],label=c('B', 'A', 'A', 'C', 'D'), size = 4)+
  annotate("text",x=c(3.36, 3.18, 3, 2.81, 2.63),y=y3[1:5],label=c('B', 'A', 'A', 'C', 'D'), size = 4)+
  annotate("text",x=c(2.36, 2.18, 2, 1.81, 1.63),y=y2[1:5],label=c('A', 'A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.36, 1.18, 1, .81, .63),y=y1[1:5],label=c('A', 'A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(6.36, 5.36),y=c(y6[1],y5[1])+.02,label=c('*'), size = 6)

#PCLRC
pr_pclrc_df <- melt(pr_curve[,pclrc_g] ,  id.vars = 'sample_size', variable.name = 'algorithm')

y1 <- as.numeric(rev(pr_curve[1,pclrc_g]) +.03)
y2 <- as.numeric(rev(pr_curve[2,pclrc_g]) +.03)
y3 <- as.numeric(rev(pr_curve[3,pclrc_g]) +.03)
y4 <- as.numeric(rev(pr_curve[4,pclrc_g]) +.03)
y5 <- as.numeric(rev(pr_curve[5,pclrc_g]) +.03)
y6 <- as.numeric(rev(pr_curve[6,pclrc_g]) +.03)

ggplot(pr_pclrc_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_p, labels = pclrc_labs)+
  labs( x = 'sample size', y='AUPR', fill = 'Algorithm') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=16)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.33, 6.13, 5.88, 5.65),y=y6[1:4],label=c('B', 'B', 'B', 'A'), size = 4)+
  annotate("text",x=c(5.33, 5.13, 4.88, 4.65),y=y5[1:4],label=c('A', 'B', 'B', 'AB'), size = 4)+
  annotate("text",x=c(4.33, 4.13, 3.88, 3.65),y=y4[1:4],label=c('A', 'C', 'C', 'B'), size = 4)+
  annotate("text",x=c(3.33, 3.13, 2.88, 2.65),y=y3[1:4],label=c('A', 'C', 'C', 'B'), size = 4)+
  annotate("text",x=c(2.33, 2.13, 1.88, 1.65),y=y2[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.33, 1.13, 0.88, 0.65),y=y1[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(5.65,4.33, 3.33),y=c(y6[4],y4[1],y3[1])+.02,label=c('*'), size = 6)


#MI
pr_mi_df <- melt(pr_curve[,c('sample_size','MI','clr', 'mrnet')] ,  id.vars = 'sample_size', variable.name = 'algorithm')

m <- c('MI','clr', 'mrnet')

y1 <- as.numeric(rev(pr_curve[1,m]) +.03)
y2 <- as.numeric(rev(pr_curve[2,m]) +.03)
y3 <- as.numeric(rev(pr_curve[3,m]) +.03)
y4 <- as.numeric(rev(pr_curve[4,m]) +.03)
y5 <- as.numeric(rev(pr_curve[5,m]) +.03)
y6 <- as.numeric(rev(pr_curve[6,m]) +.03)


ggplot(pr_mi_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + scale_fill_manual(values = clrs_i, labels = mi_labs)+
  labs( x = 'sample size', y='AUPR', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.23, 5.99, 5.77),y=y6[1:3],label=c('B', 'A', 'B'), size = 4)+
  annotate("text",x=c(5.23, 4.99, 4.77),y=y5[1:3],label=c('B', 'A', 'B'), size = 4)+
  annotate("text",x=c(4.23, 3.99, 3.77),y=y4[1:3],label=c('B', 'A', 'B'), size = 4)+
  annotate("text",x=c(3.23, 2.99, 2.77),y=y3[1:3],label=c('B', 'A', 'B'), size = 4)+
  annotate("text",x=c(2.23, 1.99, 1.77),y=y2[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.23, 0.99, 0.77),y=y1[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(5.99, 4.99, 3.99, 2.99),y=c(y6[2],y5[2],y4[2], y3[2])+.02,label=c('*'), size = 6)


 



# MCC plots -------------

#CORRELATION
mcc_corr_df <- melt(p_mcc[,corr_g] ,  id.vars = 'sample_size', variable.name = 'algorithm')

y1 <- as.numeric(rev(p_mcc[1,corr_g]) +.03)
y2 <- as.numeric(rev(p_mcc[2,corr_g]) +.03)
y3 <- as.numeric(rev(p_mcc[3,corr_g]) +.03)
y4 <- as.numeric(rev(p_mcc[4,corr_g]) +.03)
y5 <- as.numeric(rev(p_mcc[5,corr_g]) +.03)
y6 <- as.numeric(rev(p_mcc[6,corr_g]) +.03)


ggplot(mcc_corr_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_c, labels = corr_labs)+
  labs( x = 'sample size', y='MCC', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1))+ scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.36, 6.18, 6, 5.81, 5.63),y=y6[1:5],label=c('B', 'A', 'A', 'B', 'C'), size = 4)+
  annotate("text",x=c(5.36, 5.18, 5, 4.81, 4.63),y=y5[1:5],label=c('B', 'A', 'A', 'B', 'C'), size = 4)+
  annotate("text",x=c(4.36, 4.18, 4, 3.81, 3.63),y=y4[1:5],label=c('B', 'A', 'AB', 'C', 'D'), size = 4)+
  annotate("text",x=c(3.36, 3.18, 3, 2.81, 2.63),y=y3[1:5],label=c('AB', 'A', 'A', 'B', 'C'), size = 4)+
  annotate("text",x=c(2.36, 2.18, 2, 1.81, 1.63),y=y2[1:5],label=c('BC', 'AB', 'A', 'BC', 'C'), size = 4)+
  annotate("text",x=c(1.36, 1.18, 1, .81, .63),y=y1[1:5],  label=c('A', 'A', 'A', 'A', 'A'), size = 4)



#PCLRC
mcc_pclrc_df <- melt(p_mcc[,pclrc_g] ,  id.vars = 'sample_size', variable.name = 'algorithm')


y1 <- as.numeric(rev(p_mcc[1,pclrc_g]) +.03)
y2 <- as.numeric(rev(p_mcc[2,pclrc_g]) +.03)
y3 <- as.numeric(rev(p_mcc[3,pclrc_g]) +.03)
y4 <- as.numeric(rev(p_mcc[4,pclrc_g]) +.03)
y5 <- as.numeric(rev(p_mcc[5,pclrc_g]) +.03)
y6 <- as.numeric(rev(p_mcc[6,pclrc_g]) +.03)


ggplot(mcc_pclrc_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_p, labels = pclrc_labs)+
  labs( x = 'sample size', y='MCC', fill = 'Algorithms') +
  theme(plot.title = element_text( color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text( color="#666666", face="bold", size=16))+
  ylim(c(0,1))+ scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.33, 6.13, 5.88, 5.65),y=y6[1:4],label=c('A', 'C', 'C', 'B'), size = 4)+
  annotate("text",x=c(5.33, 5.13, 4.88, 4.65),y=y5[1:4],label=c('A', 'C', 'C', 'B'), size = 4)+
  annotate("text",x=c(4.33, 4.13, 3.88, 3.65),y=y4[1:4],label=c('A', 'B', 'B', 'B'), size = 4)+
  annotate("text",x=c(3.33, 3.13, 2.88, 2.65),y=y3[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(2.33, 2.13, 1.88, 1.65),y=y2[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.33, 1.13, 0.88, 0.65),y=y1[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(6.33, 5.33, 4.33),y=c(y6[1],y5[1],y4[1])+.02,label=c('*'), size = 6)



#MI
mcc_mi_df <- melt(p_mcc[,c('sample_size','MI','clr','mrnet')] ,  id.vars = 'sample_size', variable.name = 'algorithm')

m <- c('MI','clr', 'mrnet')

y1 <- as.numeric(rev(p_mcc[1,m]) +.02)
y2 <- as.numeric(rev(p_mcc[2,m]) +.02)
y3 <- as.numeric(rev(p_mcc[3,m]) +.02)
y4 <- as.numeric(rev(p_mcc[4,m]) +.02)
y5 <- as.numeric(rev(p_mcc[5,m]) +.02)
y6 <- as.numeric(rev(p_mcc[6,m]) +.02)


ggplot(mcc_mi_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + scale_fill_manual(values = clrs_i, labels = mi_labs)+
  labs( x = 'sample size', y='MCC', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.23, 5.99, 5.77),y=y6[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(5.23, 4.99, 4.77),y=y5[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(4.23, 3.99, 3.77),y=y4[1:3],label=c('B', 'A', 'A'), size = 4)+
  annotate("text",x=c(3.23, 2.99, 2.77),y=y3[1:3],label=c('B', 'A', 'A'), size = 4)+
  annotate("text",x=c(2.23, 1.99, 1.77),y=y2[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.23, 0.99, 0.77),y=y1[1:3],label=c('A', 'A', 'A'), size = 4)



#ALL
mcc_all_df <- melt(p_mcc[,c('sample_size','CORR_k', 'PCLRC_b', 'clr','GENIE')] ,  id.vars = 'sample_size', variable.name = 'algorithm')

a <- c('sample_size','CORR_k', 'PCLRC_b', 'clr','GENIE')

y1 <- as.numeric(rev(p_mcc[1,a]) +.02)
y2 <- as.numeric(rev(p_mcc[2,a]) +.02)
y3 <- as.numeric(rev(p_mcc[3,a]) +.02)
y4 <- as.numeric(rev(p_mcc[4,a]) +.02)
y5 <- as.numeric(rev(p_mcc[5,a]) +.02)
y6 <- as.numeric(rev(p_mcc[6,a]) +.02)


ggplot(mcc_all_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_a, labels = all_labs12)+
  labs( x = 'sample size', y='MCC', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.33, 6.13, 5.88, 5.65),y=y6[1:4],label=c('B', 'B', 'A', 'A'), size = 4)+
  annotate("text",x=c(5.33, 5.13, 4.88, 4.65),y=y5[1:4],label=c('B', 'B', 'A', 'A'), size = 4)+
  annotate("text",x=c(4.33, 4.13, 3.88, 3.65),y=y4[1:4],label=c('B', 'B', 'A', 'A'), size = 4)+
  annotate("text",x=c(3.33, 3.13, 2.88, 2.65),y=y3[1:4],label=c('B', 'B', 'A', 'A'), size = 4)+
  annotate("text",x=c(2.33, 2.13, 1.88, 1.65),y=y2[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.33, 1.13, 0.88, 0.65),y=y1[1:4],label=c('A', 'A', 'A', 'A'), size = 4)




#FDR plots ------------------------------

#CORRELATION
fdr_corr_df <- melt(p_fdr[,corr_g] ,  id.vars = 'sample_size', variable.name = 'algorithm')

y1 <- as.numeric(rev(p_fdr[1,corr_g]) +.03)
y2 <- as.numeric(rev(p_fdr[2,corr_g]) +.03)
y3 <- as.numeric(rev(p_fdr[3,corr_g]) +.03)
y4 <- as.numeric(rev(p_fdr[4,corr_g]) +.03)
y5 <- as.numeric(rev(p_fdr[5,corr_g]) +.03)
y6 <- as.numeric(rev(p_fdr[6,corr_g]) +.03)

ggplot(fdr_corr_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_c, labels = corr_labs)+
  labs( x = 'sample size', y='FDR', fill = 'Algorithms') +
  theme(plot.title = element_text( color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.36, 6.18, 6, 5.81, 5.63),y=y6[1:5],label=c('D', 'A', 'B', 'C', 'D'), size = 4)+
  annotate("text",x=c(5.36, 5.18, 5, 4.81, 4.63),y=y5[1:5],label=c('BC', 'A', 'A', 'B', 'C'), size = 4)+
  annotate("text",x=c(4.36, 4.18, 4, 3.81, 3.63),y=y4[1:5],label=c('B', 'A', 'A', 'B', 'C'), size = 4)+
  annotate("text",x=c(3.36, 3.18, 3, 2.81, 2.63),y=y3[1:5],label=c('B', 'A', 'A', 'B', 'B'), size = 4)+
  annotate("text",x=c(2.36, 2.18, 2, 1.81, 1.63),y=y2[1:5],label=c('A', 'A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.36, 1.18, 1, .81, .63),y=y1[1:5],  label=c('A', 'A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(6.18),y=c(y6[2])+.02,label=c('*'), size = 6)

#PCLRC
fdr_pclrc_df <- melt(p_fdr[,pclrc_g] ,  id.vars = 'sample_size', variable.name = 'algorithm')

y1 <- as.numeric(rev(p_fdr[1,pclrc_g]) +.03)
y2 <- as.numeric(rev(p_fdr[2,pclrc_g]) +.03)
y3 <- as.numeric(rev(p_fdr[3,pclrc_g]) +.03)
y4 <- as.numeric(rev(p_fdr[4,pclrc_g]) +.03)
y5 <- as.numeric(rev(p_fdr[5,pclrc_g]) +.03)
y6 <- as.numeric(rev(p_fdr[6,pclrc_g]) +.03)

ggplot(fdr_pclrc_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_p, labels = pclrc_labs)+
  labs( x = 'sample size', y='FDR', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1))+ scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.33, 6.13, 5.88, 5.65),y=y6[1:4],label=c('A', 'B', 'B', 'A'), size = 4)+
  annotate("text",x=c(5.33, 5.13, 4.88, 4.65),y=y5[1:4],label=c('A', 'B', 'B', 'A'), size = 4)+
  annotate("text",x=c(4.33, 4.13, 3.88, 3.65),y=y4[1:4],label=c('A', 'B', 'B', 'A'), size = 4)+
  annotate("text",x=c(3.33, 3.13, 2.88, 2.65),y=y3[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(2.33, 2.13, 1.88, 1.65),y=y2[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.33, 1.13, 0.88, 0.65),y=y1[1:4],label=c('A', 'A', 'A', 'A'), size = 4)


#MI
fdr_mi_df <- melt(p_fdr[,c('sample_size','MI','clr','mrnet')] ,  id.vars = 'sample_size', variable.name = 'algorithm')


m <- c('MI','clr','mrnet')

y1 <- as.numeric(rev(p_fdr[1,m]) +.03)
y2 <- as.numeric(rev(p_fdr[2,m]) +.03)
y3 <- as.numeric(rev(p_fdr[3,m]) +.03)
y4 <- as.numeric(rev(p_fdr[4,m]) +.03)
y5 <- as.numeric(rev(p_fdr[5,m]) +.03)
y6 <- as.numeric(rev(p_fdr[6,m]) +.03)

ggplot(fdr_mi_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + scale_fill_manual(values = clrs_i, labels = mi_labs)+
  labs( x = 'sample size', y='FDR', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.23, 5.99, 5.77),y=y6[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(5.23, 4.99, 4.77),y=y5[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(4.23, 3.99, 3.77),y=y4[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(3.23, 2.99, 2.77),y=y3[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(2.23, 1.99, 1.77),y=y2[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.23, 0.99, 0.77),y=y1[1:3],label=c('A', 'A', 'A'), size = 4)


#ALL
fdr_all_df <- melt(p_fdr[,c('sample_size','CORR_k', 'PCLRC_b', 'clr','GENIE')] ,  id.vars = 'sample_size', variable.name = 'algorithm')

a <- c('sample_size','CORR_k', 'PCLRC_b', 'clr','GENIE')

y1 <- as.numeric(rev(p_fdr[1,a]) +.03)
y2 <- as.numeric(rev(p_fdr[2,a]) +.03)
y3 <- as.numeric(rev(p_fdr[3,a]) +.03)
y4 <- as.numeric(rev(p_fdr[4,a]) +.03)
y5 <- as.numeric(rev(p_fdr[5,a]) +.03)
y6 <- as.numeric(rev(p_fdr[6,a]) +.03)

ggplot(fdr_all_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_a, labels = all_labs12)+
  labs( x = 'sample size', y='FDR', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5)) +
  annotate("text",x=c(6.33, 6.13, 5.88, 5.65),y=y6[1:4],label=c('C', 'C', 'B', 'A'), size = 4)+
  annotate("text",x=c(5.33, 5.13, 4.88, 4.65),y=y5[1:4],label=c('C', 'C', 'B', 'A'), size = 4)+
  annotate("text",x=c(4.33, 4.13, 3.88, 3.65),y=y4[1:4],label=c('C', 'C', 'B', 'A'), size = 4)+
  annotate("text",x=c(3.33, 3.13, 2.88, 2.65),y=y3[1:4],label=c('C', 'C', 'B', 'A'), size = 4)+
  annotate("text",x=c(2.33, 2.13, 1.88, 1.65),y=y2[1:4],label=c('C', 'C', 'B', 'A'), size = 4)+
  annotate("text",x=c(1.33, 1.13, 0.88, 0.65),y=y1[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(5.65, 4.65, 3.65, 2.65, 1.65),y=c(y6[4], y5[4], y4[4],y3[4],y2[4])+.02,label=c('*'), size = 6)




#PPV plots ------------------------------
#CORRELATION
ppv_corr_df <- melt(p_ppv[,corr_g] ,  id.vars = 'sample_size', variable.name = 'algorithm')

y1 <- as.numeric(rev(p_ppv[1,corr_g]) +.03)
y2 <- as.numeric(rev(p_ppv[2,corr_g]) +.03)
y3 <- as.numeric(rev(p_ppv[3,corr_g]) +.03)
y4 <- as.numeric(rev(p_ppv[4,corr_g]) +.03)
y5 <- as.numeric(rev(p_ppv[5,corr_g]) +.03)
y6 <- as.numeric(rev(p_ppv[6,corr_g]) +.03)

ggplot(ppv_corr_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_c, labels = corr_labs)+
  labs(x = 'sample size', y='PPV', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.36, 6.18, 6, 5.81, 5.63),y=y6[1:5],label=c('D', 'A', 'B', 'C', 'D'), size = 4)+
  annotate("text",x=c(5.36, 5.18, 5, 4.81, 4.63),y=y5[1:5],label=c('BC', 'A', 'A', 'B', 'C'), size = 4)+
  annotate("text",x=c(4.36, 4.18, 4, 3.81, 3.63),y=y4[1:5],label=c('B', 'A', 'A', 'B', 'C'), size = 4)+
  annotate("text",x=c(3.36, 3.18, 3, 2.81, 2.63),y=y3[1:5],label=c('B', 'A', 'A', 'B', 'B'), size = 4)+
  annotate("text",x=c(2.36, 2.18, 2, 1.81, 1.63),y=y2[1:5],label=c('A', 'A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.36, 1.18, 1, .81, .63),y=y1[1:5],  label=c('A', 'A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(6.18),y=c(y6[2])+.02,label=c('*'), size = 6)



#PCLRC
ppv_pclrc_df <- melt(p_ppv[,pclrc_g],  id.vars = 'sample_size', variable.name = 'algorithm')

y1 <- as.numeric(rev(p_ppv[1,pclrc_g]) +.03)
y2 <- as.numeric(rev(p_ppv[2,pclrc_g]) +.03)
y3 <- as.numeric(rev(p_ppv[3,pclrc_g]) +.03)
y4 <- as.numeric(rev(p_ppv[4,pclrc_g]) +.03)
y5 <- as.numeric(rev(p_ppv[5,pclrc_g]) +.03)
y6 <- as.numeric(rev(p_ppv[6,pclrc_g]) +.03)

ggplot(ppv_pclrc_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_p, labels = pclrc_labs)+
  labs(x = 'sample size', y='PPV', fill = 'Algorithms' ) +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.33, 6.13, 5.88, 5.65),y=y6[1:4],label=c('A', 'B', 'B', 'A'), size = 4)+
  annotate("text",x=c(5.33, 5.13, 4.88, 4.65),y=y5[1:4],label=c('A', 'B', 'B', 'A'), size = 4)+
  annotate("text",x=c(4.33, 4.13, 3.88, 3.65),y=y4[1:4],label=c('A', 'B', 'B', 'A'), size = 4)+
  annotate("text",x=c(3.33, 3.13, 2.88, 2.65),y=y3[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(2.33, 2.13, 1.88, 1.65),y=y2[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.33, 1.13, 0.88, 0.65),y=y1[1:4],label=c('A', 'A', 'A', 'A'), size = 4)


ppv_mi_df <- melt(p_ppv[,c('sample_size','MI','clr','mrnet')] ,  id.vars = 'sample_size', variable.name = 'algorithm')


m <- c('MI','clr','mrnet')

y1 <- as.numeric(rev(p_ppv[1,m]) +.02)
y2 <- as.numeric(rev(p_ppv[2,m]) +.02)
y3 <- as.numeric(rev(p_ppv[3,m]) +.02)
y4 <- as.numeric(rev(p_ppv[4,m]) +.02)
y5 <- as.numeric(rev(p_ppv[5,m]) +.02)
y6 <- as.numeric(rev(p_ppv[6,m]) +.02)

ggplot(ppv_mi_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.7) + scale_fill_manual(values = clrs_i, labels = mi_labs)+
  labs( x = 'sample size', y='PPV', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5))+
  annotate("text",x=c(6.23, 5.99, 5.77),y=y6[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(5.23, 4.99, 4.77),y=y5[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(4.23, 3.99, 3.77),y=y4[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(3.23, 2.99, 2.77),y=y3[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(2.23, 1.99, 1.77),y=y2[1:3],label=c('A', 'A', 'A'), size = 4)+
  annotate("text",x=c(1.23, 0.99, 0.77),y=y1[1:3],label=c('A', 'A', 'A'), size = 4)

ppv_all_df <- melt(p_ppv[,c('sample_size','CORR_k', 'PCLRC_b', 'clr','GENIE')] ,  id.vars = 'sample_size', variable.name = 'algorithm')

a <- c('sample_size','CORR_k', 'PCLRC_b', 'clr','GENIE')

y1 <- as.numeric(rev(p_ppv[1,a]) +.03)
y2 <- as.numeric(rev(p_ppv[2,a]) +.03)
y3 <- as.numeric(rev(p_ppv[3,a]) +.03)
y4 <- as.numeric(rev(p_ppv[4,a]) +.03)
y5 <- as.numeric(rev(p_ppv[5,a]) +.03)
y6 <- as.numeric(rev(p_ppv[6,a]) +.03)

ggplot(ppv_all_df, aes(fill=algorithm, y=value, x=as.factor(sample_size))) + 
  geom_bar(position="dodge", stat="identity", width = 0.9) + scale_fill_manual(values = clrs_a, labels = all_labs12)+
  labs( x = 'sample size', y='PPV', fill = 'Algorithms') +
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14)) +
  theme(legend.text = element_text(color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(color="#666666", face="bold", size=16))+
  ylim(c(0,1)) + scale_x_discrete(expand = c(0,.5)) +
  annotate("text",x=c(6.33, 6.13, 5.88, 5.65),y=y6[1:4],label=c('C', 'C', 'B', 'A'), size = 4)+
  annotate("text",x=c(5.33, 5.13, 4.88, 4.65),y=y5[1:4],label=c('C', 'C', 'B', 'A'), size = 4)+
  annotate("text",x=c(4.33, 4.13, 3.88, 3.65),y=y4[1:4],label=c('C', 'C', 'B', 'A'), size = 4)+
  annotate("text",x=c(3.33, 3.13, 2.88, 2.65),y=y3[1:4],label=c('C', 'C', 'B', 'A'), size = 4)+
  annotate("text",x=c(2.33, 2.13, 1.88, 1.65),y=y2[1:4],label=c('C', 'C', 'B', 'A'), size = 4)+
  annotate("text",x=c(1.33, 1.13, 0.88, 0.65),y=y1[1:4],label=c('A', 'A', 'A', 'A'), size = 4)+
  annotate("text",x=c(5.65, 4.65, 3.65, 2.65, 1.65),y=c(y6[4], y5[4], y4[4],y3[4],y2[4])+.02,label=c('*'), size = 6)




#threshold plot -------------------------

thresh_df <- melt(best_thresh[,all] ,  id.vars = 'sample_size', variable.name = 'algorithm')

ggplot(thresh_df, aes(y=value, x=as.factor(sample_size), group = algorithm), fill = clrs_all) + 
  geom_point(size = 3.5, aes(color = algorithm)) +
  geom_line(size = 1.5, aes(color = algorithm))+
  scale_color_manual(values = clrs_all)+
  labs(x = 'sample size', y='Threshold') +
  theme(plot.title = element_text(family = "mono", color="#666666", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(family = "mono", color="#666666", face="bold", size=16)) +
  theme(legend.text = element_text(family = "mono", color="#666666", face="bold", size=16))+
  theme(legend.title = element_text(family = "mono", color="#666666", face="bold", size=16))


#permutation tests -----------------------
source('~/Thesis/final_code/permutations.R')

ss <- c(5,10,20,50,100,200)


#test <- permutation_test(sim_e, ss = 100, m1 = 'MI', m2 = 'mrnet', stat = 'mcc', n = 50, ref_net = prim_e, d = .5 )

#<- function(df, ss, m1, m2, n = 100, stat, ref_net, d)

#ran fine!
stat_perm<- matrix(nrow = length(ss), ncol = 3)
stat_perm[,] <- 0
stat_perm[,1] <- ss

to_test <- c('sample_size', 'clr','GENIE')

n = 5
for(r in 1:n){
  #print(r)
  j <- 2
  for(m in to_test[2:length(to_test)]){
    i <- 1
    print(m)
    for(s in ss){
      #print(s)
      diff <- abs(pr_curve[i,m] - pr_curve[i,'PCLRC_b'])
      #print(diff)
      temp <- permutation_test(sim_e, ss = s, m1 = m, m2 = 'PCLRC_b', stat = 'pr', n = 10, ref_net = prim_e, d=diff )
      #print(temp)
      stat_perm[i,j] <- stat_perm[i,j] + temp 
      i <- i+1
    }
    j <- j+1
  }
  
}
stat_perm <- stat_perm/n

colnames(stat_perm) <- to_test

write.csv(stat_perm,'~/Thesis/final_code/full_perm/all_pr.csv',row.names = FALSE)

diff <- abs(pr_curve[i,'clr'] - pr_curve[i,'mrnet'])
temp <- permutation_test(sim_e, ss = 200, m1 = 'mrnet', m2 = 'clr', stat = 'pr', n = 1, ref_net = prim_e, d=diff )
print(temp)


