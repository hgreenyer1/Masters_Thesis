#example experimental network analysis 
#Haley Greenyer 
#packages ----------------
library(tidyverse)
library(DiffCorr)
library(infotheo) # for MI
library(minet) #for DCA
library(rlist)
library(ggplot2)
library(WGCNA) #for MI
library(dplyr)
library(data.table) #for filtering columns
library(stringr)


library(devtools)
library(igraph) #for outputting edge lists 
library(pROC) #for calculating auroc

source('~/Thesis/final_code/inference_a/PCLRC_CORR.R') #PCLRC
source('~/Thesis/r_files/CalcDiffConn_final.R') #PCLRC

# library(ggm)
# library(GENIE3)
#library(tigress)
# library(GeneNet)

#apply thresholds 
apply_thresh <- function(net, thresh){
  x <- net 
  x[abs(net) < thresh] <- 0
  x[abs(net) >= thresh] <- 1
  return(x)
}


#set working directory 
setwd("~/Thesis/r_files")
#
#Load and clean NRC data -----------------------------

#function to replace NAs with five fold less than the column min 
fix_miss <- function(df){
  out <- matrix(nrow = nrow(df), ncol = ncol(df))
  for(c in 1:ncol(df)){
    m <- min(df[,c], na.rm = TRUE)/5
    out[,c] <- replace_na(df[,c], m)
  }
  return(out)
}

nrc_raw <- read.csv('~/Thesis/NRC_data/data_fixed.csv', header=TRUE)
#remove A9947 (fasted)
nrc_raw <- nrc_raw[-62,]
#replace 0 with NA
nrc_raw[nrc_raw == 0] <- NA
#remove features that have more than 60% NAs 
nrc_full <- nrc_raw[, which(colMeans(!is.na(nrc_raw)) > 0.6)]

#replace NAs with minimum value 
nrc_full[,10:ncol(nrc_full)] <- fix_miss(nrc_full[,10:ncol(nrc_full)])

#data transformations ------------------

#checking scewness -- if between -0.5, .5 then symmetric 
#skewness(nrc_full[,40])
#visualize skew 
#hist(nrc_full[,170])
get_logNorm <- function(data){
  log_d <- log(data)
  scale_d <- scale(log_d)
  return(scale_d)
}


#log transformation to correct for skewness 
nrc_log <- log(nrc_full[10:length(nrc_full)])

#z score standardization + log
nrc_full[,10:length(nrc_full)] <- scale(nrc_log)

#subsetting --------------------------------

#fasting vs. non-fasting 
fast <- subset(nrc_full, nrc_full$Fasting.state == 'Fasting')
non_fast <- subset(nrc_full, nrc_full$Fasting.state == 'Non-fasting')

#female vs male 
female <- subset(nrc_full, nrc_full$Sex == 'Female')
male <- subset(nrc_full, nrc_full$Sex == 'Male')

#Tg vs non-Tg
tg <- subset(nrc_full, nrc_full$Genotype =='Tg')
non_tg <- subset(nrc_full, nrc_full$Genotype =='NonTg')

#smaller subsets 
f_tg <- subset(fast, fast$Genotype == 'Tg')

nf_tg <- subset(non_fast, fast$Genotype == 'Tg')

f_ntg <- subset(fast, non_fast$Genotype == 'NonTg')

nf_ntg <- subset(non_fast, non_fast$Genotype == 'NonTg')



#network inference------------------------

#SIMULATED DATA ----------
sim_e <- read.csv('~/Thesis/edoardo_sim/Steady_state.csv', header=TRUE, fileEncoding="UTF-8-BOM") 
prim_e <-as.matrix(read.csv('~/Thesis/edoardo_sim/Primary_adj.csv', header=TRUE, fileEncoding="UTF-8-BOM"))
# ter_e <-as.matrix(read.csv('~/Thesis/edoardo_sim/Tertiary_adj.csv', header=TRUE, fileEncoding="UTF-8-BOM"))

# #randomly sample from simulated data
samples=sample(dim(sim_e)[1], 20)
df_s <- sim_e[samples,]
colnames(df_s) <- colnames(sim_e)

sim_pclrc_b <- getPCLRC(df_s, m='bicor')
sim_corr_k <- abs(cor(df_s, method = 'kendall'))

diag(sim_pclrc_b) <- 0
diag(sim_corr_k) <- 0

sim_pclrc_b <- apply_thresh(sim_pclrc_b, 0.64)
sim_corr_k <- apply_thresh(sim_corr_k, 0.62)

#get igraph objects
sim_pclrc_b_ig <- graph_from_adjacency_matrix(sim_pclrc_b, mode = 'undirected', diag = FALSE)
sim_corr_k_ig <- graph_from_adjacency_matrix(sim_corr_k, mode = 'undirected', diag = FALSE)

#reference network
prim_e_ig <- graph_from_adjacency_matrix(prim_e, mode = 'undirected', diag = FALSE)


#get edgelists 
sim_pclrc_b_e <- as_edgelist(sim_pclrc_b_ig)
sim_corr_k_e <- as_edgelist(sim_corr_k_ig)

prim_e_e <- as_edgelist(prim_e_ig)



#DEGREE
sim_pclrc_b_DG <- degree(sim_pclrc_b_ig)
sim_corr_k_DG <- degree(sim_corr_k_ig)


#write csvs
write.csv(sim_pclrc_b_e, '~/Thesis/final_code/example_networks/20sim_pclrc_b_e.csv', row.names = FALSE)
write.csv(sim_corr_k_e, '~/Thesis/final_code/example_networks/20sim_corr_k_e.csv', row.names = FALSE)
write.csv(prim_e_e, '~/Thesis/final_code/example_networks/prim_e_e.csv', row.names = FALSE)

#betweenness
bt_prim <- betweenness(prim_e_ig)
bt_pclrc <- betweenness(sim_pclrc_b_ig)
bt_corr <- betweenness(sim_corr_k_ig)

write.csv(bt_prim, '~/Thesis/final_code/example_networks/between_ref.csv')
write.csv(bt_pclrc, '~/Thesis/final_code/example_networks/between_pclrc.csv')
write.csv(bt_corr, '~/Thesis/final_code/example_networks/between_corr.csv')

#closeness
c_prim <- closeness(prim_e_ig)
c_pclrc <- closeness(sim_pclrc_b_ig)
c_corr <- closeness(sim_corr_k_ig)

write.csv(c_prim, '~/Thesis/final_code/example_networks/close_ref.csv')
write.csv(c_pclrc, '~/Thesis/final_code/example_networks/close_pclrc.csv')
write.csv(c_corr, '~/Thesis/final_code/example_networks/close_corr.csv')

#eigenvector centrality 
e_prim <- eigen_centrality(prim_e_ig)
e_pclrc <- eigen_centrality(sim_pclrc_b_ig)
e_corr <- eigen_centrality(sim_corr_k_ig)

write.csv(e_prim$vector, '~/Thesis/final_code/example_networks/eigen/eigen_ref.csv')
write.csv(e_pclrc$vector, '~/Thesis/final_code/example_networks/eigen/eigen_pclrc.csv')
write.csv(e_corr$vector, '~/Thesis/final_code/example_networks/eigen/eigen_corr.csv')


#modularity(prim_e_ig)


#EDGE DENSITY 
sim_pclrc_b_ED <- edge_density(sim_pclrc_b_ig)
sim_corr_k_ED <- edge_density(sim_corr_k_ig)

#EXPERIMENTAL DATA --------------

#fasted tg and non-tg
ftg_pclrc_b <- getPCLRC(f_tg[,10:length(f_tg)], m='bicor')
ftg_corr_k <- abs(cor(f_tg[,10:length(f_tg)], method='kendall'))
fntg_pclrc_b <- getPCLRC(f_ntg[,10:length(f_tg)], m='bicor')
fntg_corr_k <- abs(cor(f_ntg[,10:length(f_tg)], method='kendall'))

diag(ftg_corr_k) <- 0
diag(ftg_pclrc_b) <- 0
diag(fntg_corr_k) <- 0
diag(fntg_pclrc_b) <- 0

#THRESHOLDING
ftg_pclrc_b <- apply_thresh(ftg_pclrc_b, 0.64)
ftg_corr_k <- apply_thresh(ftg_corr_k, 0.62)
fntg_pclrc_b <- apply_thresh(fntg_pclrc_b, 0.64)
fntg_corr_k <- apply_thresh(fntg_corr_k, 0.62)


#get igraphs
ftg_pclrc_b_ig <- graph_from_adjacency_matrix(ftg_pclrc_b, mode = 'undirected', diag = FALSE)
ftg_corr_k_ig <- graph_from_adjacency_matrix(ftg_corr_k, mode = 'undirected', diag = FALSE)
fntg_pclrc_b_ig <- graph_from_adjacency_matrix(fntg_pclrc_b, mode = 'undirected', diag = FALSE)
fntg_corr_k_ig <- graph_from_adjacency_matrix(fntg_corr_k, mode = 'undirected', diag = FALSE)

#get edgelists 
ftg_pclrc_b_e <- as_edgelist(ftg_pclrc_b_ig)
ftg_corr_k_e <- as_edgelist(ftg_corr_k_ig)
fntg_pclrc_b_e <- as_edgelist(fntg_pclrc_b_ig)
fntg_corr_k_e <- as_edgelist(fntg_corr_k_ig)


#write csvs
write.csv(ftg_pclrc_b_e, '~/Thesis/final_code/example_networks/ftg_pclrc_b_e.csv')
write.csv(ftg_corr_k_e, '~/Thesis/final_code/example_networks/ftg_corr_k_e.csv')
write.csv(fntg_pclrc_b_e, '~/Thesis/final_code/example_networks/fntg_pclrc_b_e.csv')
write.csv(fntg_corr_k_e, '~/Thesis/final_code/example_networks/fntg_corr_k_e.csv')


#betweenness

ftg_pclrc_b_B <-  betweenness(ftg_pclrc_b_ig)
ftg_corr_k_B <- betweenness(fntg_pclrc_b_ig)
fntg_pclrc_b_B <- betweenness(ftg_corr_k_ig)
fntg_corr_k_B <- betweenness(fntg_corr_k_ig)

write.csv(ftg_pclrc_b_B, '~/Thesis/final_code/example_networks/betweenness/ftg_pclrc_b_B.csv')
write.csv(fntg_pclrc_b_B, '~/Thesis/final_code/example_networks/betweenness/fntg_pclrc_b_B.csv')
write.csv(ftg_corr_k_B, '~/Thesis/final_code/example_networks/betweenness/ftg_corr_k_B.csv')
write.csv(fntg_corr_k_B, '~/Thesis/final_code/example_networks/betweenness/fntg_corr_k_B.csv')
#closeness

ftg_pclrc_b_C <-  closeness(ftg_pclrc_b_ig)
ftg_corr_k_C <- closeness(fntg_pclrc_b_ig)
fntg_pclrc_b_C <- closeness(ftg_corr_k_ig)
fntg_corr_k_C <- closeness(fntg_corr_k_ig)

write.csv(ftg_pclrc_b_C, '~/Thesis/final_code/example_networks/closeness/ftg_pclrc_b_C.csv')
write.csv(fntg_pclrc_b_C, '~/Thesis/final_code/example_networks/closeness/fntg_pclrc_b_C.csv')
write.csv(ftg_corr_k_C, '~/Thesis/final_code/example_networks/closeness/ftg_corr_k_C.csv')
write.csv(fntg_corr_k_C, '~/Thesis/final_code/example_networks/closeness/fntg_corr_k_C.csv')

#eigenvector centrality 
ftg_pclrc_b_E <-  eigen_centrality(ftg_pclrc_b_ig)
ftg_corr_k_E <- eigen_centrality(fntg_pclrc_b_ig)
fntg_pclrc_b_E <- eigen_centrality(ftg_corr_k_ig)
fntg_corr_k_E <- eigen_centrality(fntg_corr_k_ig)

write.csv(ftg_pclrc_b_E$vector, '~/Thesis/final_code/example_networks/eigen/ftg_pclrc_b_E.csv')
write.csv(fntg_pclrc_b_E$vector, '~/Thesis/final_code/example_networks/eigen/fntg_pclrc_b_E.csv')
write.csv(ftg_corr_k_E$vector, '~/Thesis/final_code/example_networks/eigen/ftg_corr_k_E.csv')
write.csv(fntg_corr_k_E$vector, '~/Thesis/final_code/example_networks/eigen/fntg_corr_k_E.csv')


#differential connectivity ---------------------
diff_corr_k <- DifferentialConnectivity(f_tg[,10:208], f_ntg[,10:208], method = 'CORR_k', thresh = 0.62)
diff_corr_k3 <- DifferentialConnectivity(f_tg[1:4,10:208], f_tg[5:9,10:208], method = 'CORR_k', thresh = 0.62)
diff_corr_k2 <- DifferentialConnectivity(f_ntg[,10:208], nf_ntg[,10:208], method = 'CORR_k', thresh = 0.62)

sum(diff_corr_k$PredictedConn)
sum(diff_corr_k2$PredictedConn)
sum(diff_corr_k3$PredictedConn)


test <- DifferentialConnectivity(f_tg[,10:208], nf_tg[,10:208], method = 'PCLRC_b', thresh = 0.62)
sum(test$PredictedConn)


diff_pclrc_b <- DifferentialConnectivity(f_tg[,10:208], f_ntg[,10:208], method = 'PCLRC_p', thresh = 0.64)
diff_pclrc_b3 <- DifferentialConnectivity(f_tg[1:4,10:208], f_tg[5:9,10:208], method = 'PCLRC_b', thresh = 0.64)
diff_pclrc_b2 <- DifferentialConnectivity(f_ntg[,10:208], nf_ntg[,10:208], method = 'PCLRC_b', thresh = 0.64)

sum(diff_pclrc_b$PredictedConn)
sum(diff_pclrc_b2$PredictedConn)
sum(diff_pclrc_b3$PredictedConn)


#export results for cytoscape 
write.csv(diff_pclrc_b,'~/Thesis/final_code/example_networks/DCA/pclrc_dca.csv')
write.csv(diff_corr_k,'~/Thesis/final_code/example_networks/DCA/corr_dca.csv')


#heatmaps
heatmap(ftg_corr_k, Colv = NA, Rowv = NA, scale="none",col = c('#F4F4F4','#423369'))
heatmap(ftg_pclrc_b, Colv = NA, Rowv = NA, scale="none",col = c('#F4F4F4','#423369'))

heatmap(fntg_pclrc_b, Colv = NA, Rowv = NA, scale="none",col = c('#F4F4F4','#423369'))
heatmap(fntg_corr_k, Colv = NA, Rowv = NA, scale="none",col = c('#F4F4F4','#423369'))

#treemaps 
library(treemap)

pclrc_dca <- as.data.frame(matrix(nrow = 3, ncol = 2))
corr_dca <- as.data.frame(matrix(nrow = 3, ncol = 2))

pclrc_dca[,1] <- c('DC = 3', 'N_DC = 119', 'DC_B = 77')
pclrc_dca[,2] <- c(3,119,77)
colnames(pclrc_dca) <- c('group','value')

corr_dca[,1] <- c('DC = 104', 'N_DC = 18', 'DC_B = 77')
corr_dca[,2] <- c(104,18,77)

colnames(corr_dca) <- c('group','value')

treemap(pclrc_dca, index = 'group', vSize = 'value', type = 'index', palette = c('#b2df8a', '#89b562', '#bc80bd'))

treemap(corr_dca, index = 'group', vSize = 'value', type = 'index', palette = c('#b2df8a', '#89b562', '#bc80bd'))
