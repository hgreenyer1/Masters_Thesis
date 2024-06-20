#function to calculate the AUROC for a given reference matrix and range of sample sizes 
library(pROC) #for calculating the auroc 
library(ggm) #pcor

source('~/Thesis/final_code/inference_a/PCLRC_CORR.R') 
#source('~/Thesis/final_code/inference_a/LOPC_R.R') #for partial correlations 
source('~/Thesis/final_code/inference_a/CLRCR.R') 
source('~/Thesis/final_code/inference_a/test_method.R') 

#function to test performance of algorithms for different sample sizes 
get_auc <- function(df, p_ref, method){
  #ARGS: 
  # df: experimental data
  # _ref: primary, secondary, and tertiary 'real' networks (binary)
  # method: name referencing selected network inference algorithm
  
  
  s <- c(5,10,20,50,100,200)
  #s <- c(5)
  
  stats <- c('sample_size','runtime','p_AUROC', 'best_thresh', 
             'p_fdr', 'p_ppv','p_mcc','alpha','KS', 'KS_p')
  
  out <- matrix(nrow = length(s), ncol = length(stats))
  
  out[,1] <- s
  
  
  #vectorize reference matrices 
  p_v <- get_vec(p_ref)
  # s_v <- get_vec(s_ref)
  # t_v <- get_vec(t_ref)
   
  r <- 1
  #loop through globally declared vector of sample sizes
  for(i in s){
    
    #print(i)
    
    #randomly sample from simulated data
    samples=sample(dim(df)[1], i)
    df_s <- df[samples,]
    
    #get predicted network from selected method
    #CORRELATION 
    if(method == 'CORR_p'){
      ptm <- proc.time()
      pred <- cor(df_s, method = 'pearson')
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
      
    }
    else if(method == 'CORR_s'){
      ptm <- proc.time()
      pred <- cor(df_s, method = 'spearman')
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
      
    }
    else if(method == 'CORR_k'){
      ptm <- proc.time()
      pred <- cor(df_s, method = 'kendall')
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
      
    }
    #PARTIAL CORRELATION 
    if(method == 'pcor'){
      ptm <- proc.time()
      pred <- ggm.estimate.pcor(df_s)
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
      
    }
    
    
    # PCLRC VARIATIONS 
    else if(method == 'PCLRC_p'){
      ptm <- proc.time()
      pred <- getPCLRC(df_s, m = 'pearson')
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
      
    }
    else if(method == 'PCLRC_s'){
      ptm <- proc.time()
      pred <- getPCLRC(df_s, m = 'spearman')
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
      
    }
    else if(method == 'PCLRC_k'){
      ptm <- proc.time()
      pred <- getPCLRC(df_s, m = 'kendall')
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
      
    }
    else if(method == 'PCLRC_b'){
      ptm <- proc.time()
      pred <- getPCLRC(df_s, m = 'bicor')
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
      
    }
    
    
    #MUTUAL INFORMATION
    else if(method == 'MI'){
      ptm <- proc.time()
      pred <- mutualInfoAdjacency(df_s)$MutualInformation 
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
    }
    
    else if(method == 'clr'){
      ptm <- proc.time()
      pred <- clr(mutualInfoAdjacency(df_s)$MutualInformation) 
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
    }
    
    else if(method == 'mrnet'){
      ptm <- proc.time()
      pred <- mrnet(mutualInfoAdjacency(df_s)$MutualInformation) 
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
    }
    
    
    #BICOR CORRELATION 
    else if(method == 'bicor'){
      ptm <- proc.time()
      pred <- bicor(df_s)
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
    }
    
    #GENIE (DECISION TREE BASED)
    else if(method == 'GENIE'){
      ptm <- proc.time()
      pred <- GENIE3(as.matrix(t(df_s)), nTrees = 100)
      tm <- proc.time() - ptm
      out[r,2] <- tm[3]
      out <- get_stats(pred, p_ref, p_v, r, out)
      r <- r+1
      next
    }

  }
  colnames(out) <- stats
  return(as.data.frame(out))
}



#function to vectorize matrix and remove redundant entries
get_vec <- function(d){
  #set repeat entries to NA
  d[upper.tri(d, diag = TRUE)] <- NA
  #vectorize
  d_v <- c(d)
  #remove NA entries 
  d_out <- d_v[!is.na(d_v)]
  return(d_out)
}

#function to apply 'best' threshold to get binary matrix 
apply_thresh <- function(net, thresh){
  x <- net 
  x[abs(net) < thresh] <- 0
  x[abs(net) >= thresh] <- 1
  return(x)
}


#set of functions to calculate tpr, fpr, tnr, fnr, fdr, ppv
#ARGS:
#net: binary (thresholded) predicted network 
#ref: one class reference network (ex: primary)

get_fdr <- function(net_t, ref){
  fp <- sum(net_t == 1 & ref == 0)
  tp <- sum(net_t == 1 & ref == 1)
  return(fp/(fp+tp))
  
}

get_ppv <- function(net_t, ref){
  fp <- sum(net_t == 1 & ref == 0)
  tp <- sum(net_t == 1 & ref == 1)
  return(tp/(fp+tp))
  
}

get_mcc <- function(net_t, ref){
  x <- mcc(preds = net_t, actuals = ref)
  return(x)
}



#function to calculate all general algorithm performance statistics 
#ARGS:
#net: weighted predicted network 
# _ref: binary reference matrices 
#p,s,t: vectorized primary, secondary, and tertiary reference matrices
#i: row # corresponding to current sample size iteration
get_stats <- function(net, p_ref, p, i, out){
  
  #vectorize predicted network 
  pred <- get_vec(net)
  
  #roc objects
  p_roc <- roc(p, pred, quiet = TRUE)
  
  #print(i)
  
  #auroc
  out[i,3] <- p_roc$auc
  
  #get best thresholds for primary reference
  p_thresh <- coords(p_roc, 'best', ret = 'threshold')
  out[i,4] <- p_thresh[1,1]
  
  #apply thresholds to get binary network 
  net_thresh <- apply_thresh(net, p_thresh[1,1])
  net_v <- get_vec(net_thresh)
  
  #get degree distribtion and test against power law 
  dg <- rowSums(net_thresh)
  pl <- fit_power_law(dg)
  
  out[i,5] <- get_fdr(net_thresh, p_ref)
  
  out[i,6] <- get_ppv(net_thresh, p_ref)
  
  out[i,7] <- get_mcc(net_v, p)

  out[i,8] <- pl$alpha
  
  out[i,9] <- pl$KS.stat
  
  out[i,10] <- pl$KS.p
    
  
  return(out)
  
}


