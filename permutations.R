

permutation_test <- function(df, ss, m1, m2, n = 100, stat, ref_net, d){
  #ARGS:
  # df: data from which to generate networks
  # ss: sample size to be tested 
  # m1,m2: methods to be tested (ex: auroc for CORR_p vs auroc for CORR_s)
  # n = number of permutations to run 
  # p = pval significance cutoff 
  # stat: statistic to be tested
  
  #get vectorized versions of the reference matrices 
   ref_v <- get_vec(ref_net)
  # 
  # 
  # 
  # #get predicted networks
  # pred1 = get_net(df, m = m1)
  # pred2 = get_net(df, m = m2)
  # 
  # 
  # #threshold networks and get stats 
  # stat1 = get_stat(pred1, ref_net, ref_v, st = stat)
  # stat2 = get_stat(pred2, ref_net, ref_v, st = stat)
  
  #get reference difference 
  #diff0 = abs(stat1 - stat2)
  diff0 = d
  # get data of appropriate sample size 
  samples=sample(dim(df)[1], ss)
  df_s <- df[samples,]
  
  DIFFP = c()
  #permutation begins
  r = dim(df_s)[1]
  
  for (ii in 1:n) {
    X1p = NULL
    X2p = NULL
    
    #shuffle samples per col 
    for (jj in 1:dim(df_s)[2])  {
      idx1p = sample(r,r)
      idx2p = sample(r,r)
      X1p = cbind(X1p,df_s[idx1p,jj])
      X2p = cbind(X2p,df_s[idx2p,jj])
    }
    
    # Calculate statistics
    #get predicted networks
    
    colnames(X1p) <- colnames(df_s)
    colnames(X2p) <- colnames(df_s)
    
    
    predX1 = get_net(X1p, m = m1)
    predX2 = get_net(X2p, m = m2)
    
    
    #threshold networks and get stats 
    statx1 = get_stat(predX1, ref_net, ref_v, st = stat)
    statx2 = get_stat(predX2, ref_net, ref_v, st = stat)

    #get permutated difference 
    diffp = abs(statx1 - statx2)
    #print(d)
    #print(diffp)
    DIFFP = append(DIFFP, diffp)
    
    
  }
  
  
  #Calculate Pvalue
  LLC = length(which(DIFFP>=diff0))
  pval = (LLC )/n
  
  return(pval)
  
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

get_fdr <- function(net_t, ref){
  fp <- as.numeric(sum(net_t == 1 & ref == 0))
  tp <- as.numeric(sum(net_t == 1 & ref == 1))
  return(fp/(fp+tp))
  
}

get_ppv <- function(net_t, ref){
  fp <- as.numeric(sum(net_t == 1 & ref == 0))
  tp <- as.numeric(sum(net_t == 1 & ref == 1))
  return(tp/(fp+tp))
  
}

get_mcc <- function(net_v, ref_v){
  x <- NULL
  x <- mcc(preds = net_v, actuals = ref_v)
  return(x)
}

get_pr <- function(net_v, ref_v){
  return(PRAUC(net_v, ref_v))
}

get_net <- function(df, m){
  x <- NULL
  #CORRELATION 
  if(m == 'CORR_p'){
    x <- cor(df, method = 'pearson')
  }
  else if(m == 'CORR_s'){
    x <- cor(df, method = 'spearman')
  }
  else if(m == 'CORR_k'){
    x <- cor(df, method = 'kendall')
  }
  #PCLRC  
  else if(m == 'PCLRC_p'){
    x <- getPCLRC(df, m = 'pearson')
  }
  else if(m == 'PCLRC_s'){
    x <- getPCLRC(df, m = 'spearman')
  }
  else if(m == 'PCLRC_k'){
    x <- getPCLRC(df, m = 'kendall')
  }
  else if(m == 'PCLRC_b'){
    x <- getPCLRC(df, m = 'bicor')
  }
  #bicor 
  else if(m == 'bicor'){
    x <- bicor(df)
  }
  
  else if(m == 'pcor'){
    x <- ggm.estimate.pcor(df)
  }
  
  else if(m == 'MI'){
    x <- mutualInfoAdjacency(df)$MutualInformation
  }

  else if(m == 'clr'){
    x <- clr(mutualInfoAdjacency(df)$MutualInformation)
  }
  else if(m == 'mrnet'){
    x <- mrnet(mutualInfoAdjacency(df)$MutualInformation) 
  }
  
  else if(m == 'GENIE'){
    x <- GENIE3(as.matrix(t(df)), nTrees = 100)
  }
  
  
  
  return(x)
}

#function to apply 'best' threshold to get binary matrix 
apply_thresh <- function(net, thresh){
  x <- net 
  x[abs(net) < thresh] <- 0
  x[abs(net) >= thresh] <- 1
  return(x)
}


get_stat <- function(net, ref_net, ref_v, st = 'auroc'){
  
  #get rocs and reference networks/vectors
  net_v <- get_vec(net)
  pROC = roc(ref_v, net_v, quiet = TRUE)
  
  #get best thresholds for primary reference
  p_thresh <- coords(pROC, 'best', ret = 'threshold')
  
  #apply thresholds to get binary network 
  net_thresh <- apply_thresh(net, p_thresh[1,1])
  net_vt <- get_vec(net_thresh)
  
  
  #get stat to return
  if(st == 'auroc'){
    return(pROC$auc)
  }
  else if(st == 'fdr'){
    return(get_fdr(net_thresh, ref_net))
  }
  else if(st == 'ppv'){
    return(get_ppv(net_thresh, ref_net))
  }
  else if(st == 'mcc'){
    return(get_mcc(net_vt, ref_v))
  }
  else if(st == 'pr'){
    return(get_pr(net_vt, ref_v))
  }
  
  
}
  

