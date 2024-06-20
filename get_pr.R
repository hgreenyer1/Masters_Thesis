#function to calculate the AUROC for a given reference matrix and range of sample sizes 
library(pROC) #for calculating the auroc 
library(ggm) #pcor

source('~/Thesis/final_code/inference_a/PCLRC_CORR.R') 
#source('~/Thesis/final_code/inference_a/LOPC_R.R') #for partial correlations 
source('~/Thesis/final_code/inference_a/CLRCR.R') 
source('~/Thesis/final_code/inference_a/test_method.R') 

#function to test performance of algorithms for different sample sizes 
get_pr <- function(df, p_ref, method){
  #ARGS: 
  # df: experimental data
  # _ref: primary, secondary, and tertiary 'real' networks (binary)
  # method: name referencing selected network inference algorithm
  
  
  s <- c(5,10,20,50,100,200)
  #s <- c(5)
  
  stats <- c('sample_size','PR')
  
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
      pred <- cor(df_s, method = 'pearson')
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
      
    }
    else if(method == 'CORR_s'){
      pred <- cor(df_s, method = 'spearman')
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
      
    }
    else if(method == 'CORR_k'){
      ptm <- proc.time()
      pred <- cor(df_s, method = 'kendall')
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
      
    }
    #PARTIAL CORRELATION 
    if(method == 'pcor'){
      pred <- ggm.estimate.pcor(df_s)
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
      
    }
    
    
    # PCLRC VARIATIONS 
    else if(method == 'PCLRC_p'){
      pred <- getPCLRC(df_s, m = 'pearson')
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
      
    }
    else if(method == 'PCLRC_s'){
      pred <- getPCLRC(df_s, m = 'spearman')
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
      
    }
    else if(method == 'PCLRC_k'){
      pred <- getPCLRC(df_s, m = 'kendall')
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
      
    }
    else if(method == 'PCLRC_b'){
      pred <- getPCLRC(df_s, m = 'bicor')
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
      
    }
    
    
    #MUTUAL INFORMATION
    else if(method == 'MI'){
      pred <- mutualInfoAdjacency(df_s)$MutualInformation 
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
    }
    
    else if(method == 'clr'){
      pred <- clr(mutualInfoAdjacency(df_s)$MutualInformation) 
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
    }
    
    else if(method == 'mrnet'){
      pred <- mrnet(mutualInfoAdjacency(df_s)$MutualInformation) 
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
    }
    
    
    #BICOR CORRELATION 
    else if(method == 'bicor'){
      pred <- bicor(df_s)
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
      r <- r+1
      next
    }
    
    #GENIE (DECISION TREE BASED)
    else if(method == 'GENIE'){
      pred <- GENIE3(as.matrix(t(df_s)), nTrees = 100)
      v <- get_vec(pred)
      out[r,2] <- PRAUC(v,p_v)
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
