#Haley Greenyer 
#edited 2022-02-03

#functions for extracting various 'connectivity' matrices from a stoichiometric matrix (rows = reactions, col = metabolites)


#takes stoichiometric matrix where columns are species and rows are reactions 
get_primary <- function(s, species){
  primary <- matrix(nrow = ncol(s), ncol = ncol(s) )
  primary[,] <- 0
  #for each reaction
  for (k in 1:nrow(s)){
    #get indices of reactant and product
    r <- which(s[k,] == -1)
    p <- which(s[k,] == 1)
    
    #continue if no product or reactant
    if(is_empty(r) | is_empty(p)){
      next
    }
    else{
      primary[r,p] <- 1
      primary[p,r] <- 1
    }
  }
  diag(primary)<-1
  colnames(primary) <- species
  rownames(primary) <- species
  return(primary)
}

#takes in primary matrix and stoichiomatric matrix 
get_secondary <- function(s, primary){
  secondary <- matrix(nrow = ncol(s), ncol = ncol(s) )
  secondary[,] <- 0
  #for each reaction
  for (k in 1:nrow(s)){
    
    #get indices of reactant and product
    r <- which(s[k,] == -1)
    p <- which(s[k,] == 1)
    
    #continue if no product or reactant
    if(is_empty(p) | is_empty(r)){
      next
    }
    else{
      #get rows where product of this reaction is a reactant 
      sec <- which(s[,p] == -1)
      if(is_empty(sec)){
        next
      }
      else{
        #get products (pp) for each reaction where p is a reactant
        for(i in sec){
          pp <- which(s[i,] == 1)
          if(is_empty(pp)){
            next
          }
          else{
            secondary[r,pp] <- 1
            secondary[pp,r] <- 1
          }
        }
      }
    }
  }
  #include links from primary
  diag(secondary)<-1
  secondary[primary == 1] <- 1
  
  #labels 
  colnames(secondary) <- colnames(primary)
  rownames(secondary) <- rownames(primary)
  
  return(secondary)
}

#takes in secondary matrix and stoichiomatric matrix 
get_tertiary <- function(s, secondary){
  tertiary <- matrix(nrow = ncol(s), ncol = ncol(s) )
  tertiary[,] <- 0
  #for each reaction
  for (k in 1:nrow(s)){
    
    #get indices of reactant and product
    r <- which(s[k,] == -1)
    p <- which(s[k,] == 1)
    
    #continue if no product or reactant
    if(is_empty(p) | is_empty(r)){
      next
    }
    else{
      #get rows where product of this reaction (p) is a reactant 
      sec <- which(s[,p] == -1)
      if(is_empty(sec)){
        next
      }
      else{
        #get products (pp) of secondary reactions
        for(i in sec){
          pp <- which(s[i,] == 1)
          if(is_empty(pp)){
            next
          }
          #get rows where secondary products (pp) are reactants 
          else{
            ter <- which(s[,pp] == -1)
            if(is_empty(ter)){
              next
            }
            else{
              #get tertiary interactions (3 step products)
              for(j in ter){
                ppp <- which(s[j,] == 1)
                if(is_empty(ppp)){
                  next
                }
                else{
                  tertiary[r,ppp] <- 1
                  tertiary[ppp,r] <- 1
                }
              }
            }
          }
        }
      }
    }
  }
  #include all links from secondary
  diag(tertiary)<-1
  tertiary[secondary == 1] <- 1
  
  #labels 
  colnames(tertiary) <- colnames(secondary)
  rownames(tertiary) <- rownames(secondary)
  
  return(tertiary)
}




