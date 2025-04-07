# Forms the labelled candidate set of treatments

candidate_trt_set <- function(Levels, K, Hypercube = TRUE)
{
  cand <- as.matrix(expand.grid(Levels))
  candl <- cbind(label(cand, Levels, K), cand)                         # creating a column of treatment labels
  cand.trt <- cbind(candl[,1], apply(as.matrix(candl[,-1]), 2, Transform))   # rescaling the factors' values to [-1, 1]
  if (!Hypercube) {cand.trt<-as.matrix(spheric(cand.trt, K))}     # spheric coords
  return (cand.trt)
}



# Forms the full candidate set of treatments polynomial terms
# cand.trt -- candidate set of treatments, the first column contains treatment labels. 
candidate_set_full = function(cand.trt, K) {
  
  cand.terms = cand.trt[, -1, drop = F]  #  linear terms
  
  ### Naming linear and generating and naming quadratic terms
  for (k in 1:K) {
    colnames(cand.terms)[k] = paste("x", as.character(k), sep = "")  # named linear terms, e.g. "x3"
    cand.terms = cbind(cand.terms, cand.terms[,k]^2)
    colnames(cand.terms)[K+k] = paste("x", as.character(k), as.character(2), sep = "") # quadratic terms: "x12", "x22", "x32", etc.
  }
  
  ### Interaction terms: linear by linear
  if (K > 1) {
    int.count = 0
    for (k in 1:(K-1)){
      for (j in (k+1):K) {
        int.count = int.count + 1                  # count the number of interaction terms
        cand.terms = cbind(cand.terms, cand.terms[,k]* cand.terms[,j])
        colnames(cand.terms)[2*K + int.count] = 
          paste("x", as.character(k), "x", as.character(j), sep = "")   # interaction terms: "x1x2", "x2x3", "x3x4", etc.
      }
    }
    
    ### Interaction terms: quadratic by linear 
    int2.count = 0
    nterms = ncol(cand.terms)
    for (k in 1:(K-1)){
      for (j in (k+1):K) {
        int2.count = int2.count + 2                # count the number of QxL interaction terms
        cand.terms = cbind(cand.terms, cand.terms[,k]^2 * cand.terms[,j], cand.terms[,k] * cand.terms[,j]^2)
        colnames(cand.terms)[nterms + int2.count - c(1,0)] = 
          c(paste("x", as.character(k), as.character(2), "x", as.character(j), sep = ""),
            paste("x", as.character(k), "x", as.character(j), as.character(2), sep = "")) # QxL  and LxQ interaction terms: "x12x2", "x1x22", etc.
      }
    }
  }
  
  ### Linear-by-linear-by-linear terms
  if (K > 2) {
    int3.count = 0
    nterms = ncol(cand.terms)
    for (k in 1:(K-2)){
      for (j in (k+1):(K-1)) {
        for (i in (j+1):K){
          int3.count = int3.count + 1                     # count the number of LxLxL interaction terms
          cand.terms = cbind(cand.terms, cand.terms[,k] * cand.terms[,j]* cand.terms[,i])
          colnames(cand.terms)[nterms + int3.count] = 
            paste("x", as.character(k), "x", as.character(j), "x", as.character(i), sep = "")   # LxLxL terms: "x1x2x3", "x1x3x4", etc.
        }
        
      }
    }
  }
  
  ### Cubic terms
  nterms = ncol(cand.terms)
  for (k in 1:K){
    cand.terms = cbind(cand.terms, cand.terms[,k]^3)
    colnames(cand.terms)[nterms + k] = paste("x", as.character(k) , as.character(3), sep = "")   # cubic terms: "x13", "x33", etc.
  }
  
  ### Fourth order terms
  
  if (K>1){
    nterms = ncol(cand.terms); int4.count = 0
    for (k in 1:(K-1)){
      for (j in (k+1):K) { 
        # Quadratic x Quadratic
        cand.terms = cbind(cand.terms, cand.terms[,k]^2*cand.terms[,j]^2)
        int4.count = int4.count + 1
        colnames(cand.terms)[nterms + int4.count] = paste0("x", as.character(k), as.character(2),
                                                    "x", as.character(j), as.character(2))
        # Cubic x Linear
        cand.terms = cbind(cand.terms, cand.terms[,k]^3*cand.terms[,j], 
                           cand.terms[,k]*cand.terms[,j]^3)
        int4.count = int4.count + 2
        colnames(cand.terms)[nterms + int4.count - c(1,0)] = 
          c(paste("x", as.character(k), as.character(3), "x", as.character(j), sep = ""),
            paste("x", as.character(k), "x", as.character(j), as.character(3), sep = ""))
      }
    }
  }
  # QxLxL, LxQxL and LxLxQ terms: "x12x2x3", "x1x22x4", etc.
  if (K > 2) {
    nterms = ncol(cand.terms); int5.count = 0
    for (k in 1:(K-2)){
      for (j in (k+1):(K-1)) {
        for (i in (j+1):K){  
          cand.terms = cbind(cand.terms, cand.terms[,k]^2*cand.terms[,j]*cand.terms[,i], 
                             cand.terms[,k]*cand.terms[,j]^2*cand.terms[,i],
                             cand.terms[,k]*cand.terms[,j]*cand.terms[,i]^2)
          int5 = int5.count + 3
          colnames(cand.terms)[nterms + int5 - c(2,1,0)] = 
            c(paste("x", as.character(k), as.character(2), "x", as.character(j), "x", as.character(i), sep = ""),
              paste("x", as.character(k), "x", as.character(j), as.character(2), "x", as.character(i), sep = ""),
              paste("x", as.character(k), "x", as.character(j), "x", as.character(i), as.character(2), sep = ""))
        }
      }
    }
  }
  
  if (K > 3) {
    nterms = ncol(cand.terms); int6.count = 0
    for (k in 1:(K-3)){
      for (j in (k+1):(K-2)) {
        for (i in (j+1):(K-1)){
          for (l in (i+1):K){   # LxLxLxL terms: "x1x2x3x4", etc.
            cand.terms = cbind(cand.terms, 
                               cand.terms[,k]*cand.terms[,j]*cand.terms[,i]*cand.terms[,i])
            int6.count = int6.count + 1
            colnames(cand.terms)[nterms + int6.count] = 
              paste("x", as.character(k), "x", as.character(j),
                    "x", as.character(i), "x", as.character(l),sep = "")
          }
        }
      }
    }
  }
  
  nterms = ncol(cand.terms)
  for (k in 1:K){
    cand.terms = cbind(cand.terms, cand.terms[,k]^4)
    colnames(cand.terms)[nterms + k] = paste0("x", as.character(k), as.character(4))
  }
  
  #########################################################################
  
  cand.terms = cbind(cand.trt[,1], 1, cand.terms)    # append the column of labels and intercept column
  colnames(cand.terms)[1:2] <- c("label", "intercept")
  
  return (cand.terms)
}
