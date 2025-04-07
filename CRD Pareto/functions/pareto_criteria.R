## Individual criteria evaluation functions.  

## LoF component criterion

LoF_DP = function(Labels, Xp, Xq, eps = 10^(-23)) {
  # assumes diagonal matrices Sigma
  Nruns = nrow(Xp); P = ncol(Xp); Q = ncol(Xq)
  M12 = crossprod(Xp,Xq)
  M_ = crossprod(Xp) 
  D_ = prod(round(eigen(M_,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D_ > eps){
    M_.inv = solve(M_)
  } else {return (NaN)}
  A = M_.inv%*%M12
  
  L0 = crossprod(Xq)-t(M12)%*%A+diag(1./tau2, nrow = Q)
  LoFD = (det(L0))^(-1.0/Q)
  
  DF = nlevels(as.factor(Labels))              
  df = Nruns-DF 
  if (df>0)
  { 
    LoFDP = LoFD*qf(1-alpha.LoF,Q, df)               # LoF (DP)   
  } else {return (NaN)}
  
  return (LoFDP)
}

####### MSE component  ####### 

MSE_DP = function(Xp, Xq, beta_q, biter = Biter, eps = 10^(-23)) {
  
  nruns = nrow(Xp); P = ncol(Xp); Q = ncol(Xq)
  Z0 = diag(1, nrow = nruns)-matrix(1/nruns, nrow=nruns, ncol=nruns) 
  M = crossprod(Xp)
  D = prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D > eps) 
  {
    Minv = solve(M)
  } else {return (NaN);}          # if M is computationally singular
  
  Evalue = rep(0, biter)
  M120 = t(Xp[,-1])%*%Z0%*%Xq    # excluding the intercept
  M12b = t(M120)%*%Minv[-1,-1]%*%M120  # Q x Q matrix for the final quadratic form
  
  MC =  mean(log(1 + colSums(beta_q*(M12b%*%beta_q)))) # MC estimation of the second part of the MSE(D)-component
  mse = (exp(MC)*Nruns/D)^(1./(P-1))  
  
  return(list(mse = mse))
}

####### DPs component ##########

DPs = function(Labels, Xp, eps = 10^(-23)) {
  
  Nruns = nrow(Xp); P = ncol(Xp);
  DF = nlevels(as.factor(Labels))              
  df = Nruns-DF                                # df - pure error degrees of freedom
  
  M = crossprod(Xp)                            # information matrix of primary terms
  D = prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D > eps) 
  {
    Minv = solve(M)
  } else {return (list (Ds=NaN, DP=NaN));}          # if M is computationally singular
  
  if (D^(1.0/(P-1))>0)
  {
    Ds = (D/Nruns)^(-1.0/(P-1))     
  } else {return (list (Ds=NaN, DP=NaN));} # Ds
  
  if (df>0)
  {
    DP = Ds*qf(1-alpha.DP,P-1,df)
  } else {DP = NaN}
  return (list(Ds = Ds, DP = DP))
}


#### Obtain a vector of criteria for constructing Pareto frontier #### 

pareto_criteria = function(Labels, Xp, Xq, beta_q, eps = 10^(-23)){

  DP = DPs(Labels, Xp)$DP
  LoF_DP = LoF_DP(Labels, Xp, Xq)
  mse = MSE_DP(Xp, Xq, beta_q)
 
  return (list(values = list(DP=DP, LoF_DP=LoF_DP, mse = mse$mse)))

}


####### Design domination #######################################

### The one with three components
is.dominated = function(p.crit1, p.crit2){

  p.1 = unlist(p.crit1$values); 
  p.2 = unlist(p.crit2$values)

  if (any(is.nan(p.2))) {return (FALSE)}
  if (any(is.nan(p.1))) {return (TRUE)}
  
  
  return ((prod(c(p.1[1:3] >= p.2[1:3])) == 1) && 
             (sum(c(p.1[1:3] > p.2[1:3])) > 0))     # p.crit2 is not worse, and for at least one  - better
}

##### Update the existing list of non-dominated designs  ###############

update.nd_list = function (list.nd, Labels, p.crit) {
  #browser()
  n = length(list.nd); 
  if (n==0) {                     # if this is an empty list -- put in there the first design
    return (list(l = list(list(labels = Labels, p.crit = p.crit)), added = TRUE))
  }
  removed_ind = c();
  removed = FALSE; 
  for (curr in 1:n){
    nd.current = list.nd[[curr]]
    if ((is.dominated(p.crit1 = p.crit, p.crit2 = nd.current$p.crit)) || 
        identical(sort(nd.current$labels), sort(Labels))){    # the new one is dominated or it is the same design
      return (list(l = list.nd, added = FALSE))                    ## no changes performed
    } else {
      if (is.dominated(p.crit1 = nd.current$p.crit, p.crit2 = p.crit)){  # if the new one is dominating
        removed_ind = c(removed_ind, curr)       # the dominated one is to be removed
        removed = TRUE                           
      }
    }
  }
  
  list.nd = append(list.nd, list(list(labels = Labels, p.crit = p.crit)))  # adding the new one to the end of the list
  if (removed) {list.nd[removed_ind] = NULL}  # if there has been any dominated designs
  
  return (list(l = list.nd, added = TRUE))
}
  
  
