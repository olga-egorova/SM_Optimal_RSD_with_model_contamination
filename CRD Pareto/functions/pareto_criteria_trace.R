## Individual criteria evaluation functions.  

## LoF component criterion, trace based 

LoF_LP = function(Labels, Xp, Xq, eps = 10^(-23)) {
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
  LoFL = Re(sum(1./eigen(L0,only.values=TRUE)$values))
  
  DF = nlevels(as.factor(Labels))              
  df = Nruns-DF 
  if (df>0)
  { 
    LoFLP = LoFL*qf(1-alpha.LoFL,1,df)/Q               # LoF (LP)   
  } else {return (NaN)}
  
  return (LoFLP)
}

##########################

## MSE component. 

MSE_LP = function(Xp, Xq, eps = 10^(-23)) {
  
  nruns = nrow(Xp); P = ncol(Xp); Q = ncol(Xq)
  M = crossprod(Xp)
  D = prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D > eps) 
  {
    Minv = solve(M)
  } else {return (NaN);}          # if M is computationally singular
  
  M12 = crossprod(Xp,Xq)
  A = Minv%*%M12 
  A0 = Minv+tau2*tcrossprod(A)
  mse = sum(diag(A0)[-1])/(P-1) 
  
  return(list(mse = mse))
}


LPs = function(Labels, Xp, eps = 10^(-23)) {
  
  Nruns = nrow(Xp); P = ncol(Xp);
  DF = nlevels(as.factor(Labels))              
  df = Nruns-DF                                # df - pure error degrees of freedom
  
  M = crossprod(Xp)                            # information matrix of primary terms
  D = prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D > eps) 
  {
    Minv = solve(M)
  } else {return (list (Ls=NaN, LP=NaN));}          # if M is computationally singular
  
  Ls<-W%*%(diag(Minv)[-1])
  
  if (df>0)
  {
    LP = Ls*qf(1-alpha.LP,1,df)
  } else {LP = NaN}
  return (list(Ls = Ls, LP = LP))
}


#### Obtain a vector of criteria for constructing Pareto frontier #### 

pareto_criteria_trace = function(Labels, Xp, Xq, eps = 10^(-23)){

  LP = LPs(Labels, Xp)$LP
  LoF_LP = LoF_LP(Labels, Xp, Xq)
  mse = MSE_LP(Xp, Xq)
 
  return (list(values = list(LP=LP, LoF_LP=LoF_LP, mse = mse$mse)))

}


####### Design domination #######################################

### The one with three components
is.dominated = function(p.crit1, p.crit2){

  p.1 = unlist(p.crit1$values); p.2 = unlist(p.crit2$values)

  if (any(is.nan(p.2))) {return (FALSE)}
  if (any(is.nan(p.1))) {return (TRUE)}
  
  return ((prod(c(p.1[1:3] >= p.2[1:3])) == 1) && 
            (sum(c(p.1[1:3] > p.2[1:3])) > 0))     # p.crit2 is not worse and better at at least one
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
  
  
