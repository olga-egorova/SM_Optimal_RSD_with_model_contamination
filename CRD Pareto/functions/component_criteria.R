## Individual criteria functions. 

## LoF component criterion.  Checks to be added
LoF_DP = function(Labels, Xp, Xq) {
  # assumes diagonal matrices Sigma
  Q = ncol(Xq)
  M12 = crossprod(Xp,Xq)
  M_ = crossprod(Xp) + diag(pmax(rep(0,P), 1./diag(Sigma_beta_p/sigma2))) # 1st stage -- diffuse prior on beta_p
  M_.inv = solve(M_)
  A = M_.inv%*%M12
  L0 = crossprod(Xq)-t(M12)%*%A+diag(1./diag(Sigma_beta_q/sigma2))
  LoFD = sigma2*(det(L0))^(-1.0/Q)
  
  DF = nlevels(as.factor(Labels))              
  df = Nruns-DF 
  if (df>0)
  { 
    LoFDP = LoFD*qf(1-alpha.LoF,Q, df)               # LoF (DP)   
  } else {return (NaN)}
  
  return(LoFDP)
}

##########################

## MSE component. Checks to be added

MSE_DP = function(Xp, Xq, sigma2, biter = Biter) {
  
  Nruns = nrow(Xp); P = ncol(Xp); Q = ncol(Xq)
  Z0 = diag(1,Nruns)-matrix(1/Nruns,nrow=Nruns,ncol=Nruns) 
  M = crossprod(Xp)
  D = prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  Minv = solve(M)
  
  Tvalue = 0
  M120 = t(Xp[,-1])%*%Z0%*%Xq    # excluding the intercept
  
  for (j in 1:biter)            # MC estimation of the second part of the MSE(D)-component
  {
    beta2 = rnorm(Q, mean=0, sd=sqrt(diag(Sigma_beta_q/sigma2)))
    M12b = M120%*%beta2
    Evalue = t(M12b)%*%(Minv[-1,-1])%*%M12b                          
    Tvalue = Tvalue+log(1+Evalue)                              
  }  
  MC = Tvalue/Biter
  mse = (exp(MC)*Nruns/D)^(1./(P-1))  
  
  return(mse)
}


DP = function(Labels, Xp, Xq) {
  
  DF = nlevels(as.factor(Labels))              
  df = Nruns-DF                                # df - pure error degrees of freedom
  
  M = crossprod(Xp)                       # information matrix of primary terms
  D = prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D>eps) 
  {
    Minv = solve(M)
  } else {return (list (Ds=0, DP=0));} # if M is computationally singular
  
  if (D^(1.0/(P-1))>0)
  {
    D = (D/Nruns)^(-1.0/(P-1))     
  } else {return (list (Ds=0, DP=0));} # Ds
  
  if (df>0)
  {
    DP = Ds*qf(1-alpha.DP,P-1,df)
  }
  return (list(Ds = Ds, DP = DP))
  
}


compound.criteria_d = function(Labels, Xp, Xq, eps = 10^(-23)) {
  
  Ds = 0; DP = 0; LoFDP = 0; mse = 0;
  Nruns = nrow(Xp); P = ncol(Xp); Q = ncol(Xq)
  DF = nlevels(as.factor(Labels))              
  df = Nruns-DF                                # df - pure error degrees of freedom
  
  M = crossprod(Xp)                            # information matrix of primary terms
  D = prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D > eps) 
  {
    Minv<-solve(M)
  } else {return (list (eval=0, compound=10^6));}    # if M is computationally singular
  
  if (kappa.DP>0)  ## DP component
  {
    if (D^(1.0/(P-1))>0)
    {
      Ds = (D/Nruns)^(-1.0/(P-1))             
    } else {return (list (eval=0, compound=10^6));}  # Ds 
    
    if (df>0)
    {
      DP = Ds*qf(1-alpha.DP,P-1,df)
    } else {return (list (eval=0, compound=10^6));} # if df=0   
  }
  
  if (kappa.LoF>0)      # LoF(DP) component
  {
    if (df > 0) {
      M12 = crossprod(Xp,Xq)
      M_ = M + diag(max(0, 1./diag(Sigma_beta_p/sigma2)),P)
      M_.inv = solve(M_)
      A = M_.inv%*%M12
      L0 = crossprod(Xq)-t(M12)%*%A+diag(1./diag(Sigma_beta_q/sigma2))
      LoFD = sigma2*(det(L0))^(-1.0/Q)
      LoFDP = LoFD*qf(1-alpha.LoF,Q, df) 
    } else {return (list (eval=0, compound=10^6));} 
  }                           
  
  if (kappa.mse>0)     # MSE(D) component
  {
    Tvalue = 0
    Z0 = diag(1,Nruns)-matrix(1/Nruns,nrow=Nruns,ncol=Nruns) 
    M120 = t(Xp[,-1])%*%Z0%*%Xq    # excluding the intercept
    
    for (j in 1:Biter){            # MC estimation of the second part of the MSE(D)-component{
      beta2 = rnorm(Q, mean=0, sd=sqrt(diag(Sigma_beta_q/sigma2)))
      M12b = M120%*%beta2
      Evalue = t(M12b)%*%(Minv[-1,-1])%*%M12b                          
      Tvalue = Tvalue+log(1+Evalue)                              
    }  
    MC = Tvalue/Biter
    mse = (exp(MC)*Nruns/D)^(1./(P-1))  
  }
  
  compound = (DP^kappa.DP)*(LoFDP^kappa.LoF)*(mse^kappa.mse)
  list (eval=1, Ds = Ds, DP=DP, LoFDP=LoFDP, mse=mse, df=df, compound=compound)
  
}

