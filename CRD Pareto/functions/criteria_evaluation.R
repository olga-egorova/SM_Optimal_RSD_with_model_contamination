# Evaluating component criteria: DP, LoF(DP), mse(D)

##### Unwrapping the "Search" object

criteria_values = function(S, eps = 10^(-23)) {
  
  Xp = S$Xp; Xq = S$Xq; Labels = S$Labels;
  Ds = 0; DP = 0; LoFDP = 0; mse = 0;
  Nruns = nrow(Xp); P = ncol(Xp); Q = ncol(Xq)
  DF = nlevels(as.factor(Labels))              
  df = Nruns-DF                                # df - pure error degrees of freedom
  
  M = crossprod(Xp)                            # information matrix of primary terms
  D = prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D > eps) 
  {
    Minv<-solve(M)
  } else {return (list (eval=0));}    # if M is computationally singular
  
  if (D^(1.0/(P-1))>0)
  {
    Ds = (D/Nruns)^(-1.0/(P-1))             
  } else {Ds=0.0}  # Ds 
    
  if (df>0)
  {
    DP = Ds*qf(1-alpha.DP,P-1,df)
  } else { DP=0.0} # if df=0   

    if (df > 0) {
      M12 = crossprod(Xp,Xq)
      M_ = M + diag(max(0, 1./diag(Sigma_beta_p/sigma2)),P)
      M_.inv = solve(M_)
      A = M_.inv%*%M12
      L0 = crossprod(Xq)-t(M12)%*%A+diag(1./diag(Sigma_beta_q/sigma2))
      LoFD = sigma2*(det(L0))^(-1.0/Q)
      LoFDP = LoFD*qf(1-alpha.LoF,Q, df) 
    } else { LoFDP=0.0 } 
                             
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

  return (list (weights = S$weights, compound = S$compound, df = S$df, DP = DP, LoFDP = LoFDP, mse = mse))
    
}

## Constructing the matrix of relative efficiencies

criteria_efficiencies = function(SS) {
  n = length(SS)
  summary = matrix(0, nrow = n, ncol = n)
  weights = list()
  for (i in 1:n) {
    weights[[i]] = SS[[i]]$weights
    for (j in 1:n) {
      if (i==j) {summary[i,j] = 100}
    else{
      k.DP = as.numeric(SS[[j]]$weights[1]); k.LoF = as.numeric(SS[[j]]$weights[2]); 
      k.mse = as.numeric(SS[[j]]$weights[3]);
      Si = criteria_values(SS[[i]])
      summary[i,j] = SS[[j]]$compound/((Si$DP^k.DP)*(Si$LoFDP^k.LoF)*(Si$mse^k.mse))*100
      }
    }
  }
  return (list(relative_efficiency = round(summary, digits = 3), weights = weights))
}

##### Criteria evaluation from primary and potential matrices and labels

criteria_evaluation = function(Labels, Xp, Xq, eps = 10^(-20)){
  
  Ds = 0; DP = 0; LoFDP = 0; mse = 0;
  Nruns = nrow(Xp); P = ncol(Xp); Q = ncol(Xq)
  DF = nlevels(as.factor(Labels))              
  df = Nruns-DF                                # df - pure error degrees of freedom
  
  M = crossprod(Xp)                            # information matrix of primary terms
  D = prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D > eps) 
  {
    Minv<-solve(M)
  } else {return (list (eval=0));}    # if M is computationally singular
  
  if (D^(1.0/(P-1))>0)
  {
    Ds = (D/Nruns)^(-1.0/(P-1))             
  } else {Ds=0.0}  # Ds 
  
  if (df>0)
  {
    DP = Ds*qf(1-alpha.DP,P-1,df)
  } else { DP=0.0} # if df=0   
  
  if (df > 0) {
    M12 = crossprod(Xp,Xq)
    M_ = M + diag(pmax(rep(0,P), 1./diag(Sigma_beta_p/sigma2)))
    M_.inv = solve(M_)
    A = M_.inv%*%M12
    L0 = crossprod(Xq)-t(M12)%*%A+diag(1./diag(Sigma_beta_q/sigma2))
    LoFD = sigma2*(det(L0))^(-1.0/Q)
    LoFDP = LoFD*qf(1-alpha.LoF,Q, df) 
  } else { LoFDP=0.0 } 
  
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
  
  compound = (DP^kappa.DP)*(LoFDP^kappa.LoF)*(mse^kappa.mse)
  
  return (list (weights = c(kappa.DP, kappa.LoF, kappa.mse), compound = compound, 
                df=df, DP = DP, LoFDP = LoFDP, mse = mse))
}

