##### Obtaining criteria values for the designs

criteria.values.PB.mse.b<-function(X1,X2,eps=10^-23,Biter=1000)    
{
  DP<-0; LoFDP<-0; LP<-0; LoFLP<-0;
  df<-df.residual(aov(y~b+as.factor(X1[,1])))
  DF<-Nruns-df
  tau<-tau2^0.5
  
  M<-t(X1[,-1])%*%Z%*%X1[,-1]                 # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))          #/Nruns
  if (D>eps) 
  {
    Minv<-solve(M)
  } else {return (eval=0)}
  if (D^(1.0/P)>0)
  {
    Ds<-D^(-1.0/P)                        # Ds component
  } else {return (eval=0)}
  
  #Ls<-W%*%(diag(Minv))                    # Ls component
  
  X1.b<-cbind(B,X1[,-1])                      # adding columns corresponding to blocks, removing labels
  X2.b<-X2[,-1]
  M12.b<-crossprod(X1.b,X2.b)
  M.b<-crossprod(X1.b)
  Mb.inv<-solve(M.b)
  A<-Mb.inv%*%M12.b
  L0<-crossprod(X2.b)-t(M12.b)%*%A+diag(1./tau2,nrow=Q)
  LoFD<-(det(L0))^(-1.0/Q)                                          # LoF (D)
  #LoFL<-1./(sum(diag(L0))/Q)                                        # LoF (L)
  
  if (df>0)
  {
    DP<-Ds*qf(1-alpha.DP,P,df)                                 # DP component
    #LP<-Ls*qf(1-alpha.LP,1,df)                                 # LP component
    LoFDP<-LoFD*qf(1-alpha.LoF,Q,df)                           # LoF (DP)   
    L0.inv.trace<-sum(1./eigen(L0,only.values=TRUE)$values)
    #LoFLP<-L0.inv.trace*qf(1-alpha.LoFL,1,df)/Q                # LoF (LP)
  } else {DP<-0;LoFDP<0; LoFLP<-0}
  
  A22<-A[(Nblocks+1):(Nblocks+P),]
  Mb.inv22<-Mb.inv[(Nblocks+1):(Nblocks+P),(Nblocks+1):(Nblocks+P)]	
  #mseL<-sum(diag(Mb.inv22+tau2*tcrossprod(A22)))/P                   # MSE(L)

  Tvalue<-0
  M12<-t(X1[,-1])%*%Z%*%X2.b
  MM<-t(M12)%*%Minv%*%M12

  for (j in 1:Biter)                # random beta-s, MC for the expectation
  {
      beta2<-rnorm(Q,mean=0,sd=tau)
      Evalue<-t(beta2)%*%MM%*%beta2                         
      Tvalue<-Tvalue+log(1+Evalue)                              
  }  
  MC<-Tvalue/Biter
  mseD<-(exp(MC)/D)^(1./P)                                    # MSE(D)   
  
  beta2<-rep(tau,Q) 
  Tvalue<-(1+t(beta2)%*%MM%*%beta2) 
  mseD_point <- (Tvalue/D)^(1./P)   
    
  list (DP=DP, LoFDP=LoFDP, mseD=mseD, mseD_point = mseD_point)
}

### With full LoF component (DPs optimality for potential terms)

criteria.values.full.mse.b<-function(X1,X2,eps=10^-23,Biter=1000)      
{
  DP<-0; LoFDP<-0; LP<-0; LoFLP<-0;
  df<-df.residual(aov(y~b+as.factor(X1[,1])))
  DF<-Nruns-df
  tau<-tau2^0.5
  
  M<-t(X1[,-1])%*%Z%*%X1[,-1]                 # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))          #/Nruns
  if (D>eps) 
  {
    Minv<-solve(M)
  } else {return (eval=0)}
  if (D^(1.0/P)>0)
  {
    Ds<-D^(-1.0/P)                        # Ds component
  } else {return (eval=0)}
  
  Ls<-W%*%(diag(Minv))                    # Ls component
  
  X1.b<-cbind(B,X1[,-1])                      # adding columns corresponding to bloks, removing labels
  X2.b<-X2[,-1]
  M12.b<-crossprod(X1.b,X2.b)
  M.b<-crossprod(X1.b)
  Mb.inv<-solve(M.b)
  A<-Mb.inv%*%M12.b
  L0<-crossprod(X2.b)-t(M12.b)%*%A         # +diag(1./tau2,nrow=Q)
  LoFD<-(det(L0))^(-1.0/Q)                                          # LoF (D)
  LoFL<-1./(sum(diag(L0))/Q)                                        # LoF (L)
  
  if (df>0)
  {
    DP<-Ds*qf(1-alpha.DP,P,df)                                 # DP component
    LP<-Ls*qf(1-alpha.LP,1,df)                                 # LP component
    LoFDP<-LoFD           #*qf(1-alpha.LoF,Q,df)                           # LoF (DP)   
    L0.inv.trace<-sum(1./eigen(L0,only.values=TRUE)$values)
    LoFLP<-L0.inv.trace   #*qf(1-alpha.LoFL,1,df)/Q                # LoF (LP)
  } else {DP<-0;LoFDP<0; LoFLP<-0}
  
  A22<-A[(Nblocks+1):(Nblocks+P),]
  Mb.inv22<-Mb.inv[(Nblocks+1):(Nblocks+P),(Nblocks+1):(Nblocks+P)]  
  mseL<-sum(diag(Mb.inv22+tau2*tcrossprod(A22)))/P                   # MSE(L)
  
  Tvalue<-0
  M12<-t(X1[,-1])%*%Z%*%X2.b
  MM<-t(M12)%*%Minv%*%M12
  
  for (j in 1:Biter)                # random beta-s, MC for the expectation
  {
    beta2<-rnorm(Q,mean=0,sd=tau)
    Evalue<-t(beta2)%*%MM%*%beta2                         
    Tvalue<-Tvalue+log(1+Evalue)                              
  }  
  MC<-Tvalue/Biter
  mseD<-(exp(MC)/D)^(1./P)                                    # MSE(D)   
  
  list (DP=DP, LoFDP=LoFDP, mseD=mseD, LP=LP, LoFLP=LoFLP, mseL=mseL)
}


output<-function(S.GD_PB)
{
  list1<-c(S.GD_PB$df, S.GD_PB$compound, S.GD_PB$mse, S.GD_PB$time)
  list(out=list1, values=criteria.values.PB.mse.b(S.GD_PB$X1, S.GD_PB$X2))
}

### Plotting the designs

# Obtaining pure design matrix, no lables, transformed to [-1,1]
extract.design<-function(X1)
{
  labels<-as.vector(X1[,1])                            # extract labels
  cand<-as.matrix(expand.grid(Levels))
  candl<-cbind(label.b(cand),cand)            # labelling
  index<-rep(0,Nruns)
  for (i in 1:Nruns)
  {
    index[i]<-which(candl[,1]==labels[i])
  }
  design<-apply(cand[index,],2,transform.b)
}


