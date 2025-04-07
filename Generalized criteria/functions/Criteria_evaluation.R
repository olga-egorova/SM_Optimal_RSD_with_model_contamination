## Functions calculating different criteria values

candidate_set_full<-function(cand)          # full candidate matrix, with labels
{
  cand.full<-cbind(cand,potential.matrix(cand)[,-1])
  return (cand.full)
}

criteria.values.mse<-function(X1,X2,eps=10^-23, Biter=10^4)      
  # X1,X2 - primary and potential matrices, with treatment labels
{
  Nruns <- nrow(X1)
  Z0<-diag(1,Nruns)-matrix(1/Nruns,nrow=Nruns,ncol=Nruns)  
  
  DP<-0; LoFDP<-0; LP<-0; LoFLP<-0;
  DF<-nlevels(as.factor(X1[,1]))              
  df<-Nruns-DF                                # df - pure error degrees of freedom
  tau<-tau2^0.5
  
  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))          
  if (D>eps) 
  {
    Minv<-solve(M)
  } else {return (eval=0)}
  if (D^(1.0/(P-1))>0)
  {
    Ds<-(D/Nruns)^(-1.0/(P-1))                        # Ds component
  } else {return (eval=0)}
  
  Ls<-W%*%(diag(Minv)[-1])                            # Ls component
  
  M12<-crossprod(X1[,-1],X2[,-1])
  A<-Minv%*%M12
  L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)
  LoFD<-(det(L0))^(-1.0/Q)                                          # LoF (D)
  LoFL<-1./(sum(diag(L0))/Q)                                        # LoF (L)
  
  if (df>0)
  {
    DP<-Ds*qf(1-alpha.DP,P-1,df)                               # DPs component
    LP<-Ls*qf(1-alpha.LP,1,df)                                 # LPs component
    LoFDP<-LoFD*qf(1-alpha.LoF,Q,df)                           # LoF (DP)   
    L0.inv.trace<-sum(1./eigen(L0,only.values=TRUE)$values)
    LoFLP<-L0.inv.trace*qf(1-alpha.LoFL,1,df)/Q                # LoF (LP)
  } else {DP<-0;LoFDP<-0; LoFLP<-0}
  
  mseL<-sum(diag(Minv+tau2*tcrossprod(A))[-1])/(P-1)                   # MSE(Ls)
  
  Tvalue<-0
  M120<-t(X1[,-c(1,2)])%*%Z0%*%X2[,-1]
  for (j in 1:Biter)                              #  MC for the expectation in the MSE(Ds) component
  {
    beta2<-rnorm(Q,mean=0,sd=tau)
    M12b<-M120%*%beta2
    Evalue<-t(M12b)%*%(Minv[-1,-1])%*%M12b                          
    Tvalue<-Tvalue+log(1+Evalue)                              
  }  
  MC<-Tvalue/Biter
  mseD<-(exp(MC)*Nruns/D)^(1./(P-1))                           # MSE(Ds)  
  
  MM0 <-t(M120)%*%Minv[-1,-1]%*%M120
  beta2<-rep(tau,Q)                                        # prior point estimate = rep(tau,Q)                        
  Tvalue<-(1+t(beta2)%*%MM0%*%beta2)                              
  mse.point<-(Tvalue*Nruns/D)^(1./(P-1))                         # MSE(D)_s_point
  
  list (df = df, Ds=Ds, DP=DP, LoFDP=LoFDP, mseD=mseD, mse.point = mse.point, Ls=Ls, LP=LP, LoFLP=LoFLP, mseL=mseL)
}



### Obtain not orthogonal designs and simple criteria values: Ds, DP, L, LP

criteria.values<-function(X1.orth)     ## "orthogonal" matrix of primary terms, with labels
{
  cand<-candidate_set(Levels)
  labels<-as.vector(X1.orth[,1])       # labels
  index<-rep(0,Nruns)
  for (i in 1:Nruns){
    index[i]<-which(cand[,1]==labels[i])
  }
  X1<-cand[index,]                     # extended matrix, not orthonormalized, with labels
  
  M<-crossprod(X1[,-1])
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,6))               
  Ds<-(D/Nruns)^(1.0/(P-1))
  Minv<-solve(M)
  Ls<-1.0/(W%*%(diag(Minv)[-1]))        # tr(W(X'X)^-1) - weighted A  
  DF<-nlevels(as.factor(X1[,1]))        
  df<-Nruns-DF                          # df = "pure error" degrees of freedom
  if (df>0){
    DP<-Ds/qf(1-alpha.DP,P-1,df)
    LP<-Ls/qf(1-alpha.LP,1,df)
    LoF<-ifelse ((Nruns-P-df)>0,1.0/qf(1-alpha.LoF,Nruns-P-df,df),0)
  } else {DP<-0; LP<-0; LoF<-0;}
  compound<-(Ds^kappa.Ds)*(Ls^kappa.Ls)*(DP^kappa.DP)*(LP^kappa.LP)*(DF^kappa.DF)*(LoF^kappa.LoF)
  
  list (X1=X1, Ds=Ds, Ls=Ls, DP=DP, LP=LP, df=df, LoF=LoF, compound=compound)
}

##### GD, GDP GL, GLP criteria evaluation

criteria.values.G<-function(X1,X2,eps=10^-23)             # X1 - orthonormalised matrix
{
  DP<-0; LoFDP<-0; LP<-0; LoFLP<-0;
  DF<-nlevels(as.factor(X1[,1]))              # d.f. = N-number of unique design points
  df<-Nruns-DF                                # df - pure error degrees of freedom
  
  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))/Nruns
  if (D>eps) 
  {
    Minv<-solve(M)
  } else {return (eval=0)}
  if (D^(1.0/(P-1))>0)
  {
    Ds<-D^(-1.0/(P-1))                        # Ds component
  } else {return (eval=0)}
  Ls<-W%*%(diag(Minv)[-1])                    # Ls component
  
  M12<-crossprod(X1[,-1],X2[,-1])
  A<-Minv%*%M12
  L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)
  LoFD<-(det(L0))^(-1.0/Q)                                          # LoF (D)
  LoFL<-1./(sum(diag(L0))/Q)                                        # LoF (L)
    
  A0<-crossprod(A)+diag(1,nrow=Q)
  biasD<-prod(round(eigen(A0,only.values=TRUE)$values,10))^(1.0/Q)  # bias (D)
  biasL<-sum(diag(A0))/Q                                            # bias (L)
    
  if (df>0)
  {
    DP<-Ds*qf(1-alpha.DP,P-1,df)                               # DP component
    LP<-Ls*qf(1-alpha.LP,1,df)                                 # LP component
    LoFDP<-LoFD*qf(1-alpha.LoF,Q,df)                           # LoF (DP)   
    L0.inv.trace<-sum(1./eigen(L0,only.values=TRUE)$values)
    LoFLP<-L0.inv.trace*qf(1-alpha.LoFL,1,df)/Q                 # LoF (LP)
  } else {DP<-0;LoFDP<0; LoFLP<-0}
  
  list (Ds=Ds,DP=DP,LoFD=LoFD,LoFDP=LoFDP,biasD=biasD, Ls=Ls,LP=LP,LoFL=LoFL,LoFLP=LoFLP,biasL=biasL)
}

# Output function

output<-function(S.G)
{
  list1<-c(S.G$df, S.G$compound, S.G$time)
  list(out=list1, values=criteria.values.G(S.G$X1.orth, S.G$X2.orth))
}

# Obtaining pure design matrix, no labels, transformed to [-1,1]

extract.design<-function(X1)
{
  labels<-as.vector(X1[,1])                          # extract labels
  cand<-as.matrix(expand.grid(Levels))
  candl<-cbind(label(cand),cand)                     # labelling
  index<-rep(0,Nruns)
  for (i in 1:Nruns)
  {
    index[i]<-which(candl[,1]==labels[i])
  }
  design<-apply(cand[index,],2,transform.b)
}