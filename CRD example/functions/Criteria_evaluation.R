### Obtain individual criteria values: Ds, DP, L, LP

criteria.values<-function(X1)     
{
  M<-crossprod(X1[,-1])
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,6))               
  Ds<-(D/Nruns)^(1.0/(P-1))
  Minv<-solve(M)
  Ls<-1.0/(W%*%(diag(Minv)[-1]))        # tr(W(X'X)^-1) - weighted A  
  DF<-nlevels(as.factor(X1[,1]))        # df = pure error degrees of freedom
  df<-Nruns-DF
  if (df>0){
    DP<-Ds/qf(1-alpha.DP,P-1,df)
    LP<-Ls/qf(1-alpha.LP,1,df)
    LoF<-ifelse ((Nruns-P-df)>0,1.0/qf(1-alpha.LoF,Nruns-P-df,df),0)
  } else {DP<-0; LP<-0; LoF<-0;}
  compound<-(Ds^kappa.Ds)*(Ls^kappa.Ls)*(DP^kappa.DP)*(LP^kappa.LP)*(DF^kappa.DF)*(LoF^kappa.LoF)
  
  list (Ds=Ds, Ls=Ls, DP=DP, LP=LP, df=df, LoF=LoF, compound=compound)
}


### Evaluation of MSE-based criteria functions

criteria.values.mse<-function(X1,X2,eps=10^-23, Biter=100)      # X1,X2 - primary and potential matrices, with treatment labels
{
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
  } else {DP<-0;LoFDP<0; LoFLP<-0}
  
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
    
  list (Ds=Ds, DP=DP, LoFDP=LoFDP, mseD=mseD, Ls=Ls, LP=LP, LoFLP=LoFLP, mseL=mseL)
}

### Obtaining a point estimate of the MSE(Ds)-component

MSE.point<-function(S.mse)
{
  M12<-t(S.mse$X1[,-c(1,2)])%*%Z0%*%(S.mse$X2[,-1])
  M<-t(S.mse$X1[,-1])%*%S.mse$X1[,-1]
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))
  Minv<-solve(M)
  MM<-t(M12)%*%Minv[-1,-1]%*%M12
  beta2<-rep(tau,Q)                                        # prior point estimate = rep(tau,Q)                        
  Tvalue<-(1+t(beta2)%*%MM%*%beta2)                              
  mse.point<-(Tvalue*Nruns/D)^(1./(P-1))                         # MSE(D)_s_point
  return (mse.point)
}

### Output function

output<-function(S.mse)
{
  list1<-c(S.mse$df, S.mse$compound, S.mse$mse, S.mse$time)
  list(out=list1, 
       values=criteria.values.mse(S.mse$X1, S.mse$X2), 
       mse.point=MSE.point(S.mse))
}





