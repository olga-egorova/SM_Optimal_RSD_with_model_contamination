### Functions used to find unblocked MSE(LP)-optimal designs (nearly optimal)

### Minimizing Components: LP, LoF(LP), mse(L)

criteria.mseL<-function(X1,X2,eps=10^-23)      # X1, X2 -- matrices of primary and potential terms, both with labels
{                                             
  Ls<-0; LP<-0; LoF<-0; mse<-0;
  DF<-nlevels(as.factor(X1[,1]))              
  df<-Nruns-DF                                # df - pure error degrees of freedom
  
  M<-crossprod(X1[,-1])                       # information matrix of primary terms
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8)) 
  if (D>eps) 
  {
    Minv<-solve(M)
  } else {return (list (Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}
  
  if ((kappa.Ls>0)||(kappa.LP>0))
  {
   Ls<-W%*%(diag(Minv)[-1])                   # Ls
  } 
  if (kappa.LP>0)
  {
    if (df>0)
    {
      LP<-Ls*qf(1-alpha.LP,1,df)              # LPs
    } else {return (list (Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}
  }
  if (((kappa.LoF>0) && (df>0))||(kappa.mse>0))                        # check for A calculation
  {
    M12<-crossprod(X1[,-1],X2[,-1])
    A<-Minv%*%M12
  }
  if (kappa.LoF>0)
  {
    if (df>0)
    {
      L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)            # dispersion matrix + Iq/tau2
      L0.inv.trace<-Re(sum(1./eigen(L0,only.values=TRUE)$values))       # trace of the inverse matrix
      LoF<-L0.inv.trace*qf(1-alpha.LoFL,1,df)/Q
    } else {return (list (Ls=0, LP=0, LoF=0, bias=0, df=df, compound=10^6));}     
  }
  if (kappa.mse>0)
  {
    A0<-Minv+tau2*tcrossprod(A)
    mse<-sum(diag(A0)[-1])/(P-1)                                        # averaged trace of the A'A+Iq
  }
  compound<-LP^kappa.LP*LoF^kappa.LoF*mse^kappa.mse
  list (Ls=Ls, LP=LP, LoF=LoF, mse=mse, df=df, compound=compound)
}

### Swapping treatments 

swap.mseL<-function(X1,X2,cand.full) 
{
  Xcrit<-criteria.mseL(X1,X2)
  Xcomp<-Xcrit$compound
  search<-0
  n<-nrow(cand.full)
  for (l in 1:Nruns)
  {
    move<-ifelse(l==1||((l>1)&&(X1[l,1]!=X1[(l-1),1])),1,0)
    if (move == 1){
      Xc1<-X1 
      Xc2<-X2
      for (i in 1:n)
      {
        if (X1[l,1]!=cand.full[i,1])         # look at labels
        {
          Xc1[l,]<-cand.full[i,1:(P+1)]
          Xc2[l,]<-cand.full[i,c(1,(P+2):(P+Q+1))]
          Ccrit<-criteria.mseL(X1=Xc1,X2=Xc2)
          Ccomp<-Ccrit$compound
          if (Xcomp>Ccomp)                ### if the new design is better (minimizing)
          {
            X1<-Xc1; X2<-Xc2
            Xcomp<-Ccomp
            search<-1
          }
        }
      }
    }
  }
  list (X1=X1, X2=X2, compound=Xcomp, search=search, crit=criteria.mseL(X1,X2))
}

### Search function 

Search.mseL<-function()
{
  start_time<-Sys.time()
  cand<-candidate_set(Levels)                       # form the candidate set of treatments, primary terms
  cand.full<-candidate_set_full(cand)               # candidate set, potential terms
  crit.values<-matrix(0,ncol=1, nrow=Nstarts)
  for (k in 1:Nstarts)
  {
    initial<-initial.full(cand.full)
    X1<-initial$X1
    X2<-initial$X2
    if (k==1) 
    {
      crit.opt<-criteria.mseL(X1,X2)$compound
      X1.opt<-X1; X2.opt<-X2
    }
    s<-1
    while (s==1)
    {
      Xs<-swap.mseL(X1,X2,cand.full)
      X1<-Xs$X1; X2<-Xs$X2;
      s<-Xs$search
    }
    crit.values[k]<-Xs$compound       # track the change of criterion values
    if (crit.opt>Xs$compound)
    {
      X1.opt<-Xs$X1; X2.opt<-Xs$X2
      crit.opt<-Xs$compound
    }
  }
  finish_time<-Sys.time()
  time<-finish_time-start_time
  criteria.opt<-criteria.mseL(X1.opt,X2.opt)
  list (time=time, X1=X1.opt, X2=X2.opt, df=criteria.opt$df, Ls=criteria.opt$Ls, LP=criteria.opt$LP,
        LoF=criteria.opt$LoF, mse=criteria.opt$mse, compound=criteria.opt$compound, path=crit.values,
        kappa = c(kappa.LP, kappa.LoF, kappa.mse), tau2 = tau2)
}
