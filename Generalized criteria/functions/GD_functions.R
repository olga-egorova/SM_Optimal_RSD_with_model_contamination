### Functions used to find GD-optimal designs (nearly optimal)

### Minimising Components: Ds, LoF(D), bias(D)

criteria.GD<-function(X1, X2, eps=10^-23)        # GD criteria evaluation; X1,X2 -   matrices of primary and potential terms, with labels
{
  Ds<-0; LoF<-0; bias<-0;
  DF<-nlevels(as.factor(X1[,1]))             
  df<-Nruns-DF                                # df - pure error degrees of freedom
  
  M<-crossprod(X1[,-1])
  D<-prod(round(eigen(M,symmetric=TRUE,only.values=TRUE)$values,8))/Nruns
  
  if (D>eps) 
  {
    Minv<-solve(M)
  } else {return (list (eval=0, Ds=0, DP=0, LoF=0, bias=0, df=df, compound=10^6));} # if M comp. singular
  
  if (kappa.Ds>0)
  {
    if (D^(1.0/(P-1))>0)
    {
      Ds<-D^(-1.0/(P-1))
    } else {return (list (eval=0, Ds=0, LoF=0, bias=0,compound=10^6));} # Ds  
  }
  
  if ((kappa.LoF>0) ||(kappa.bias>0))  # check for A calculation
  {
    M12<-crossprod(x=X1[,-1],y=X2[,-1])
    A<-Minv%*%M12
  }
  if (kappa.LoF>0)
  {
      L0<-crossprod(X2[,-1])-t(M12)%*%A+diag(1./tau2,nrow=Q)
      LoF<-(det(L0))^(-1.0/Q)                         
  }
  if (kappa.bias>0)
  {
    A0<-crossprod(A)+diag(1,nrow=Q)
    bias<-(det(A0))^(1.0/Q)
  }
  compound<-Ds^kappa.Ds*LoF^kappa.LoF*bias^kappa.bias
  list (eval=1,Ds=Ds, LoF=LoF, bias=bias, df=df, compound=compound)
}

### Swapping treatments 

swap.GD<-function(X1,X2,cand.full.orth) 
{
  Xcrit<-criteria.GD(X1,X2)
  Xcomp<-Xcrit$compound
  search<-0
  n<-nrow(cand.full.orth)
  for (l in 1:Nruns)
  {
    move<-ifelse(l==1||((l>1)&&(X1[l,1]!=X1[(l-1),1])),1,0)
    if (move == 1){
      Xc1<-X1 
      Xc2<-X2
      for (i in 1:n)
      {
        if (X1[l,1]!=cand.full.orth[i,1])         # look at labels
        {
          Xc1[l,]<-cand.full.orth[i,1:(P+1)]
          Xc2[l,]<-cand.full.orth[i,c(1,(P+2):(P+Q+1))]
          Ccrit<-criteria.GD(X1=Xc1,X2=Xc2)
          Ccomp<-Ccrit$compound
          if (Xcomp>Ccomp)                ### if the new design is better (minimising)
          {
            X1<-Xc1; X2<-Xc2
            Xcomp<-Ccomp
            search<-1
          }
        }
      }
    }
  }
  list (X1=X1, X2=X2, compound=Xcomp, search=search, crit=criteria.GD(X1,X2))
}

### Search function 

Search.GD<-function()
{
  start_time<-Sys.time()
  cand<-candidate_set(Levels)                       # form the candidate set of treatments, primary terms, not(!) orthonormal
  cand.full.orth<-candidate_set.orth(cand)          # candidate set, potential terms
  crit.values<-matrix(0,ncol=1, nrow=Nstarts)
  for (k in 1:Nstarts)
  {
    initial<-initial.orth(cand.full.orth)
    X1<-initial$X1
    X2<-initial$X2
    if (k==1) 
    {
      crit.opt<-criteria.GD(X1,X2)$compound
      X1.opt<-X1; X2.opt<-X2
    }
    s<-1
    while (s==1)
    {
      Xs<-swap.GD(X1,X2,cand.full.orth)
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
  criteria.opt<-criteria.GD(X1.opt,X2.opt)
  list (time=time, X1.orth=X1.opt, X2.orth=X2.opt, df=criteria.opt$df, Ds=criteria.opt$Ds,
        LoF=criteria.opt$LoF, bias=criteria.opt$bias, compound=criteria.opt$compound, path=crit.values)
}

