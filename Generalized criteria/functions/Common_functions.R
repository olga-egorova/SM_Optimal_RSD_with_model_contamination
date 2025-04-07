### Functions common for all Generalised criteria (Example 1) 

### Forming the candidate set as all possible combinations of factors levels
candidate_set<-function(Levels)
{
  cand<-as.matrix(expand.grid(Levels))
  candl<-cbind(label(cand),cand)                           # labelling
  candlt<-cbind(candl[,1],apply(candl[,-1],2,transform))   # rescaling
  if (Cubic=='N') {candlt<-as.matrix(spheric(candlt))}     # spheric coords
  candset<-extend(candlt)                                  #extended matrix
  return (candset)
}

label<-function(treat_set)                  # treat_set is not scaled to [-1,1]
{
  if (K>1){
    n<-c(rep(1,K-1),length(Levels[[K]]))        
    levelprod<-rep(1,K)
    for (i in (K-1):1)                    
    {
      n[i]<-length(Levels[[i]])             # n[i] - number of levels of the i-th factor
      levelprod[i]<-levelprod[i+1]*n[i+1]   # levelprod[i]=n[i+1]*...*n[K]
    }
    
    N<-nrow(treat_set)
    Label<-rep(1,N)
    m<-as.numeric(lapply(Levels, min))
    for (i in 1:N)                          # labelling
    {
      for (j in 1:K)
      {
        Label[i]<-Label[i]+(treat_set[i,j]-m[j])*levelprod[j]
      }
    }
  } else {Label<-seq(1:length(Levels[[1]]))}
  return(matrix(Label,ncol=1))
}

transform<-function(x)            # transformation to [-1,1]   
{
  a<-min(x)
  b<-max(x)
  xc<-(2*x-(a+b))/(b-a)           #to [-1,1]
  return (xc)
}

spheric<-function(candlt)          # cand - treatment matrix, with labels
{
  treat<-candlt[,-1]
  cand<-treat[apply(treat^2,1,sum)==0,]
  for (i in 1:K)
  {
    cand<-rbind(cand, sqrt(K/i)*treat[apply(treat^2,1,sum)==i,])
  }
  return (cbind(candlt[,1],cand))
}

### Extended design matrix, with labels

extend<-function(candlt)            # treatment matrix (with labels) -> extended design matrix (with labels)
{
  treat<-as.matrix(candlt[,-1])
  X<-matrix(1,nrow=nrow(candlt),ncol=1)  # intercept
  X<-cbind(X,treat[,Parameters[1:K]==1])
  count<-ncol(X)-1
  for (i in 1:K)
  {
    count<-count+1
    if (Parameters[count]>0) X<-cbind(X,treat[,i]^2)
  }
  if (K>=2)
  {
    for (i in 1:(K-1))
    {
      for (j in (i+1):K)
      {
        count<-count+1
        if (Parameters[count]>0) X<-cbind(X,treat[,i]*treat[,j])
      }
    }
  }
  X<-cbind(candlt[,1],X)
  return (X)
}


### Generating initial design

initial<-function(cand,Nruns=Nruns)          # cand - matrix with labels
{
  eps<-10^(-6)
  det<-0
  while (det<eps)
  {
    index<-sample(1:nrow(cand),size=Nruns,replace=TRUE)
    X<-cand[index,]
    det<-round(prod(eigen(t(X[,-1])%*%X[,-1],symmetric=TRUE,only.values=TRUE)$values),6)
  }
  list(X=X, det=det)
}

# Obtaining matrix of potential terms only, not orthonormalised

potential.matrix<-function(cand)              # 5 factors, 3 levels each. primary - quadratic model.
{
  cand.X2<-matrix(0,nrow=nrow(cand),ncol=Q)
  if (potential.terms==1)                     # all L*Q terms
  {
    cand.X2[,1]<-cand[,3]*cand[,9]; cand.X2[,2]<-cand[,3]*cand[,10]; cand.X2[,3]<-cand[,3]*cand[,11]; cand.X2[,4]<-cand[,3]*cand[,12];
    cand.X2[,5]<-cand[,4]*cand[,8]; cand.X2[,6]<-cand[,4]*cand[,10]; cand.X2[,7]<-cand[,4]*cand[,11]; cand.X2[,8]<-cand[,4]*cand[,12];
    cand.X2[,9]<-cand[,5]*cand[,8]; cand.X2[,10]<-cand[,5]*cand[,9]; cand.X2[,11]<-cand[,5]*cand[,11]; cand.X2[,12]<-cand[,5]*cand[,12];
    cand.X2[,13]<-cand[,6]*cand[,8]; cand.X2[,14]<-cand[,6]*cand[,9]; cand.X2[,15]<-cand[,6]*cand[,10]; cand.X2[,16]<-cand[,6]*cand[,12];
    cand.X2[,17]<-cand[,7]*cand[,8]; cand.X2[,18]<-cand[,7]*cand[,9]; cand.X2[,19]<-cand[,7]*cand[,10]; cand.X2[,20]<-cand[,7]*cand[,11];
  }
  if (potential.terms==2)             # all L*L*L terms
  {
    cand.X2[,1]<-cand[,3]*cand[,17]; cand.X2[,2]<-cand[,3]*cand[,18]; cand.X2[,3]<-cand[,3]*cand[,19];
    cand.X2[,4]<-cand[,3]*cand[,20]; cand.X2[,5]<-cand[,3]*cand[,21]; cand.X2[,6]<-cand[,3]*cand[,22];
    cand.X2[,7]<-cand[,4]*cand[,20]; cand.X2[,9]<-cand[,4]*cand[,21]; cand.X2[,9]<-cand[,4]*cand[,22];
    cand.X2[,10]<-cand[,5]*cand[,22];  
  }
  if (potential.terms==3)
  {
    cand.X2[,1]<-cand[,3]*cand[,9]; cand.X2[,2]<-cand[,3]*cand[,10]; cand.X2[,3]<-cand[,3]*cand[,11]; cand.X2[,4]<-cand[,3]*cand[,12];   # Q*L
    cand.X2[,5]<-cand[,4]*cand[,8]; cand.X2[,6]<-cand[,4]*cand[,10]; cand.X2[,7]<-cand[,4]*cand[,11]; cand.X2[,8]<-cand[,4]*cand[,12];
    cand.X2[,9]<-cand[,5]*cand[,8]; cand.X2[,10]<-cand[,5]*cand[,9]; cand.X2[,11]<-cand[,5]*cand[,11]; cand.X2[,12]<-cand[,5]*cand[,12];
    cand.X2[,13]<-cand[,6]*cand[,8]; cand.X2[,14]<-cand[,6]*cand[,9]; cand.X2[,15]<-cand[,6]*cand[,10]; cand.X2[,16]<-cand[,6]*cand[,12];
    cand.X2[,17]<-cand[,7]*cand[,8]; cand.X2[,18]<-cand[,7]*cand[,9]; cand.X2[,19]<-cand[,7]*cand[,10]; cand.X2[,20]<-cand[,7]*cand[,11];
    
    cand.X2[,21]<-cand[,3]*cand[,17]; cand.X2[,22]<-cand[,3]*cand[,18]; cand.X2[,23]<-cand[,3]*cand[,19];                               # L*L*L
    cand.X2[,24]<-cand[,3]*cand[,20]; cand.X2[,25]<-cand[,3]*cand[,21]; cand.X2[,26]<-cand[,3]*cand[,22];
    cand.X2[,27]<-cand[,4]*cand[,20]; cand.X2[,28]<-cand[,4]*cand[,21]; cand.X2[,29]<-cand[,4]*cand[,22];
    cand.X2[,30]<-cand[,5]*cand[,22];
  }
  X2<-cbind(cand[,1],cand.X2)        # adding column with labels
  return (X2)
}

### Orthonormalisation

candidate_set.orth<-function(cand)     # returns candidate set for potential terms, with labels, orthogonalised
{                                      # cand - full extended candidate set, primary terms (full quadratic model)
  potential<-potential.matrix(cand)                                    # potential -- orthogonalised matrix for extended model, labelled
  cand.full<-cbind(cand[,-1],potential[,-1])                           # extended model matrix, no labels
  cand.full.orth<-cbind(cand[,1],orthonormalization(cand.full,basis=FALSE))        # orthonormalisation, adding labels
  return (cand.full.orth)
}

initial.orth<-function(cand.full.orth)          # cand.full.orth -- extended model matrix, orthonormalised, with labels
{                                               # returns initial primary and potential matrices
  eps<-10^(-23)
  det<-0
  while (det<eps)
  {
    index<-sample(1:nrow(cand.full.orth),size=Nruns,replace=TRUE)
    X1<-cand.full.orth[index,1:(P+1)]
    X2<-cand.full.orth[index,c(1,(P+2):(P+Q+1))]
    det<-prod(round(eigen(t(X1[,-1])%*%X1[,-1],symmetric=TRUE,only.values=TRUE)$values,6))   # det of the information matrix
  }
  list(X1=X1, X2=X2, det=det)               # X1, X2 -- matrices of primary and potential terms, both with labels
}
