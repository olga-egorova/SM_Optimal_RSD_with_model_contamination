
### Forming the candidate set as all possible combinations of factors levels
candidate_set.b<-function(Levels)
{
  cand<-as.matrix(expand.grid(Levels))
  candl<-cbind(label.b(cand),cand)                           # labelling
  candlt<-cbind(candl[,1],apply(candl[,-1],2,transform.b))   # rescaling
  if (Cubic=='N') {candlt<-as.matrix(spheric.b(candlt))}     # spheric coords
  candset<-extend.b(candlt)                                  #extended matrix
  return (candset)
}

candidate_set_full.b<-function(cand)          # full candidate matrix, with labels  
{
  cand.full<-cbind(cand,potential.matrix.b(cand)[,-1])
  return (cand.full)
}


label.b<-function(treat_set)                  # treat_set is not scaled to [-1,1]
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

transform.b<-function(x)            # transformation to [-1,1]   
{
  a<-min(x)
  b<-max(x)
  xc<-(2*x-(a+b))/(b-a)           #to [-1,1]
  return (xc)
}

spheric.b<-function(candlt)          # cand - treatment matrix, with labels
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

extend.b<-function(candlt)            # treatment matrix (with labels) -> extended design matrix (with labels)
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
  X<-X[,-1]
  X<-cbind(candlt[,1],X)              # removing intercept 
  return (X)
}


### Generating initial design, ensuring that it contains two center points per block

initial.full.b<-function(cand.full)          # cand.full - matrix with labels
{
  eps<-10^(-6)
  det<-0
  while (det<eps)
  {
    index<-c(sample(1:nrow(cand.full),size=Nruns/Nblocks-CP,replace=TRUE),rep(CP_ind,CP),sample(1:nrow(cand.full),size=Nruns/Nblocks-CP,replace=TRUE),rep(CP_ind,CP))
    X1<-cand.full[index,1:(P+1)]
    X2<-cand.full[index,c(1,(P+2):(P+Q+1))]
    det<-round(prod(eigen(t(X1[,-1])%*%X1[,-1],symmetric=TRUE,only.values=TRUE)$values),6)
  }
  list(X1=X1, X2=X2, det=det)
}


# Obtaining matrix of potential terms only

potential.matrix.b<-function(cand.b)    # 3 factors, 5 levels each. primary - quadratic model. potential -- all up to ^3  
{
  cand.X2<-matrix(0,nrow=nrow(cand.b),ncol=Q)
  cand<-cbind(cand.b[,1],matrix(rep(1,nrow(cand.b)),ncol=1),cand.b[,-1])
  if (potential.terms==1)                     # all C,L*Q and L*L*L terms
  {
    cand.X2[,1]<-cand[,3]*cand[,6]; cand.X2[,2]<-cand[,4]*cand[,7]; cand.X2[,3]<-cand[,5]*cand[,8];    # cubic terms
    cand.X2[,4]<-cand[,6]*cand[,4]; cand.X2[,5]<-cand[,6]*cand[,5];                                   # Q*L terms
    cand.X2[,6]<-cand[,7]*cand[,3]; cand.X2[,7]<-cand[,7]*cand[,5]; 
    cand.X2[,8]<-cand[,8]*cand[,3]; cand.X2[,9]<-cand[,8]*cand[,4];
    cand.X2[,10]<-cand[,9]*cand[,5];                                                                  # L*L*L term
  }
  X2<-cbind(cand.b[,1],cand.X2)        # adding column with labels
  return (X2)
}
