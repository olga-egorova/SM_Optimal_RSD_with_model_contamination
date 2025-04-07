#Create a label column 
label<-function(treat_set, Levels, K)                  # treat_set is not scaled to [-1,1]
{
  if (K>1){
    n <- c(rep(1,K-1),length(Levels[[K]]))
    levelprod <- rep(1,K)
    for (i in (K-1):1)
    {
      n[i] <- length(Levels[[i]])             # n[i] - number of levels of the i-th factor
      levelprod[i] <- levelprod[i+1]*n[i+1]   # levelprod[i]=n[i+1]*...*n[K]
    }

    N <- nrow(treat_set)
    Label <- rep(1, N)
    m <- as.numeric(lapply(Levels, min))
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
