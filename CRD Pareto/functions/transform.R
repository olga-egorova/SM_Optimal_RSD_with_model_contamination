# transformation to [-1,1]
Transform<-function(x)           
{
  a<-min(x)
  b<-max(x)
  xc<-(2*x-(a+b))/(b-a)
  return (xc)
}
