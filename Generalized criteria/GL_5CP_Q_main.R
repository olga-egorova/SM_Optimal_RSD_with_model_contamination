##########################################################################################
## Main script -- obtaining GL-optimal designs, adding 5 center points and 
## evaluating DP-, LoF(DP)-, MSE(D)- and LP-, LoF(LP)- and MSE(L)-components
##########################################################################################

source("functions/GL_functions.R")
source("functions/Common_functions.R")
source("functions/Criteria_evaluation.R")
library(far)


Nstarts<-100     # Number of random starts of the search



##### Set up with 5 added center points -- hence, initial 35 runs.

# Five factors, three levels each (243 candidate points). 
# GL-optimal designs (Ls+LoF(L)+bias(L))
Levels<-list(1:3,1:3,1:3,1:3,1:3)          # Levels of each factor
Cubic<-'Y'                                 # 'Y' - cubic, 'N' - spheric coordinates

K <- 5           # Number of factors
Ncenter <- 5     # Number of center points
center.label <- (3^K+1)/2
Nruns <- 40 - Ncenter        # Number of runs of the experiment

# Model, P and Q
Parameters<-c(rep(1,K),rep(1,K), rep(1,K*(K-1)/2))
P<-sum(Parameters)+1                       # number of all parameters, including intercept
potential.terms<-3
if (potential.terms==3) {Q<-30;}

tau2<-1.0/Q         # scaling parameter of the prior variance, squared
tau<-sqrt(tau2)

# Weights for parameters (Ls-optimality)
W<-matrix(1,nrow=1,ncol=P-1)
for (i in (K+1):(2*K))
{
  if (Parameters[i]==1)
  {W[sum(Parameters[1:i])]<-W[sum(Parameters[1:i])]/4}
}
W<-matrix(W/sum(W),nrow=1)

# Confidence levels

MC<-'Y'          # 'Y' - Yes, 'N' - No
prob.LP<-0.95; prob.LoFL<-0.95;
prob.LP<-ifelse(MC=='Y',prob.LP^(1./(P-1)),prob.LP)            # Correction for multiple comparisons             
prob.LoFL<-ifelse(MC=='Y',prob.LoFL^(1./Q),prob.LoFL)

alpha.DP<-0.05; alpha.LP<-1-prob.LP; alpha.LoF<-0.05; alpha.LoFL<-1-prob.LoFL; alpha.bias<-0.05;
# Weights put to zero
kappa.DP<-0.0; kappa.Ds<-0; kappa.LP<-0

#### Matrix of weights: each row contains a combination of weights. 
#### GL weights (kappa.Ls, kappa.LoF(L), kappa.bias(L)), sum of the weights = 1

m.kappa <- matrix(c(1,0,0,     ## First three rows are used to calculating efficiency values
                    0,1,0,
                    0,0,1,
                    0.5, 0.5, 0,
                    0.5, 0, 0.5,
                    0, 0.5, 0.5,
                    rep(1./3, 3),
                    0.5, 0.25, 0.25,
                    0.25, 0.5, 0.25,
                    0.25, 0.25, 0.5),
                  ncol = 3, byrow = TRUE)


## optional -- initiating pseudo-random number generation for each weight combination
seeds = c(16633, 141621,169983, 38278,  53079, 64567, 184522, 53912, 148296,  32416)    # for tau2 = 1./Q

S.GL <- vector(mode = "list", length = nrow(m.kappa))

## output: weights, P.E. df, criteria values
crit.names = c("criterion", "CPs", "tau2", "df", "LP", "LoFLP", "mseL", "compound")

## MSE(LP) optimal values
MSE_optim.Q <- read.csv("output/MSELP_Q_values.csv", header = TRUE)[, c(7:9, 14)]
df.output <- data.frame(matrix(NA, ncol = length(crit.names) + ncol(m.kappa), 
                               nrow = nrow(m.kappa)))
colnames(df.output) <- c("kappa.1", "kappa.2", "kappa.3", crit.names)



############# Search and design evaluation #################



for (k in 1:nrow(m.kappa)){   # 
  
  if (sum(m.kappa[k,]) != 1.0) {
    print(paste("Row" ,k,": sum of the weights should be equal to 1."))
    next
  }
  
  kappa.Ls <- m.kappa[k,1]; kappa.LoF <- m.kappa[k,2]; kappa.bias <- m.kappa[k,3]
  
  set.seed(seeds[k])
  S.current <- Search.GL()
  
  S.GL[[k]] <- S.current
  print(paste("Last iteration time: ", S.current$time))
  
  current.labels <- c(S.current$X1.orth[,1], 
                      rep(center.label, Ncenter))
  m.cand <- candidate_set(Levels)
  
  ## Add center points, evaluate the design
  
  current.labels <- c(S.current$X1.orth[,1], 
                      rep(center.label, Ncenter))
  m.cand <- candidate_set_full(candidate_set(Levels))
  m.X <- matrix(nrow=0, ncol = ncol(m.cand))
  
  for (ind in 1:length(current.labels)){
    m.X <- rbind(m.X, m.cand[m.cand[,1] %in% current.labels[ind],])
  }
  
  X1<-m.X[,1:(P+1)]
  X2<-m.X[,c(1,(P+2):(P+Q+1))]
  
  out <- criteria.values.mse(X1 = X1, X2 = X2) 
  compound <- prod(c(out$LP, out$LoFLP, out$mseL)^m.kappa[k,1:3])
  
  df.output[k,] <- c(m.kappa[k,], "GL", "5CP", "1/q", 
                     out$df, 
                     MSE_optim.Q[k,]/c(out$LP, out$LoFLP, out$mseL, compound))
  
  
  save.image("output/workspaces/GL_5CP_Q.RData")
  write.csv(df.output, "output/tables/GL_5CP_Q.csv")
  
}


