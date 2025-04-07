
##########################################################################################
####################### Setting up the blocked design search ############################
##########################################################################################


# Criteria: DPs + LoF(DP) + MSE(D)

rm(list=ls())    # Clear the environment
source("functions/Common_functions.R")
source("functions/MSEB_functions.R")
source("functions/Criteria_evaluation.R")



Nstarts <- 150     # Number of random starts of the search


# Three factors, five levels each (125 candidate points)

K<-3             # Number of factors
Nruns<-36        # Number of runs of the experiment
Nblocks<-2       # Number of blocks

Biter <- 500       # MC iterations

# Factors
Levels<-list(1:5,1:5,1:5)          # Levels of each factor
Cubic<-'Y'  
                               # 'Y' - cubic, 'N' - spheric coordinates
# Equal blocks
blocksize<-18
BlockSize<-matrix(rep(blocksize,Nblocks),ncol=1)
b<-as.factor(rep(1:Nblocks,BlockSize))
B<-as.matrix(model.matrix(~-1+b))
Z<-diag(rep(1,Nruns))-B%*%solve(crossprod(B))%*%t(B)
y<-rnorm(Nruns)   # needed further to calculate df 

# Model, P and Q
Parameters<-c(rep(1,K),rep(1,K), rep(1,K*(K-1)/2))
P<-sum(Parameters)                                   # number of all parameters, no intercept (as we have blocks)
potential.terms<-1
if (potential.terms==1) {Q<-10;}                     # number of potential terms
alpha.DP<-0.05; alpha.LoF<-0.05; 

# Weights to be set to zero
kappa.Ds<-0; kappa.Ls<-0; kappa.LP<-0; kappa.DF<-0;



##########################################################################################
####################### Criteria parameters, search procedure ############################
##########################################################################################


#### scaling parameter of the prior variance, squared; tau2 = 1 or tau2 = 1./Q in the manuscript. 
tau2<-1.0                                 
tau<-tau2^0.5


#### Matrix of weights: each row contains a combination of weights. 
#### MSE(D) weights (kappa.DP, kappa.LoF, kappa.mse), sum of the weights = 1

m.kappa <- matrix(c(rep(1./3, 3),
                    0.4, 0.2, 0.4,
                    0.25, 0.25, 0.5,
                    1,0,0,           ## Last three rows are used to 
                    0,1,0,           ## calculate efficiency values
                    0,0,1), ncol = 3, byrow = TRUE)

## optional -- initiating pseudo-random number generation for each weight combination
seeds = c(465, 9682, 2226, 8890, 659, 8462)


## output: weights, P.E. df, criteria values
crit.names = c("df", "lof", "DP", "LoFDP", "mseD", "compound")
df.output <- data.frame(matrix(NA, ncol = length(crit.names) + ncol(m.kappa), 
                               nrow = nrow(m.kappa)))
colnames(df.output) <- c("kappa.DP", "kappa.LoFDP", "kappa.mseD", crit.names)

S.mseB <- vector(mode = "list", length = nrow(m.kappa))


##########################################################################################
####################### Criteria parameters, search procedure ############################
##########################################################################################


for (k in 1:nrow(m.kappa)){   # might choose a subset of weights
  
  if (sum(m.kappa[k,]) != 1.0) {
    print(paste("Row" ,k,": sum of the weights should be equal to 1."))
    next
  }
  
  kappa.DP <- m.kappa[k,1]; kappa.LoF <- m.kappa[k,2]; kappa.mse <- m.kappa[k,3]
  
  set.seed(seeds[k])
  S.current <- Search.GD_PB.b()
  print(paste("Last iteration time: ", S.current$time))
  out <- output(S.current)
  S.mseB[[k]] <- S.current
  
  df.output[k,] <- c(m.kappa[k,], S.current$df, 
                     Nruns - P - Nblocks - S.current$df,
                     out$values$DP, out$values$LoFDP, out$values$mseD_point,
                     S.current$compound)

  save.image("output/workspaces/MSEB_search.RData")
  write.csv(df.output, "output/tables/MSEB_values.csv")
  
}


############################################
############# Saving the designs ###########
############################################

design_list = vector(mode = "list", length = nrow(m.kappa)-3)

d <- 1
for (k in 1:length(S.mseB)) {
  S = S.mseB[[k]]
  if(!is.null(S)) {
    design_list[[d]] <- data.frame(S$X1[,2:(K+1)])
    d<- d+1
  }
}

### Writing all the designs in one file
writexl::write_xlsx(design_list, "output/designs/MSEDB_designs.xlsx")

