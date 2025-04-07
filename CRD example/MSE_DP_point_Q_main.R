### Completely randomized design search.
### Five factors, three levels each (243 candidate points). 

# MSE(DPs)-optimal designs (DPs+LoF(DP)+mse(D)), with point prior

rm(list=ls())    # Clear the environment

K<-5             # Number of factors
Nruns<-40        # Number of runs of the experiment
Nstarts<-100     # Number of random starts of the search
Biter <- 500

Z0<-diag(1,Nruns)-matrix(1/Nruns,nrow=Nruns,ncol=Nruns)  

# Factors
Levels<-list(1:3,1:3,1:3,1:3,1:3)          # Levels of each factor
Cubic<-'Y'                                 # 'Y' - cubic, 'N' - spheric coordinates

# Model, P and Q
Parameters<-c(rep(1,K),rep(1,K), rep(1,K*(K-1)/2))
P<-sum(Parameters)+1                       # number of all parameters, including intercept
potential.terms<-3
if (potential.terms==1) {Q<-20;}           # number of potential terms
if (potential.terms==2) {Q<-10;}
if (potential.terms==3) {Q<-30;}

# Weights for parameters (for LPs-optimality evaluation)
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

alpha.DP<-0.05; alpha.LP<-1-prob.LP; alpha.LoF<-0.05; alpha.LoFL<-1-prob.LoFL; 
kappa.Ds<-0;

##########################################################################################
####################### Criteria parameters, search procedure ############################
##########################################################################################

source("functions/Common_functions.R")
source("functions/MSE.D_functions.R")
source("functions/MSE.D_point_functions.R")
source("functions/Criteria_evaluation.R")

#### scaling parameter of the prior variance, squared; tau2 = 1 or tau2 = 1./Q in the manuscript. 
tau2<-1.0/Q                                 
tau<-tau2^0.5

#### Matrix of weights: each row contains a combination of weights. 
#### MSE(D) weights (kappa.DP, kappa.LoF, kappa.mse), sum of the weights = 1

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
seeds = c(74361, 95159, 118292, 92904, 60251, 96266, 121351, 150162, 78249, 137645)    # for tau2 = 1./Q

S.mseD.point <- vector(mode = "list", length = nrow(m.kappa))

## output: weights, P.E. df, criteria values
crit.names <- c("df", "lof", 
               "DP", "LoFDP", "mseD.point", "mseD", "compound", 
               "LP", "LoFLP", "mseL")

df.output <- data.frame(matrix(NA, ncol = length(crit.names) + ncol(m.kappa), 
                               nrow = nrow(m.kappa)))
colnames(df.output) <- c("kappa.DP", "kappa.LoFDP", "kappa.mseD", crit.names)



for (k in 1:nrow(m.kappa)){   
  
  if (sum(m.kappa[k,]) != 1.0) {
    print(paste("Row" ,k,": sum of the weights should be equal to 1."))
    next
  }
  
  kappa.DP <- m.kappa[k,1]; kappa.LoF <- m.kappa[k,2]; kappa.mse <- m.kappa[k,3]
  
  set.seed(seeds[k])
  S.current <- Search.mseD.point()
  print(paste("Last iteration time: ", S.current$time))
  out = output(S.current)
  
  S.mseD.point[[k]] <- S.current
  
  df.output[k,] <- c(m.kappa[k,],
                     S.current$df, Nruns - P - S.current$df,
                     out$values$DP, out$values$LoFDP, out$mse.point,
                     out$values$mseD, S.current$compound,
                     out$values$LP, out$values$LoFLP, out$values$mseL)
  
  save.image("output/workspaces/MSEDP_point_Q.RData")
  write.csv(df.output, "output/tables/MSEDP_point_Q_values.csv")
  
}

############################################
############# Saving the designs ###########
############################################

design_list = vector(mode = "list", length = nrow(m.kappa)-3)

d <- 1
for (k in 1:length(S.mseD.point)) {
  S = S.mseD.point[[k]]
  if(!is.null(S)) {
    design_list[[d]] <- data.frame(S$X1[,3:(K+2)])
    d <- d+1
  }
}

### Writing all the designs in one file
writexl::write_xlsx(design_list, "output/designs/MSEDP_point_Q_designs.xlsx")




