############### Generating initial random design ##########################

initial.full = function(cand.full) {
  eps = 10^(-6); det = 0; df = 0;
  
  while ((det<eps) || (df == 0)) # making sure the initial design contains replicates
  {
    index = sample(1:nrow(cand.full), size=Nruns, replace=TRUE)
    Xp =  as.matrix(cand.full[index, 2:(P+1)])  
    det = round(prod(eigen(crossprod(Xp),symmetric=TRUE,only.values=TRUE)$values),6)
    Labels = c(cand.full[index, 1])
    df = nrow(Xp) - nlevels(as.factor(Labels))
  }
  Xq = as.matrix(cand.full[index, (P+2):(P+Q+1)])
  
  list(Labels = Labels, Xp=Xp, Xq = Xq, det=det)
}

################ Swapping treatments #######################################

## cand.full -- matrix of [labels, cand primary terms, cand potential terms]; 1+P+Q columns

pareto.swap.full.L = function(Labels, Xp, Xq, cand.full, list.nd)  
{
  P = ncol(Xp); Q = ncol(Xq)
  Xpareto = pareto_criteria_trace(Labels, Xp, Xq)
  
  search = 0
  n = nrow(cand.full)
  for (l in 1:Nruns)  
  {
    move = ifelse(l == 1||((l>1)&&(Labels[l] != Labels[(l-1)])),1,0)
    if (move == 1){
      Lc = Labels; Xcp = Xp; Xcq = Xq;
      for (i in 1:n)
      {
        if (Labels[l] != cand.full[i,1])         # look at labels
        {
          Xcp[l,] = cand.full[i, 2:(P+1)]
          Xcq[l,] = cand.full[i,(P+2):(P+Q+1)]
          Lc[l] = cand.full[i,1]
          Cpareto = pareto_criteria_trace(Labels = Lc, Xp = Xcp, Xq = Xcq)
          
          if (!is.dominated(p.crit1 = Cpareto, p.crit2 = Xpareto)){ # if the new one is not dominated 
            ## Compare with the list of non-dominated designs (should do here or only if new is accepted?)
            upd = update.nd_list(list.nd, Labels = Lc, p.crit = Cpareto)
            list.nd = upd$l
            if (upd$added) {search = 1}  # if the new design has been added to the list
            
            if (is.dominated(p.crit1 = Xpareto, p.crit2 = Cpareto))   ### if the new design is better (minimising) -- new "threshold"
            {
              Xp = Xcp; Xq = Xcq; Labels = Lc; 
              Xpareto = Cpareto;
              #search = 1
            }
          } 
        }
      }
    }
  }
  list (Labels = Labels, Xp=Xp, Xq=Xq, search=search, p.crit = Xpareto, list.nd = list.nd)
} 

############# Pareto search function ##############################

Search.pareto.L = function(cand.full)
{
  
  cand.full = as.matrix(cand.full)
  list.nd = list()
  
  start_time = Sys.time()
  for (k in 1:Nstarts)
  {
    initial = initial.full(cand.full)
    Labels = initial$Labels
    Xp = initial$Xp; Xq = initial$Xq
    if (k==1) 
    {
      pareto.opt = pareto_criteria_trace(Labels, Xp, Xq)
      # initiate the set of non-dominant designs and criteria values
      list.nd = update.nd_list(list.nd, Labels = Labels, p.crit = pareto.opt)$l
    }
    s = 1
    while (s==1)
    {
      Xs = pareto.swap.full.L(Labels, Xp, Xq, cand.full, list.nd)
      Labels = Xs$Labels; Xp = Xs$Xp; Xq = Xs$Xq;
      list.nd = Xs$list.nd
      s = Xs$search
    }
  }
  finish_time = Sys.time()
  time = finish_time - start_time
  
  list.designs = list()
  for (l in 1:length(list.nd)){
    # make cand. full a data table
    list.designs = append(list.designs, list(append(list.nd[[l]], 
                                                    list(X.full = cand.full[list.nd[[l]]$labels,]))))
  }
  
  return (list (time=time, designs = list.designs, list.nd = list.nd))
}
