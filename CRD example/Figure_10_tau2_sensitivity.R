
################################################################################################
### Generating a figure with relative efficiency values of MSE(DP)-optimal (compound) designs,
### for different values of tau2
################################################################################################

library(ggtern)
library(tidyr)
library(dplyr)
rm(list=ls()) 

### Load the table containing the criteria values of the designs, for both values of tau2

MSEDP.output <- read.csv("output/tables/MSEDP_point_values.csv", header = TRUE)
MSEDP_Q.output <- read.csv("output/tables/MSEDP_point_Q_values.csv", header = TRUE)

### Optimal values of the component criteria
crit.optim <- list("DP" = MSEDP.output[1, "DP"],
                   "LoFDP" = MSEDP.output[2, "LoFDP"],
                   "MSE_D" = MSEDP.output[3, "mseD.point"])
crit.optim.Q <- list("DP" = MSEDP_Q.output[1, "DP"],
                     "LoFDP" = MSEDP_Q.output[2, "LoFDP"],
                     "MSE_D" = MSEDP_Q.output[3, "mseD.point"])


crit.optim <- as.numeric(Re(unlist(crit.optim)))
crit.optim.Q <- as.numeric(Re(unlist(crit.optim.Q)))

m.kappa <- as.matrix(MSEDP.output[, 2:4])
relative_names <- c("DP", "LoFDP", "mseD", "compound")
df_relative <- data.frame(matrix(NA, ncol = length(relative_names) + ncol(m.kappa), 
                                 nrow = nrow(m.kappa)))
df_relative_Q <- data.frame(matrix(NA, ncol = length(relative_names) + ncol(m.kappa), 
                                 nrow = nrow(m.kappa)))

colnames(df_relative) <- colnames(df_relative_Q) <-
  c("kappa.DP", "kappa.LoFDP", "kappa.mseD", relative_names)


###########################################################################################
################# Efficiency values of (tau2 = 1)-optimal designs w.r.t. (tau2 = 1/q) criteria
###########################################################################################

load("output/workspaces/MSEDP_point.RData")

## set tau2 to 1/q -- to obtain the relative efficiency values
tau2 <- 1./Q
tau <- sqrt(tau2)

for (k in 1:nrow(m.kappa)) {
  S = S.mseD.point[[k]]
  
  if(!is.null(S)) {
    
    kappa.DP <- m.kappa[k,1]; kappa.LoF <- m.kappa[k,2]; kappa.mse <- m.kappa[k,3]
    if (!all(m.kappa[k,] == S$kappa)) {print("Kappa-s"); next}

    eval <- criteria.mseD.point(S$X1, S$X2)  # design evaluation
    out <- output(S)
    eff <- as.numeric(Re(c(crit.optim.Q*100/c(out$values$DP, out$values$LoFDP, out$mse.point),
                           MSEDP_Q.output$compound[k]*100/eval$compound)))
    eff[which(eff == "Inf")] <- 0.00
    df_relative[k,] <- c(m.kappa[k,], round(eff, 2))
    }
}


###########################################################################################
################# Efficiency values of (tau2 = 1/q)-optimal designs w.r.t. (tau2 = 1) criteria
###########################################################################################

load("output/workspaces/MSEDP_point_Q.RData")

## set tau2 to 1 -- to obtain the relative efficiency values
tau2 <- 1.0
tau <- sqrt(tau2)

for (k in 1:nrow(m.kappa)) {
  S = S.mseD.point[[k]]
  
  if(!is.null(S)) {
    
    kappa.DP <- m.kappa[k,1]; kappa.LoF <- m.kappa[k,2]; kappa.mse <- m.kappa[k,3]
    if (!all(m.kappa[k,] == S$kappa)) {print("Kappa-s"); next}
    
    eval <- criteria.mseD.point(S$X1, S$X2)  # design evaluation
    out <- output(S)
    eff <- as.numeric(Re(c(crit.optim*100/c(out$values$DP, out$values$LoFDP, out$mse.point),
                           MSEDP.output$compound[k]*100/eval$compound)))
    eff[which(eff == "Inf")] <- 0.00
    df_relative_Q[k,] <- c(m.kappa[k,], round(eff, 2))
  }
}

df_relative$tau2 = "1"                # the 'original \tau^2'
df_relative_Q$tau2 = "1/q"
df_relative <- rbind(df_relative, df_relative_Q)

df_relative <- dplyr::rename(df_relative, "DPs" = "DP", "MSE_Ds" = "mseD")

###################### Plotting ###############################


# reshaping the data for plotting
data_long <- tidyr::gather(df_relative, key = "variable", value = "value", DPs, LoFDP, MSE_Ds, compound)

data_long <- dplyr::rename(data_long, "DPs" = "kappa.DP", "LoFDP" = "kappa.LoFDP",
                           "MSE_Ds" = "kappa.mseD")

  
# Create a function to generate simplex plots
generate_simplex_plot <- function(data, criterion, t_val) {
  data_filtered <- data %>%
    filter(tau2 == t_val) %>%
    filter(pull(., variable) == criterion)
  
    l_value <- ifelse(t_val == "1", "(A)", "(B)")

    ggtern(data = data_filtered, aes(x = DPs, y = LoFDP,
                         z = MSE_Ds, label = sprintf("%.2f", value))) +
    geom_mask() +
    geom_point(data=subset(data_filtered, value!=0),
               aes(shape = 'Not NA', size = value, color = -value), alpha = 0.9) +
    scale_fill_viridis_c()+
    geom_point(data=subset(data_filtered, value == 0), 
               aes(shape = 'NA'), color = "grey30", size = 1, alpha = 0.8) + 
    scale_shape_manual(values=c('NA'= 19, 'Not NA'=19)) +
    scale_size(range = c(0, 11), limits = c(0.0, 102)) +
    geom_label(hjust=+1.0, vjust=-0.25) + 
    theme(legend.position = "none",
          tern.axis.title.T = element_text(size = 12, hjust = -.1, vjust = 0.5),
          tern.axis.title.L = element_text(size = 12, hjust = .4, vjust = 1.25),
          tern.axis.title.R = element_text(size = 12, hjust = 0.75, vjust = 1.25),
          plot.title = element_text(size=14, hjust=0.15, vjust = -9, face = "bold")) +
    labs(title = paste(l_value, criterion, sep = " ")) + 
    theme_showarrows() +
    theme_clockwise() 

}

# Create a grid of simplex plots
gg <- vector(mode = "list", length = 8)

# For each criterion and each value of tau2
criteria <- c("DPs", "LoFDP", "MSE_Ds", "compound")
t_values <- c("1", "1/q")
k<-1

for (t_val in t_values) {
  for (criterion in criteria) {
    gg[[k]] <- generate_simplex_plot(data_long, criterion, t_val)
    k <- k+1
  }
}

gg_all <- ggtern::grid.arrange(grobs = gg, 
                               nrow = 4, ncol = 2, as.table = FALSE)

## Saving the plot
ggsave(gg_all, filename ="output/figures/Figure_relative_efficiencies.jpeg", height = 16, width = 9)

