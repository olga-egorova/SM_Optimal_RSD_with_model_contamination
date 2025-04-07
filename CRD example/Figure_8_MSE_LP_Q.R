################################################################################################
### Generating Table with efficiency values of MSE(LP)-optimal (compound) designs
################################################################################################

library(ggtern)
rm(list=ls()) 

### Optional -- run the search for the optimal designs
### This can take some time -- up to 10 hours
### Alternatively -- skip this and move to the next line

#source("MSE_LP_Q_main.R")

### Load the table containing the criteria values of the designs
MSELP_Q.output <- read.csv("output/tables/MSELP_Q_values.csv", header = TRUE)
MSEDP_Q.output <- read.csv("output/tables/MSEDP_point_Q_values.csv", header = TRUE)

### Optimal values of the component criteria
crit.optim.Q <- list("LP" = MSELP_Q.output[1, "LP"],
                     "LoFLP" = MSELP_Q.output[2, "LoFLP"],
                     "MSE_L" = MSELP_Q.output[3, "mseL"],
                     "DP" = MSEDP_Q.output[1, "DP"],
                     "LoFDP" = MSEDP_Q.output[2, "LoFDP"],
                     "MSE_D" = MSEDP_Q.output[3, "mseD.point"]
                   )

crit.optim.Q <- as.numeric(Re(unlist(crit.optim.Q)))

### Matrix of efficiency values in %
df.efficiency.Q <- Re(as.matrix(MSELP_Q.output[,c(1:11, 13)]))

for (k in 1: nrow(df.efficiency.Q)) {
 df.efficiency.Q[k, 7:12] <- ifelse(df.efficiency.Q[k, 7:12] == 0,
                                    0, crit.optim.Q*100./df.efficiency.Q[k, 7:12])
}

### The table of designs' efficiency values and d.f. distribution
print(df.efficiency.Q[,-1], digits = 4)

### Table with efficiency values
write.csv(df.efficiency.Q[,-1], "output/tables/MSELP_eff_Q.csv")
print(xtable::xtable(df.efficiency.Q), file = "output/tables/MSELP_Q_eff.tex")

##########################################################################################
####################### Simplex plots ####################################################

df_eff <- read.csv("output/tables/MSELP_eff_Q.csv", header = TRUE, sep=",")[1:10, c(2:4, 7:12)]
names(df_eff) <- c("LPs", "LoFLP", "MSE_Ls", "LPs-eff", "LoFLP-eff", "MSE_Ls-eff",
                   "DPs-eff", "LoFDP-eff", "MSE_Ds-eff")

gg <- vector(mode = "list", length = 6)

for (k in 1:6) {
  
  des <- df_eff[1:10, c(1:3, 3+k)]
  
  names(des) <- c("LPs", "LoFLP", "MSE_Ls", "y")
  
  gg[[k]] <- ggtern(data = des, aes(x=LPs, y=LoFLP,
                                    z = MSE_Ls, label = sprintf("%.2f", y))) +
    geom_mask() +
    geom_point(data=subset(des, y!=0),
               aes(shape = 'Not NA', size = y, color = -y), alpha = 0.9) +
    scale_fill_viridis_c()+
    geom_point(data=subset(des, y==0), 
               aes(shape = 'NA'), color = "grey30", size = 1, alpha = 0.8) + 
    scale_shape_manual(values=c('NA'= 19, 'Not NA'=19)) +
    scale_size(range = c(0, 11), limits = c(0.0, 101)) +
    geom_label(hjust=+1.0, vjust=-0.25) + 
    theme(legend.position = "none",
          tern.axis.title.T = element_text(size = 12, hjust = -.1, vjust = 0.5),
          tern.axis.title.L = element_text(size = 12, hjust = .4, vjust = 1),
          tern.axis.title.R = element_text(size = 12, hjust = 0.7, vjust = 1),
          plot.title = element_text(size=14, hjust=0.15, vjust = -9, face = "bold")) +
    labs(title = names(df_eff[k+3])) + 
    theme_showarrows() +
    theme_clockwise() 
  
}

gg_all <- ggtern::grid.arrange(grobs=gg, 
                               labels = c("LP", "LoFLP", "MSE_L",
                                          "DP", "LoFDP", "MSE_D"),
                               nrow = 3, ncol = 2, as.table = FALSE)

ggsave(gg_all, filename ="output/figures/Figure_MSELP_Q.jpeg", height = 13, width = 11)
