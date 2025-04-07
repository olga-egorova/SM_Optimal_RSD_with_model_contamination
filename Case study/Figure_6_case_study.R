################################################################################################
### Generating a figure with efficiency values of MSE(DP)-optimal designs for the blocked experiment
################################################################################################

library(glue)
rm(list=ls()) 

### Optional -- run the search for the optimal designs. These can be run in parallel.
### This can take some time -- up to 5 hours each.
### Alternatively -- skip this and move to the next line

#source("MSEB_main.R")
#source("MSEB_Q_main.R")
#source("MSEB_CP_main.R")
#source("MSEB_CP_Q_main.R")


### Load the table containing the criteria values of the designs

MSEB.output <- read.csv("output/tables/MSEB_values.csv", header = TRUE)
MSEB_Q.output <- read.csv("output/tables/MSEB_Q_values.csv", header = TRUE)
MSEB_CP.output <- read.csv("output/tables/MSEB_CP_values.csv", header = TRUE)
MSEB_CP_Q.output <- read.csv("output/tables/MSEB_CP_Q_values.csv", header = TRUE)

## tau2 = 1
MSEB_eff <- rbind(MSEB_CP.output, MSEB_CP.output)
MSEB_eff$efficiency <- c(rep("No fixed CP", 6), rep("Fixed CP", 6))
MSEB_eff[1:6, c("DP", "LoFDP", "mseD")] <- t(apply(MSEB_eff[1:6, c("DP", "LoFDP", "mseD")], 
                                                   MARGIN = 1, 
                                                   FUN = function(x) diag(as.matrix(MSEB.output[4:6, c("DP", "LoFDP", "mseD")]))/x*100)
)  ## "No CP" individual efficiency values
MSEB_eff[7:12, c("DP", "LoFDP", "mseD")] <- t(apply(MSEB_eff[7:12, c("DP", "LoFDP", "mseD")], 
                                                   MARGIN = 1, 
                                                   FUN = function(x) (diag(as.matrix(MSEB_CP.output[4:6, c("DP", "LoFDP", "mseD")])))/x*100)
) ## "CP" individual efficiency values
## compound efficiency values
MSEB_eff$compound <- c(MSEB.output$compound, MSEB_CP.output$compound)/MSEB_eff$compound*100 


## tau2 = 1/Q
MSEB_Q_eff <- rbind(MSEB_CP_Q.output, MSEB_CP_Q.output)
MSEB_Q_eff$efficiency <- c(rep("No fixed CP", 6), rep("Fixed CP", 6))
MSEB_Q_eff[1:6, c("DP", "LoFDP", "mseD")] <- t(apply(MSEB_Q_eff[1:6, c("DP", "LoFDP", "mseD")], 
                                                     MARGIN = 1, 
                                                     FUN = function(x) diag(as.matrix(MSEB_Q.output[4:6, c("DP", "LoFDP", "mseD")]))/x*100)
)  ## "No CP" individual efficiency values
MSEB_Q_eff[7:12, c("DP", "LoFDP", "mseD")] <- t(apply(MSEB_Q_eff[7:12, c("DP", "LoFDP", "mseD")], 
                                                      MARGIN = 1, 
                                                      FUN = function(x) (diag(as.matrix(MSEB_CP_Q.output[4:6, c("DP", "LoFDP", "mseD")])))/x*100)
) ## "CP" individual efficiency values
## compound efficiency values
MSEB_Q_eff$compound <- c(MSEB_Q.output$compound, MSEB_CP_Q.output$compound)/MSEB_Q_eff$compound*100 

MSEB_eff$tau2 = "1"
MSEB_Q_eff$tau2 = "1/q"

df_eff <- rbind(MSEB_eff, MSEB_Q_eff)
df_eff[, c("DP", "LoFDP", "mseD", "compound")] <- 
  round(df_eff[, c("DP", "LoFDP", "mseD", "compound")], 2)


### Table with efficiency values of the designs

print(df_eff, digits = 4)

write.csv(df_eff, "output/tables/MSEB_efficiency.csv")
print(xtable::xtable(df_eff, file = "output/tables/MSEB_efficiency.tex"))




#################################################################################
####################### Plot ####################################################
#################################################################################



df_eff <- read.csv("output/tables/MSEB_efficiency.csv", header = TRUE, sep=",")

df_eff <- tidyr::pivot_longer(df_eff, cols = c("DP", "LoFDP", "mseD", "compound"),
                              names_to = "criterion", 
                              values_to = "efficiency_value")

df_eff$weights = paste("[", df_eff$kappa.DP, ",", 
                       df_eff$kappa.LoFDP, ",", 
                       df_eff$kappa.mseD, "]", sep = "")
df_eff$weights[which(abs(df_eff$kappa.DP - 1./3) < 0.01)] <- "[1/3, 1/3, 1/3]"
df_eff$weights[which(abs(df_eff$kappa.DP - 0.25) < 0.01)] <- "[1/4, 1/4, 1/2]"
df_eff$weights[which(abs(df_eff$kappa.DP - 0.4) < 0.01)] <- "[2/5, 1/5, 2/5]"

df_plot <- df_eff[df_eff$X %in% 1:3,]


ggplot2::ggplot(df_plot, aes(x = criterion, y = efficiency_value, 
                             color = weights, shape = weights, group = weights)) +
  geom_line(linewidth = 1.1, alpha = 0.6) +
  geom_point(size = 5, alpha = 0.6) +
  facet_grid(glue('tau^2*" = {tau2}"') ~ glue('"{efficiency}"'),
             labeller = label_parsed) + 
  labs(x = "", y = "Efficiency, %", color = "Weights\n [DP, LoFDP, MSE(D)]", shape = "Weights\n [DP, LoFDP, MSE(D)]") +
  scale_color_brewer(palette = "Set1") +  # Colorblind-friendly palette
  theme_minimal(base_size = 16) +  # Increase base font size for better readability
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle x-axis labels for readability
    axis.text.y.right = element_text(angle = 90, vjust = 0, hjust = 1),  # Angle y-axis labels for readability
    axis.title = element_text(size = 18),  # Increase size of axis titles
    axis.text = element_text(size = 14),   # Increase size of axis ticks
    strip.text = element_text(size = 16),  # Increase size of facet labels
    plot.title = element_text(size = 20, face = "bold"),  # Increase size of plot title
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
    legend.title = element_text(size = 14),  # Increase size of legend title
    legend.text = element_text(size = 14)   # Increase size of legend text
   )

ggsave(filename ="output/Figure_MSEDB.jpeg", height = 7, width = 12)



