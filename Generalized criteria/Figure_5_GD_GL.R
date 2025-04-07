
#### Efficiency values of GD and GL-optimum designs, across the values of tau2 and 
#### numbers of added center points


library(data.table)
library(ggplot2)
library(ggpattern)
library(viridis)


#source("GD_4CP_main.R"); source("GD_4CP_Q_main.R")
#source("GD_5CP_main.R"); source("GD_5CP_Q_main.R")

#source("GL_4CP_main.R"); source("GL_4CP_Q_main.R")
#source("GL_5CP_main.R"); source("GL_5CP_Q_main.R")


#### Either run the lines above to search for the GD- and GL-optmal designs, or
#### read the criteria values:


GD_4CP <- read.csv("output/tables/GD_4CP.csv", header = TRUE, sep = ",")
GD_4CP_Q <- read.csv("output/tables/GD_4CP_Q.csv", header = TRUE, sep = ",")
GD_5CP <- read.csv("output/tables/GD_5CP.csv", header = TRUE, sep = ",")
GD_5CP_Q <- read.csv("output/tables/GD_5CP_Q.csv", header = TRUE, sep = ",")

GL_4CP <- read.csv("output/tables/GL_4CP.csv", header = TRUE, sep = ",")
GL_4CP_Q <- read.csv("output/tables/GL_4CP_Q.csv", header = TRUE, sep = ",")
GL_5CP <- read.csv("output/tables/GL_5CP.csv", header = TRUE, sep = ",")
GL_5CP_Q <- read.csv("output/tables/GL_5CP_Q.csv", header = TRUE, sep = ",")


##### Combining the data for plotting

GD_all <- data.table(rbind(GD_4CP, GD_4CP_Q, GD_5CP, GD_5CP_Q))
GL_all <- data.table(rbind(GL_4CP, GL_4CP_Q, GL_5CP, GL_5CP_Q))

dt_compound <- data.table(rbind(GD_all[, -c(9:11)], GL_all[, -c(9:11)]))
dt_compound$compound <- round(Re(dt_compound$compound)*100, 2)
dt_compound$setting <- paste(dt_compound$CPs, ", ", "τ² = ", dt_compound$tau2, sep = "")


############################################################
#################### Making the plot #######################
############################################################


ggplot(dt_compound, aes(x = setting, y = compound, 
                        color = setting)) +
  geom_boxplot_pattern(aes(pattern = setting, pattern_fill = setting),
                       pattern_density = 0.35) + 
  facet_wrap(
    ~ criterion, 
    scales = "free_x", 
    ncol = 2,  
    labeller = labeller(criterion = c("GD" = "GD-optimal designs", "GL" = "GL-optimal designs"))  
  ) +
  theme_minimal(base_size = 18) +  
  theme(
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),   # Increase x-axis values size
    axis.text.y = element_text(size = 14),   # Increase y-axis values size
    #panel.grid = element_blank(),  # Remove gridlines
    legend.position = "none", 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Add a black border around the plot
  ) + 
  labs(
    x = "", 
    y = "Efficiency, %"  
  ) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#228B22", "#8B0000")) + 
  scale_pattern_manual(values = c("stripe", "wave", "crosshatch", "weave")) + 
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#228B22", "#8B0000")) +  
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))


ggsave(filename = "output/Figure_GD_GL.jpeg", width = 8, height = 6)
