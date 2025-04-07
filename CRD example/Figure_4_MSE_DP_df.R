################################################################################################
### Generating a figure with pure error and lack-of-fit d.f. of MSE(DP)-optimal (compound) designs
################################################################################################

library(ggtern)
library(data.table)
rm(list=ls()) 

### Load the table containing the criteria values of the designs

MSEDP.output <- read.csv("output/tables/MSEDP_point_values.csv", header = TRUE)
MSEDP.output.Q <- read.csv("output/tables/MSEDP_point_Q_values.csv", header = TRUE)

df_df<- cbind(Re(as.matrix(MSEDP.output[, 2:6])), 
              Re(as.matrix(MSEDP.output.Q[, 5:6])))
names(df_df) <- c("DPs", "LoFDP", "MSE_Ds", "(a) PE", "(a) LoF", "(b) PE", "(b) LoF")


gg <- vector(mode = "list", length = 4)

for (k in 1:4) {
  
  des <- data.table(df_df[1:10, c(1:3, 3+k)])
  
  names(des) <- c("DPs", "LoFDP", "MSE_Ds", "y")
  
  gg[[k]] <- ggtern(data = des, aes(x=DPs, y=LoFDP,
                                    z = MSE_Ds, label = y)) +
    geom_mask() +
    geom_point(data=subset(des, y!=0),
               aes(shape = 'Not NA', size = y, color = -y), alpha = 0.9) +
    scale_fill_viridis_c()+
    geom_point(data=subset(des, y==0), 
               aes(shape = 'NA'), color = "grey30", size = 1, alpha = 0.8) + 
    scale_shape_manual(values=c('NA'= 19, 'Not NA'=19)) +
    scale_size(range = c(0, 11), limits = c(0.0, 22)) +
    geom_label(hjust=+1.0, vjust=-0.25) + 
    theme(legend.position = "none",
          tern.axis.title.T = element_text(size = 12, hjust = -.1, vjust = 0.5),
          tern.axis.title.L = element_text(size = 12, hjust = .45, vjust = 1.4),
          tern.axis.title.R = element_text(size = 12, hjust = 0.95, vjust = 1.4),
          plot.title = element_text(size=12, hjust=0.05, vjust = -9, face = "bold")) +
    labs(title = names(df_df[k+3])) + 
    theme_showarrows() +
    theme_clockwise() 
  
}

gg_all <- ggtern::grid.arrange(grobs=gg, 
                               nrow = 2, ncol = 2, 
                               as.table = FALSE)


ggsave(gg_all, filename ="output/figures/Figure_MSE_DP_df.jpeg", height = 6, width = 10)

