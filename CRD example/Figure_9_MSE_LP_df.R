################################################################################################
### Generating a figure with pure error and lack-of-fit d.f. of MSE(LP)-optimal (compound) designs
################################################################################################

library(ggtern)
rm(list=ls()) 

### Load the table containing the criteria values of the designs

MSELP.output <- read.csv("output/tables/MSELP_values.csv", header = TRUE)
MSELP.output.Q <- read.csv("output/tables/MSELP_Q_values.csv", header = TRUE)

df_df <- data.frame(cbind(Re(as.matrix(MSELP.output[, 2:6])), Re(as.matrix(MSELP.output.Q[, 5:6]))))
names(df_df) <- c("LPs", "LoFLP", "MSE_Ls", "(a) PE", "(a) LoF", "(b) PE", "(b) LoF")
gg <- vector(mode = "list", length = 4)

for (k in 1:4) {
  
  des <- df_df[1:10, c(1:3, 3+k)]
  
  names(des) <- c("LPs", "LoFLP", "MSE_Ls", "y")
  
  gg[[k]] <- ggtern(data = des, aes(x=LPs, y=LoFLP,
                                    z = MSE_Ls, label = y)) +
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
          tern.axis.title.T = element_text(size = 12, hjust = -.1, vjust = 0.35),
          tern.axis.title.L = element_text(size = 12, hjust = .45, vjust = 1.42),
          tern.axis.title.R = element_text(size = 12, hjust = 0.95, vjust = 1.42),
          plot.title = element_text(size=12, hjust=0.05, vjust = -9, face = "bold")) +
    labs(title = names(df_df[k+3])) + 
    theme_showarrows() +
    theme_clockwise() 
  
}

gg_all <- ggtern::grid.arrange(grobs=gg, 
                               nrow = 2, ncol = 2, 
                               as.table = FALSE)


ggsave(gg_all, filename ="output/figures/Figure_MSE_LP_df.jpeg", height = 6, width = 10)

