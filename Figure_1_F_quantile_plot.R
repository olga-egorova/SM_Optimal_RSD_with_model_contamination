
###### Plot the inverse F-quantiles for different values of pure error and  
###### lack-of-fit degrees of freedom  


library(data.table)
library(ggplot2)


np0<-15; d0<-1:(np0-1)
np1<-19; d1<-1:(np1-1)
np2<-25; d2<-1:(np2-1)
np3<-30; d3<-1:(np3-1)
prob.LoF<-0.95;

n0<-rep(np0,length(d0)); n1<-rep(np1,length(d1)); 
n2<-rep(np2,length(d2)); n3<-rep(np3,length(d3))

df<-data.table("d" = c(d0, d1, d2, d3),
               "np" = c(n0, n1, n2, n3))

df$lof<-1./qf(prob.LoF, df$np-df$d,df$d)
df$np = as.factor(df$np)

axis_title_size = 18
axis_title_family = "mono"
axis_text_family = "mono"
strip_font_family = "sans"
axis_text_size = 15
ticks_length = -1.5
legend_title_size = 16

plot_width = 10
plot_height = 6

ggplot(df, aes(x=d, y = lof)) +
  geom_line(aes(colour = np, linetype = np), lwd = 1) +
  xlab(paste0(expression(d),", pure error d.f.")) +
  ylab(paste(expression(1/F(n-p-d,d)))) +
  scale_linetype_manual(values =c(5,3,1,4),
                        labels=c("n-p = 15", "n-p = 19", "n-p = 25", "n-p = 30")) +
  scale_color_manual(
                     labels=c("n-p = 15", "n-p = 19", "n-p = 25", "n-p = 30"),
                     values=c("#0072B2","#009E73","#D55E00", "#999999")) + 
  labs(color  = "Residual d.f.", linetype = "Residual d.f.") +
  theme_bw() +
  theme(aspect.ratio=0.75, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.border = element_rect(colour = "black", fill = "NA", size = 1.0),
        axis.text.x = element_text(size = axis_text_size, family = axis_text_family, colour = "black",
                                   margin = unit(c(t = 3.0, r = 0, b = 2.5, l = 0), "mm")),
        # adjust X- and Y-axis title
        axis.title.y = element_text(family = axis_title_family, size = axis_title_size), 
        axis.title.x = element_text(family = axis_title_family, size = axis_title_size),
        # adjust Y-axis labels
        axis.text.y.left = element_text(size = axis_text_size, family = axis_text_family, colour = "black",
                                        margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")),
        axis.ticks.length = unit(ticks_length, "mm"), #axis.ticks.length.y.left = unit(ticks_length, "mm"),
        legend.position = c(0.75, 0.25),
        strip.background = element_rect(color = "black", fill="white", size = 1.0),
        strip.text = element_text(family = strip_font_family, size = 13),
        legend.text = element_text(family = axis_text_family, size = axis_title_size-3),
        legend.spacing.y = unit(0.35, 'in'),
        legend.title = element_blank())

ggsave("F_quantile_plot.jpeg", width = plot_width, height = plot_height, units = "in")

