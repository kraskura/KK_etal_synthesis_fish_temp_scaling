
# Description: 
#   Supplemental material:
#     - histograms 
#     - correlations between MR metrics
# ********************************************
# ********************************************



# Figure 1 ----
# 
library(here)

data<-readxl::read_excel("./Data/Summary_scaling_b_mmr_rmr_mar2026.xlsx")
data<-data %>% 
  mutate(Order = row_number(desc(MMR)),
         ID = paste(Species_common, Test_C, Study, "")) %>% 
  pivot_longer(cols = c(MMR, RMR), names_to = "MR", values_to = "b")

# data$Running_ID<-factor(data$Running_ID)
# data$star<-NA
# for(i in 1:nrow(data)){
#   if(is.na(data$Order_MMR[i])){
#     data$Order_MMR[i]<- data[which(data$Running_ID[i] == data$Running_ID)[1],"Order_MMR"]
#   }
#   if(!data$notes.mechan[i] == "REPEAT" & !data$notes.mechan[i] == ""){
#     data$star[i]<-"•"
#   }
# }

# data

plot <- ggplot(data = data, mapping = aes(x = b, y=factor(Order), size = MR,
                                    Color = Test_C, fill = Test_C,
                                    shape = MR, group = ID))+
  geom_vline(xintercept = c(1, 0.75), lty = "dashed", color="black", linewidth= 0.3)+
  geom_path(color = "grey",
            arrow = arrow(length = unit(0.1, "cm"),
                          angle = 1), 
            show.legend = F)+
  geom_point(alpha=1)+
  xlim(0.01, 1.55)+
  geom_text(aes(label = Study, x = 0.01), family = "Helvetica", size=3, hjust = 0)+
  scale_color_viridis_c(option = "A")+
  scale_fill_viridis_c(option = "A", name = "T, ºC")+
  xlab(expression(Slope~value~(italic(b))))+
  scale_shape_manual(values = c(21, 21), name = "")+
  scale_size_manual(values = c(3, 1), name = "")+
    scale_y_discrete(
      breaks = data$Order, 
      labels = data$Species_common)+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 9.5),
        axis.text.x = element_text( color = "black", size =10),
        axis.title.y = element_blank(),
        legend.position = c(0.88, 0.1),
        axis.line.y=element_line(colour = 'black',size=0.5),
        axis.line.x=element_line(colour = 'black',size=0.5),
        axis.ticks.y=element_line(size=0.5),
        axis.ticks.x=element_line(size=0),
        text=element_text(size=12,  family="Helvetica"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.margin = margin(-0.3,0,0,0, unit="cm"), legend.box = "horizontal")+
  guides(size = "legend", colour = "legend") 

  plot
  
ggsave(filename = paste("./Figures/Figure1_main.png", sep=""),
       width = 6, height = 7.5, units = "in")
 

# Figure 2 ------------
# conceptual figures
data_sim<-read.csv(here("./Data/Conceptual_fig_data.csv"))

psim<-ggplot(data_sim[data_sim$simID==3,], aes(x=log(BW_kg), y=log(pred.mmr.mgO2min) ))+
  geom_abline(slope = data_sim[data_sim$simID==3,"slope.MMR"][1],
              intercept = data_sim[data_sim$simID==3,"int.MMR"][1],
              color = "black", linewidth=0.7)+ # MMR
  geom_abline(slope = data_sim[data_sim$simID==3,"slope.SMR"][1],
              intercept = data_sim[data_sim$simID==3,"int.SMR"][1],
              color = "#7e959b", linewidth=0.7)+ # MMR
  lims(x = c(-4.5, 4.5), y = c(-4.5, 4.5))+
  annotate("text",  x = -4.3, y = 4.4, label = "A", fontface =2,
           size=4.5, hjust=1, family="Helvetica", color = "black")+
  annotate("text",  x = 0.8, y = 3.5, label = expression(paste("MMR")),
       size=5, hjust=0, family="Helvetica", angle = 42, color = "black")+
  annotate("text",  x = 1, y = 0.9, label = expression(paste("RMR")),
           size=5, hjust=0, family="Helvetica", angle = 33, color = "#7e959b")+
  annotate("text",  x = -3.5, y = -4.2, label = expression(paste("OPTIMAL TEMPERATURES")),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = -3.8, y = 4.3, label = expression(paste("HYPOTHESIS")),
           size=4, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = -3.8, y = 3.6, label = expression(italic(b)[MMR]~">"~italic(b)[RMR]),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 1.75, y = -0.8,
           label = expression(italic(b)[MMR]~"="~"1.00"),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 1.75, y = -1.5,
           label = expression(italic(b)[RMR]~"="~0.89),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 1.75, y = -2.2,
           label = expression(italic(b)[AS]~"="~1.05),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 1.75, y = -2.9,
           label = expression(italic(b)[FAS]~"="~0.11),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  xlab(expression(italic(ln)*Body~mass~(g)))+
  ylab(expression(italic(ln)*MR~(mg~O[2]~h^-1)))+
  theme(axis.text.y=element_text(size=12, colour= 'black'),
		axis.text.x=element_text(size=12, colour= 'black'),
		axis.line.y=element_line(colour = 'black',linewidth=0.5),
		axis.line.x=element_line(colour = 'black',linewidth=0.5),
		axis.ticks.y=element_line(linewidth=0.5),
		panel.background = element_blank(),
		axis.ticks.x.bottom = element_line(linewidth=0.5, colour = "black"),
	  axis.title.y=element_text(size=12),
		axis.title.x=element_text(size=12),
		panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))
# psim

# conceptual Warm
psim_w<-ggplot(data_sim[data_sim$simID==1,], aes(x=log(BW_kg), y=log(pred.mmr.mgO2min) ))+
  geom_abline(slope = data_sim[data_sim$simID==1,"slope.MMR"][1],
              intercept = data_sim[data_sim$simID==1,"int.MMR"][1],
              color = cols.amr[2], linewidth=0.7)+ # MMR
  geom_abline(slope = data_sim[data_sim$simID==1,"slope.SMR"][1],
              intercept = data_sim[data_sim$simID==1,"int.SMR"][1],
              color = cols.rmr[2], linewidth=0.7)+ # MMR
  lims(x = c(-4.5, 4.5), y = c(-4.5, 4.5))+
  annotate("text",  x = -4.3, y = 4.4, label = "B", fontface =2,
           size=4.5, hjust=1, family="Helvetica", color = "black")+
  annotate("text",  x = 0.8, y = 3.3, label = expression(paste("MMR")),
       size=5, hjust=0, family="Helvetica", angle = 32, color = cols.amr[2])+
  annotate("text",  x = 1, y = 0.9, label = expression(paste("RMR")),
           size=5, hjust=0, family="Helvetica", angle = 47, color = cols.rmr[2])+
  annotate("text",  x = -3.3, y = -4.2, label = expression(paste("WARM TEMPERATURES")),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = -3.8, y = 4.3, label = expression(paste("HYPOTHESIS")),
           size=4, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = -3.8, y = 3.6, label = expression(italic(b)[MMR]~"<"~italic(b)[RMR]),
           size=4, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 1.75, y = -0.8,
           label = expression(italic(b)[MMR]~"="~0.75),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 1.75, y = -1.5,
           label = expression(italic(b)[RMR]~"="~0.95),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 1.75, y = -2.2,
           label = expression(italic(b)[AS]~"="~0.69),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 1.75, y = -2.9,
           label = expression(italic(b)[FAS]~"="~-0.14),
           size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
  xlab(expression(italic(ln)*Body~mass~(g)))+
  ylab(expression(italic(ln)*MR~(mg~O[2]~h^-1)))+
  theme(axis.text.y=element_text(size=12, colour= 'black'),
		axis.text.x=element_text(size=12, colour= 'black'),
		axis.line.y=element_line(colour = 'black',linewidth=0.5),
		axis.line.x=element_line(colour = 'black',linewidth=0.5),
		axis.ticks.y=element_line(linewidth=0.5),
		panel.background = element_blank(),
		axis.ticks.x.bottom = element_line(linewidth=0.5, colour = "black"),
	  axis.title.y=element_text(size=12),
		axis.title.x=element_text(size=12),
		panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))

# plot_layout(widths = c(2, 1), heights = unit(c(5, 1), c('cm', 'null')))

# save the plots
mainplot<-
  cowplot::plot_grid(psim, psim_w,
                   align = "hv", ncol = 2)

ggsave(filename = paste("./Figures/Figure2.png", sep=""),
       plot=mainplot, width = 6.8, height = 3.3, units = "in")

