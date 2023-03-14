
# import data, source local scripts ---------
# source("./R/phylo_mixed_model.R") # needed libraries called within 
# source("./R/nonPhylo_mixed_models.R") # for comparison supplemental 

source("./R/get_data_species_ecologies.R")

# colors --------
cols.as<<-c("#265F73", "#007E66", "#00C5A3")
cols.fas<<-c("#395200", "#89A000", "yellow")
cols.rmr<<-c("#C70039", "#FF6D7C", "#FFA3AC")
cols.amr<<-c("#00749F","#00A8D6", "#9CE9FF")
cols<<-c("#00749F","#C70039","#00A8D6","#FF6D7C", "#9CE9FF","#FFA3AC", "#00C5A3", "#265F73")# AMR -rmr- AMR dark - rmr dark - light - as - fas
 

# local setup ------
set.seed(51423)
# MMR and AMR used interchangeably throughout 

library(ggplot2)
library(ggpubr)
library(cowplot)
library(forcats)
library(here)

  
# Supplemental -----------
## Histograms of fish at different size suppl ----------
  
plot_hist_amr<-ggplot(data.amr, aes(x=BW_g)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  # geom_vline(xintercept = 300, color="grey", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = 1000, color="grey", lty="dashed", linewidth=0.5)+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.amr[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,300)+
  scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_amr, x_title = "Body mass (g)", y_title = "Frequency", print = FALSE)
plot_hist_amr <- plot_hist_amr + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))

plot_hist_rmr<-ggplot(data.rmr, aes(x=BW_g)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  # geom_vline(xintercept = 300, color="grey", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = 1000, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.rmr[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,300)+
  scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_rmr, x_title = "Body mass (g)", y_title = "Frequency", print = FALSE)
plot_hist_rmr <- plot_hist_rmr + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                       axis.text.x = element_text(size=21,  family="Arial"),
                                       axis.title.y = element_text(size=25,  family="Arial"),
                                       axis.text.y = element_text(size=21,  family="Arial"))

plot_hist_fas<-ggplot(data.fas, aes(x=BW_g)) +
  # geom_vline(xintercept = mean(d_bind.fas$MLEslope[d_bind.fas$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  # geom_vline(xintercept = 300, color="grey", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = 1000, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.fas[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,300)+
  scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_fas, x_title = "Body mass (g)", y_title = "Frequency", print = FALSE)
plot_hist_fas <- plot_hist_fas + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                       axis.text.x = element_text(size=21,  family="Arial"),
                                       axis.title.y = element_text(size=25,  family="Arial"),
                                       axis.text.y = element_text(size=21,  family="Arial"))


## Sum plots: histograms of fish at different temperatures suppl ----------

plot_hist_amr_T<-ggplot(data.amr, aes(x=tempTest)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  # geom_vline(xintercept = 300, color="grey", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = 15, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.amr[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  # ylim(0,300)+
  # scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_amr_T, x_title = expression(Temperature~degree*C), y_title = "Frequency", print = FALSE)
plot_hist_amr_T <- plot_hist_amr_T + theme(axis.title.x = element_text(size=21,  family="Arial"),
                                       axis.text.x = element_text(size=21,  family="Arial"),
                                       axis.title.y = element_text(size=21,  family="Arial"),
                                       axis.text.y = element_text(size=21,  family="Arial"))

plot_hist_rmr_T<-ggplot(data.rmr, aes(x=tempTest)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  # geom_vline(xintercept = 300, color="grey", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = 15, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.rmr[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  # ylim(0,300)+
  # scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_rmr_T, x_title =  expression(Temperature~degree*C), y_title = "Frequency", print = FALSE)
plot_hist_rmr_T <- plot_hist_rmr_T + theme(axis.title.x = element_text(size=21,  family="Arial"),
                                       axis.text.x = element_text(size=21,  family="Arial"),
                                       axis.title.y = element_text(size=21,  family="Arial"),
                                       axis.text.y = element_text(size=21,  family="Arial"))

plot_hist_fas_T<-ggplot(data.fas, aes(x=tempTest)) +
  # geom_vline(xintercept = mean(d_bind.fas$MLEslope[d_bind.fas$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  # geom_vline(xintercept = 300, color="grey", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = 15, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.fas[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  # ylim(0,300)+
  # scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_fas_T, x_title =  expression(Temperature~degree*C), y_title = "Frequency", print = FALSE)
plot_hist_fas_T <- plot_hist_fas_T + theme(axis.title.x = element_text(size=21,  family="Arial"),
                                       axis.text.x = element_text(size=21,  family="Arial"),
                                       axis.title.y = element_text(size=21,  family="Arial"),
                                       axis.text.y = element_text(size=21,  family="Arial"))

hist.plot<-cowplot:::plot_grid(plot_hist_amr, plot_hist_amr_T, 
                               plot_hist_rmr, plot_hist_rmr_T,
                               plot_hist_fas, plot_hist_fas_T,
                               align = "hv",
                               axis = "l",
                               nrow = 3,
                               ncol = 2) 
save_plot(filename = paste("./Figures/Suppl_Size_hist", Sys.Date(), ".png", sep=""),
          hist.plot, base_width = 12, base_height = 20, units = "in")



# Plots: Correlations metabolic rates ---------
MRcorrel1.small<-ggplot(data.amr[data.amr$BW_g<100,], aes(x=RMR, y=AMR, group=test_category3, color = tempTest, shape=test_category3, fill=tempTest)) +
  geom_point(alpha=0.6, pch=21,  show.legend = FALSE, size=2)+
  scale_fill_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_color_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_shape_manual(values= c(22,23))+
  facet_grid(.~test_category3)+
  ylim(0, 75)
ggformat(MRcorrel1.small, title = "",y_title = expression(MMR~(mg~O[2]~h^-1)), x_title = expression(RMR~(mg~O[2]~h^-1)), print = F)

MRcorrel1.large<-ggplot(data.amr, aes(x=RMR, y=AMR, group=test_category, fill=tempTest, shape=test_category)) +
  geom_point(alpha=0.6,shape=21,  show.legend = FALSE, size=3)+
  geom_point(data.amr[data.amr$species=="Somniosus microcephalus",],  mapping=aes(x=RMR, y=AMR, group=test_category, fill=tempTest, shape=test_category), shape=23,  show.legend = FALSE, size=3, stroke=2)+
  scale_color_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_fill_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_shape_manual(values= c(22,23,24))+
  facet_grid(.~test_category3)+
  ylim(0, 3500)
ggformat(MRcorrel1.large, title = "", y_title = expression(MMR~(mg~O[2]~h^-1)), x_title = expression(RMR~(mg~O[2]~h^-1)), print = F)

# inset_plot<-ggdraw(MRcorrel1.large + theme_half_open( font_size =20)) + theme(legend.position="none") +
#   draw_plot(MRcorrel1.small + theme(legend.position="none") , 0.57, .58, .42, .4)

cowplot:::plot_grid(MRcorrel1.large, MRcorrel1.small,
                         align = "hv",
                         axis = "l",
                         nrow = 2,
                         ncol = 1
                         # labels = "AUTO",
                         # label_x = c(0.33, 0.2),
                         # label_y = c(0.97, 0.97, 0.97, 0.97),
                         # label_size = 17,
                    ) %>%
ggsave(filename = paste("./Figures/AMR_RMR_correl", Sys.Date(),".png"), width = 5, height = 7)

AScorrel1.small<-ggplot(data.amr[data.amr$BW_g<100,],
                        aes(x=AMR, y=AS, group=test_category3, color = BW_g,
                            shape=test_category3, fill=BW_g)) +
  geom_point(alpha=0.6, pch=21,  show.legend = FALSE, size=2)+
  # scale_fill_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  # scale_color_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_shape_manual(values= c(22,23))+
  facet_grid(.~test_category3)+
  geom_abline(slope = 1, intercept = 0, lty = "dashed")
ggformat(AScorrel1.small, title = "",x_title = expression(MMR~(mg~O[2]~h^-1)), y_title = expression(AS~(mg~O[2]~h^-1)), print = F)

AScorrel1.large<-ggplot(data.amr,
                        aes(y=AS, x=AMR, group=test_category, fill=tempTest, shape=test_category)) +
  geom_point(alpha=0.6,shape=21,  show.legend = FALSE, size=3)+
  geom_point(data.amr[data.amr$species=="Somniosus microcephalus",],  mapping=aes(x=AS, y=AMR, group=test_category, fill=tempTest, shape=test_category), shape=23,  show.legend = FALSE, size=3, stroke=2)+
  # scale_color_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  # scale_fill_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_shape_manual(values= c(22,23,24))+
  facet_grid(.~test_category3)+
  scale_x_continuous(breaks = c(0, 1000, 3000, 5000))+
  geom_abline(slope = 1, intercept = 0, lty = "dashed")
ggformat(AScorrel1.large, title = "", x_title = expression(MMR~(mg~O[2]~h^-1)), y_title = expression(AS~(mg~O[2]~h^-1)), print = F)

cowplot:::plot_grid(AScorrel1.large, AScorrel1.small,
                    align = "hv",
                    axis = "l",
                    nrow = 2,
                    ncol = 1) %>%
  # fig_part4.1
  ggsave(filename = "./Figures/AMR_AS_correl_march232022.png", width = 5, height = 7)

AScorrel1.small.rmr<-ggplot(data.amr[data.amr$BW_g<100,],
                            aes(x=RMR, y=AS, group=test_category3,
                                color = BW_g, shape=test_category3,
                                fill=BW_g)) +
  geom_point(alpha=0.6, pch=21,  show.legend = FALSE, size=2)+
  # scale_fill_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  # scale_color_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_shape_manual(values= c(22,23))+
  facet_grid(.~test_category3)+
  geom_abline(slope = 1, intercept = 0, lty = "dashed")
ggformat(AScorrel1.small.rmr, title = "",x_title = expression(RMR~(mg~O[2]~h^-1)), y_title = expression(AS~(mg~O[2]~h^-1)), print = F)

AScorrel1.large.rmr<-ggplot(data.amr,
                            aes(y=AS, x=RMR, group=test_category,
                                fill=BW_g, shape=test_category)) +
  geom_point(alpha=0.6,shape=21,  show.legend = FALSE, size=3)+
  # scale_color_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  # scale_fill_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_shape_manual(values= c(22,23,24))+
  facet_grid(.~test_category3)+
  scale_x_continuous(breaks = c(0, 1000, 3000, 5000))+
  geom_abline(slope = 1, intercept = 0, lty = "dashed")
ggformat(AScorrel1.large.rmr, title = "", x_title = expression(RMR~(mg~O[2]~h^-1)), y_title = expression(AS~(mg~O[2]~h^-1)), print = F)

cowplot:::plot_grid(AScorrel1.large, AScorrel1.small,
                    align = "hv",
                    axis = "l",
                    nrow = 2,
                    ncol = 1) %>%
  # fig_part4.1
  ggsave(filename = "./Figures/AMR_AS_correl_march232022.png", width = 5, height = 7)

