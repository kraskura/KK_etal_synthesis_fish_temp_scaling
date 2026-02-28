# Author: Krista Kraskura
# Description: 
#   Supplemental material:
#     - histograms 
#     - correlations between MR metrics
# ********************************************
# ********************************************

# import data, source local scripts ---------
# source("./R/phylo_mixed_model.R") # need for scaling slopes 
# source("./R/nonPhylo_mixed_models.R") # for comparison supplemental 
set.seed(51423)
# MMR and AMR used interchangeably throughout 
library(here)
source(here("R/get_data_temp.R"))
source(here("R/setup.R")) # data set, libraries, colors, etc. 

data.amr$chambSizeRatio<-data.amr$chamber_vol_L/(data.amr$BW_g/1000)
data.rmr$chambSizeRatio<-data.rmr$chamber_vol_L/(data.rmr$BW_g/1000)

# Supplemental -----------
## Experimental checks -------
### type of MMR ------
ggplot(data.amr.test, aes(x = lnBWg, y = lnAMR, color = MMR_method, alpha = tempTest))+
  geom_point(size = 1)+
  theme_classic()+
  scale_color_discrete()+
  ggtitle("AMR warm")

ggplot(data.amrER, aes(x = lnBWg, y = lnAMR, color = MMR_method, alpha = tempTest))+
  geom_point(size = 1)+
  theme_classic()+
  scale_color_discrete()+
  ggtitle("AMR ecologically relev")

### the size of chamber  -------
ggplot(data.amr, aes(y = chambSizeRatio, x = BW_g/1000,
                     color = MMR_method))+
  geom_point(size = 1)+
  theme_classic()+
  ggtitle("AMR all")

ggplot(data.rmr, aes(y = chambSizeRatio, x = BW_g/1000))+
  geom_point(size = 1)+
  theme_classic()+
  ggtitle("RMR all")

ggplot(data.amr, aes(x = chambSizeRatio,
                     fill = MMR_method))+
  geom_histogram(bins = 100)+
  theme_classic()+
  ggtitle("AMR all")

ggplot(data.rmr, aes(x = chambSizeRatio))+
  geom_histogram(bins = 100)+
  theme_classic()+
  ggtitle("RMR all")

ggplot(data.amr, aes(x = lnBWg, y = lnAMR, color = chambSizeRatio))+
  geom_point(size = 1, alpha = 0.4)+
  theme_classic()+
  scale_color_viridis_c()+
  ggtitle("AMR all")

ggplot(data.rmr, aes(x = lnBWg, y = lnRMR, color = chambSizeRatio))+
  geom_point(size = 1, alpha = 0.4)+
  theme_classic()+
  scale_color_viridis_c()+
  ggtitle("RMR all")

## Histograms of fish at different size suppl ----------
plot_hist_amr<-ggplot(data.amr, aes(x=BW_g,
                                    fill = test_category, 
                                    color = test_category)) +
  geom_vline(xintercept = 1000, color="grey", lty="dashed")+
  geom_histogram(show.legend = FALSE, bins=30, alpha = 1 ,
                 color = "white", linewidth=0.1)+
  scale_fill_manual(values = c(cols.amr[2], "black", cols.amr[3]))+
  ylim(0,110)+
  scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_amr, x_title = "Body mass (g)", y_title = "Frequency", print = FALSE)
plot_hist_amr <- plot_hist_amr + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                             axis.text.x = element_text(size=12,  family="Helvetica"),
                                             axis.title.y = element_text(size=12,  family="Helvetica"),
                                             axis.text.y = element_text(size=12,  family="Helvetica"))

plot_hist_rmr<-ggplot(data.rmr, aes(x=BW_g ,
                                    fill = test_category, 
                                    color = test_category)) +
  geom_vline(xintercept = 1000, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=30, alpha = 1 ,
                 color = "grey20",
                 linewidth = 0.1)+
  scale_fill_manual(values = c(cols.rmr[1], "black", cols.rmr[3]))+
  ylim(0,120)+
  scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_rmr, x_title = "Body mass (g)", y_title = "Frequency", print = FALSE)
plot_hist_rmr <- plot_hist_rmr + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                       axis.text.x = element_text(size=12,  family="Helvetica"),
                                       axis.title.y = element_text(size=12,  family="Helvetica"),
                                       axis.text.y = element_text(size=12,  family="Helvetica"))

plot_hist_fas<-ggplot(data.fas, aes(x=BW_g,
                                    fill = test_category, 
                                    color = test_category)) +
  geom_vline(xintercept = 1000, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=30, alpha = 1,
                 color = "white",
                 linewidth = 0.1)+
  scale_fill_manual(values = c(cols.fas[1], "black", cols.fas[3]))+
  ylim(0,100)+
  scale_x_continuous(limits=c(0, 3500))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_fas, x_title = "Body mass (g)", y_title = "Frequency", print = FALSE)
plot_hist_fas <- plot_hist_fas + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                       axis.text.x = element_text(size=12,  family="Helvetica"),
                                       axis.title.y = element_text(size=12,  family="Helvetica"),
                                       axis.text.y = element_text(size=12,  family="Helvetica"))
plot_hist_fas

## Histograms of fish at different temperatures suppl (deg C) ----------
plot_hist_amr_T<-ggplot(data.amr, aes(x=tempTest,
                                    fill = test_category, 
                                    color = test_category)) +
  geom_vline(xintercept = 15, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=30, alpha = 1,
                 color = "white",
                 linewidth = 0.1)+
  scale_fill_manual(values = c(cols.amr[2], "black", cols.amr[3]))+
  # ylim(0,300)+
  # scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_amr_T, x_title = expression(Temperature~degree*C), y_title = "Frequency", print = FALSE)
plot_hist_amr_T <- plot_hist_amr_T + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                       axis.text.x = element_text(size=12,  family="Helvetica"),
                                       axis.title.y = element_text(size=12,  family="Helvetica"),
                                       axis.text.y = element_text(size=12,  family="Helvetica"))

plot_hist_rmr_T<-ggplot(data.rmr, aes(x=tempTest,
                                    fill = test_category, 
                                    color = test_category)) +
  geom_vline(xintercept = 15, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=30, color = "grey30",
                 linewidth = 0.1)+
  scale_fill_manual(values = c(cols.rmr[1], "black", cols.rmr[3]))+
  # ylim(0,300)+
  # scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_rmr_T, x_title =  expression(Temperature~degree*C), y_title = "Frequency", print = FALSE)
plot_hist_rmr_T <- plot_hist_rmr_T + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                       axis.text.x = element_text(size=12,  family="Helvetica"),
                                       axis.title.y = element_text(size=12,  family="Helvetica"),
                                       axis.text.y = element_text(size=12,  family="Helvetica"))

plot_hist_fas_T<-ggplot(data.fas, aes(x=tempTest,
                                    fill = test_category, 
                                    color = test_category)) +
  geom_vline(xintercept = 15, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=30,color = "white",
                 linewidth = 0.1)+
  scale_fill_manual(values = c(cols.fas[1], "black", cols.fas[3]))+
  # ylim(0,300)+
  # scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_fas_T, x_title =  expression(Temperature~degree*C), y_title = "Frequency", print = FALSE)
plot_hist_fas_T <- plot_hist_fas_T + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                       axis.text.x = element_text(size=12,  family="Helvetica"),
                                       axis.title.y = element_text(size=12,  family="Helvetica"),
                                       axis.text.y = element_text(size=12,  family="Helvetica"))


hist.plot<-cowplot:::plot_grid(plot_hist_amr, plot_hist_amr_T, 
                               plot_hist_rmr, plot_hist_rmr_T,
                               plot_hist_fas, plot_hist_fas_T,
                               align = "hv",
                               axis = "l",
                               nrow = 3,
                               ncol = 2, label = "AUTO") 
ggsave(filename = paste("./Figures/Suppl_Size_hist.png", sep=""),
          hist.plot, width = 13, height = 9, units = "in")

## Plots: Correlations metabolic rates ---------
MRcorrel1.small<-ggplot(data.amr[data.amr$BW_g<1000,], aes(x=RMR, y=AMR, group=test_category3, color = tempTest, shape=test_category3, fill=tempTest)) +
  geom_point(alpha=1, pch=1,  show.legend = TRUE, size=1)+
  scale_fill_viridis_c(option = "C")+
  scale_color_viridis_c(option = "C")+
  scale_shape_manual(values= c(22,23))+
  geom_abline(slope = 1, intercept = 0, lty = "dashed")+
  facet_grid(.~test_category3)
  # ylim(0, 75)+
  # xlim(0, 60)
ggformat(MRcorrel1.small, title = "",y_title = expression(MMR~(mg~O[2]~h^-1)), x_title = expression(RMR~(mg~O[2]~h^-1)), print = F, size_text = 12)

MRcorrel1.large<-ggplot(data.amr[data.amr$BW_g>=1000,], aes(x=RMR, y=AMR, group=test_category, fill=tempTest, color = tempTest, shape=test_category)) +
  geom_point(alpha=1,shape=1,  show.legend = TRUE, size=1)+
  geom_point(data.amr[data.amr$species=="Somniosus microcephalus",],
             mapping=aes(x=RMR, y=AMR, group=test_category,
                         color = tempTest, 
                         fill=tempTest,
                         shape=test_category), shape=23,
             show.legend = FALSE, size=2, stroke=1)+
  geom_abline(slope = 1, intercept = 0, lty = "dashed")+
  scale_fill_viridis_c(option = "C")+
  scale_color_viridis_c(option = "C")+
  scale_shape_manual(values= c(22,23,24))+
  facet_grid(.~test_category3)
ggformat(MRcorrel1.large, title = "", y_title = expression(MMR~(mg~O[2]~h^-1)), x_title = expression(RMR~(mg~O[2]~h^-1)), print = F, size_text = 12)

# inset_plot<-ggdraw(MRcorrel1.large + theme_half_open( font_size =20)) + theme(legend.position="none") +
#   draw_plot(MRcorrel1.small + theme(legend.position="none") , 0.57, .58, .42, .4)

cowplot:::plot_grid(MRcorrel1.small,MRcorrel1.large,
                         align = "hv",
                         axis = "l",
                         nrow = 2,
                         ncol = 1, labels = c("A (< 1000 g)", "B (>= 1000 g)")) %>%
ggsave(filename = paste("./Figures/suppl_MMRvsRMR.png"), width = 7, height = 7)

AScorrel1.small<-ggplot(data.amr[data.amr$BW_g<1000,],
                        aes(x=AMR, y=AS, group=test_category3,
                            color = tempTest,
                            shape=test_category3, fill=tempTest)) +
  geom_point(alpha=1, pch=1,  show.legend = TRUE, size=1)+
  scale_fill_viridis_c(option = "C")+
  scale_color_viridis_c(option = "C")+
  scale_shape_manual(values= c(22,23))+
  facet_grid(.~test_category3)+
  geom_abline(slope = 1, intercept = 0, lty = "dashed")
ggformat(AScorrel1.small, title = "",x_title = expression(MMR~(mg~O[2]~h^-1)), y_title = expression(AS~(mg~O[2]~h^-1)), print = F, size_text = 12)

AScorrel1.large<-ggplot(data.amr,
                        aes(y=AS, x=AMR, group=test_category,
                            fill=tempTest, shape=test_category,
                            color = tempTest)) +
  geom_point(alpha=1,shape=1,  show.legend = TRUE, size=1)+
  geom_point(data.amr[data.amr$species=="Somniosus microcephalus",],  mapping=aes(x=AS, y=AMR, group=test_category, fill=tempTest, shape=test_category), shape=23,  show.legend = FALSE, size=2, stroke=1)+
  scale_fill_viridis_c(option = "C")+
  scale_color_viridis_c(option = "C")+
  scale_shape_manual(values= c(22,23,24))+
  facet_grid(.~test_category3)+
  scale_x_continuous(breaks = c(0, 1000, 3000, 5000))+
  geom_abline(slope = 1, intercept = 0, lty = "dashed")
ggformat(AScorrel1.large, title = "", x_title = expression(MMR~(mg~O[2]~h^-1)), y_title = expression(AS~(mg~O[2]~h^-1)), print = F, size_text = 12)

cowplot:::plot_grid(AScorrel1.small,AScorrel1.large,
                    align = "hv",
                    axis = "l",
                    nrow = 2,
                    ncol = 1, labels = c("A (< 1000 g)", "B (>= 1000 g)")) %>%
  # fig_part4.1
  ggsave(filename = "./Figures/suppl_MMRvsAS.png", width = 7, height = 7)

AScorrel1.small.rmr<-ggplot(data.amr[data.amr$BW_g<1000,],
                            aes(x=RMR, y=AS, group=test_category3,
                                color = tempTest, shape=test_category3,
                                fill=tempTest)) +
  geom_point(alpha=1, pch=1,  show.legend = TRUE, size=1)+
  scale_fill_viridis_c(option = "C")+
  scale_color_viridis_c(option = "C")+
  scale_shape_manual(values= c(22,23))+
  facet_grid(.~test_category3)+
  geom_abline(slope = 1, intercept = 0, lty = "dashed")
ggformat(AScorrel1.small.rmr, title = "",x_title = expression(RMR~(mg~O[2]~h^-1)), y_title = expression(AS~(mg~O[2]~h^-1)), print = F, size_text = 12)

AScorrel1.large.rmr<-ggplot(data.amr,
                            aes(y=AS, x=RMR, group=test_category,
                                fill=tempTest, shape=test_category, 
                                color = tempTest)) +
  geom_point(alpha=1,shape=1,  show.legend = TRUE, size=1)+
  geom_point(data.as[data.as$species=="Somniosus microcephalus",],  mapping=aes(x=AS, y=RMR, group=test_category, fill=tempTest, shape=test_category), shape=23,  show.legend = FALSE, size=2, stroke=1)+
  scale_fill_viridis_c(option = "C")+
  scale_color_viridis_c(option = "C")+
  scale_shape_manual(values= c(22,23,24))+
  facet_grid(.~test_category3)+
  scale_x_continuous(breaks = c(0, 1000, 3000, 5000))+
  geom_abline(slope = 1, intercept = 0, lty = "dashed")
ggformat(AScorrel1.large.rmr, title = "", x_title = expression(RMR~(mg~O[2]~h^-1)), y_title = expression(AS~(mg~O[2]~h^-1)), print = F, size_text = 12)

cowplot:::plot_grid(AScorrel1.small.rmr,AScorrel1.large.rmr,
                    align = "hv",
                    axis = "l",
                    nrow = 2,
                    ncol = 1, labels = c("A (< 1000 g)", "B (>= 1000 g)")) %>%
  ggsave(filename = "./Figures/suppl_RMRvsAS.png", width = 7, height = 7)


# fas in pelagic fish 
# summary(data.fasER[data.fasER$DemersPelag == "pelagic","FAS"])
# summary(data.fas.test[data.fas.test$DemersPelag == "pelagic","FAS"])


