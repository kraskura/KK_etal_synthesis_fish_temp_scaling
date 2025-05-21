# last run: april 27 2025

# import data, source local scripts ---------
# source("./R/phylo_mixed_model.R") # needed libraries called within 
# source("./R/nonPhylo_mixed_models.R") # for comparison supplemental 
set.seed(51423)
# MMR and AMR used interchangeably throughout 

library(ggplot2)
library(ggpubr)
library(cowplot)
library(forcats)
library(here)
library(weathermetrics)
library(ggformat2)
  
source("./R/get_data_temp.R")
source("./R/colors_themes.R")

# Data for models ------
data.list<-get_data_temp(data.amr = here("Data", "Fish_AMR_temp_dataset_mar2022.csv"),
                                 data.rmr = here("Data","Fish_RMR_temp_dataset_mar2022.csv"),
                                 ecology.data = here("Data", "Kraskura_species_ecologies_mar2022.csv"),
                                 onlyTop.above = TRUE,
                                 save.FishBase.species.data = F,
                                 calc_mass_specific = FALSE)

# data.list.NO.SALMON<-get_data_temp(data.amr = here("Data", "Fish_AMR_temp_dataset_mar2022.csv"),
#                                  data.rmr = here("Data","Fish_RMR_temp_dataset_mar2022.csv"),
#                                  ecology.data = here("Data", "Kraskura_species_ecologies_mar2022.csv"),
#                                  onlyTop.above = TRUE, save.FishBase.species.data = F,
#                                  calc_mass_specific = FALSE)

data.amrAC<-data.frame(data.list[1])
data.rmrAC<-data.frame(data.list[2])
data.amrAM<-data.frame(data.list[3])
data.rmrAM<-data.frame(data.list[4])
data.amrER<-data.frame(data.list[5])
data.rmrER<-data.frame(data.list[6])

data.asAC<-data.frame(data.list[7])
data.fasAC<-data.frame(data.list[8])
data.asAM<-data.frame(data.list[9])
data.fasAM<-data.frame(data.list[10])
data.asER<-data.frame(data.list[11])
data.fasER<-data.frame(data.list[12])

data.amr<-data.frame(data.list[13])
data.rmr<-data.frame(data.list[14])
dataMR<-data.frame(data.list[15])
data.as<-data.frame(data.list[16])
data.fas<-data.frame(data.list[17])

# the warm ones
data.amr.test<-rbind(data.amrAC, data.amrAM)
data.rmr.test<-rbind(data.rmrAC, data.rmrAM)
data.fas.test<-rbind(data.fasAC, data.fasAM)
data.as.test<-rbind(data.asAC, data.asAM)
data.fas.test<-data.fas.test[c(!is.na(data.fas.test$FAS) & is.finite(data.fas.test$FAS)) , ]
data.as.test<-data.as.test[c(!is.na(data.as.test$lnAS) & is.finite(data.as.test$lnAS)) , ]

# Get model data set specific phylogentics model matrixes -------
data.rmrER<-droplevels(data.rmrER)
data.rmr.test<-droplevels(data.rmr.test)
data.amrER<-droplevels(data.amrER)
data.amr.test<-droplevels(data.amr.test)
data.asER<-droplevels(data.asER)
data.as.test<-droplevels(data.as.test)
data.fasER<-droplevels(data.fasER)
data.fas.test<-droplevels(data.fas.test)

# test categories
data.fas$test_category3 <- "warm"
data.fas[data.fas$test_category == "ecol_relev", "test_category3"] <- "optimal"
data.as$test_category3 <- "warm"
data.as[data.as$test_category == "ecol_relev", "test_category3"] <- "optimal"
data.amr$test_category3 <- "warm"
data.amr[data.amr$test_category == "ecol_relev", "test_category3"] <- "optimal"
data.rmr$test_category3 <- "warm"
data.rmr[data.rmr$test_category == "ecol_relev", "test_category3"] <- "optimal"


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
  facet_wrap(.~test_category)
ggformat(plot_hist_amr, x_title = "Body mass (g)", y_title = "Frequency", print = FALSE)
plot_hist_amr <- plot_hist_amr + theme(axis.title.x = element_text(size=12,  family="Arial"),
                                             axis.text.x = element_text(size=12,  family="Arial"),
                                             axis.title.y = element_text(size=12,  family="Arial"),
                                             axis.text.y = element_text(size=12,  family="Arial"))

plot_hist_rmr<-ggplot(data.rmr, aes(x=BW_g)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  # geom_vline(xintercept = 300, color="grey", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = 1000, color="grey", lty="dashed", lwd=0.5)+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.rmr[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,300)+
  scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category)
ggformat(plot_hist_rmr, x_title = "Body mass (g)", y_title = "Frequency", print = FALSE)
plot_hist_rmr <- plot_hist_rmr + theme(axis.title.x = element_text(size=12,  family="Arial"),
                                       axis.text.x = element_text(size=12,  family="Arial"),
                                       axis.title.y = element_text(size=12,  family="Arial"),
                                       axis.text.y = element_text(size=12,  family="Arial"))

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
plot_hist_fas <- plot_hist_fas + theme(axis.title.x = element_text(size=12,  family="Arial"),
                                       axis.text.x = element_text(size=12,  family="Arial"),
                                       axis.title.y = element_text(size=12,  family="Arial"),
                                       axis.text.y = element_text(size=12,  family="Arial"))


## Histograms of fish at different temperatures suppl (deg C) ----------

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
plot_hist_amr_T <- plot_hist_amr_T + theme(axis.title.x = element_text(size=12,  family="Arial"),
                                       axis.text.x = element_text(size=12,  family="Arial"),
                                       axis.title.y = element_text(size=12,  family="Arial"),
                                       axis.text.y = element_text(size=12,  family="Arial"))

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
plot_hist_rmr_T <- plot_hist_rmr_T + theme(axis.title.x = element_text(size=12,  family="Arial"),
                                       axis.text.x = element_text(size=12,  family="Arial"),
                                       axis.title.y = element_text(size=12,  family="Arial"),
                                       axis.text.y = element_text(size=12,  family="Arial"))

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
plot_hist_fas_T <- plot_hist_fas_T + theme(axis.title.x = element_text(size=12,  family="Arial"),
                                       axis.text.x = element_text(size=12,  family="Arial"),
                                       axis.title.y = element_text(size=12,  family="Arial"),
                                       axis.text.y = element_text(size=12,  family="Arial"))


hist.plot<-cowplot:::plot_grid(plot_hist_amr, plot_hist_amr_T, 
                               plot_hist_rmr, plot_hist_rmr_T,
                               plot_hist_fas, plot_hist_fas_T,
                               align = "hv",
                               axis = "l",
                               nrow = 3,
                               ncol = 2) 
ggsave(filename = paste("./Figures/Suppl_Size_hist.png", sep=""),
          hist.plot, width = 13, height = 9, units = "in")


## Histograms of fish at different temperatures suppl (Arrhenius Temp) ----------

plot_hist_amr_TK<-ggplot(data.amr, aes(x=tempTestK1000)) +
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.amr[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  # ylim(0,300)+
  # scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_amr_TK, x_title = expression(Temperature~(1000~K^-1)), y_title = "Frequency", print = FALSE)
plot_hist_amr_TK <- plot_hist_amr_TK + theme(axis.title.x = element_text(size=12,  family="Arial"),
                                       axis.text.x = element_text(size=12,  family="Arial"),
                                       axis.title.y = element_text(size=12,  family="Arial"),
                                       axis.text.y = element_text(size=12,  family="Arial"))

plot_hist_rmr_TK<-ggplot(data.rmr, aes(x=tempTestK1000)) +
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.rmr[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  # ylim(0,300)+
  # scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_rmr_TK, x_title =  expression(Temperature~(1000~K^-1)), y_title = "Frequency", print = FALSE)
plot_hist_rmr_TK <- plot_hist_rmr_TK + theme(axis.title.x = element_text(size=12,  family="Arial"),
                                       axis.text.x = element_text(size=12,  family="Arial"),
                                       axis.title.y = element_text(size=12,  family="Arial"),
                                       axis.text.y = element_text(size=12,  family="Arial"))

plot_hist_fas_TK<-ggplot(data.fas, aes(x=tempTestK1000)) +
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.fas[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  # ylim(0,300)+
  # scale_x_continuous(limits=c(0, 10000))+
  facet_wrap(.~test_category3)
ggformat(plot_hist_fas_TK, x_title =  expression(Temperature~(1000~K^-1)), y_title = "Frequency", print = FALSE)
plot_hist_fas_TK <- plot_hist_fas_TK + theme(axis.title.x = element_text(size=12,  family="Arial"),
                                       axis.text.x = element_text(size=12,  family="Arial"),
                                       axis.title.y = element_text(size=12,  family="Arial"),
                                       axis.text.y = element_text(size=12,  family="Arial"))

hist.plot<-cowplot:::plot_grid(plot_hist_amr, plot_hist_amr_TK, 
                               plot_hist_rmr, plot_hist_rmr_TK,
                               plot_hist_fas, plot_hist_fas_TK,
                               align = "hv",
                               axis = "l",
                               nrow = 3,
                               ncol = 2) 


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
ggsave(filename = paste("./Figures/suppl_x1.png"), width = 7, height = 7)

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
  ggsave(filename = "./Figures/suppl_MMRx2.png", width = 7, height = 7)

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

cowplot:::plot_grid(AScorrel1.large.rmr, AScorrel1.small.rmr,
                    align = "hv",
                    axis = "l",
                    nrow = 2,
                    ncol = 1) %>%
  # fig_part4.1
  ggsave(filename = "./Figures/suppl_RMRx3.png", width = 7, height = 7)


# fas in pelagic fish 
summary(data.fasER[data.fasER$DemersPelag == "pelagic","FAS"])
summary(data.fas.test[data.fas.test$DemersPelag == "pelagic","FAS"])


