
set.seed(51423)

library(ggplot2)
library(ggpubr)
library(here)
library(cowplot)
library(forcats)
# MMR and AMR used interchangeably throughout 

cols.as<-c("#265F73", "#007E66", "#00C5A3")
cols.fas<-c("#395200", "#89A000", "yellow")
cols.rmr<-c("#C70039", "#FF6D7C", "#FFA3AC")
cols.amr<-c("#00749F","#00A8D6", "#9CE9FF")
cols<-c("#00749F","#C70039","#00A8D6","#FF6D7C", "#9CE9FF","#FFA3AC", "#00C5A3", "#265F73")# AMR -rmr- AMR dark - rmr dark - light - as - fas

# FUNCTION BELOW ---------
# Sys.Date()
date.imports<-substr(list.files("./Data_exports/Phylo/")[1], start = 17, stop = 26)
getModelCIs<-FALSE

# ^ this is computationally heavy task, to get saved data: "NO" or 2
# to get compute now: "YES" or 1 
source("./R/phylo_mixed model.R")
source("./R/get_data_for_figures.R")

# General scaling plots ------
AMRmodel_plot1<-ggplot(data=data.amrER, aes(x=lnBWg, y=lnAMR)) +
  geom_point(alpha=0.9,  size=2, pch=19, color="grey70")+
  geom_line(data=data.plotAMRint_ER,
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000_inC, color = tempTestK1000_inC),
            linewidth=0.3, lty=1,alpha=0.8, show.legend=FALSE) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_point(data=data.amr.test, aes(x=lnBWg, y=lnAMR, fill=tempTest),
             alpha=0.9,  size=2, pch=21, show.legend = FALSE)+
  geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTestK1000,2)==3.39,],
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="#002E53", linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.5, y = 11.5, label = bquote(Optimal:~italic(b)[MMR] == change~with~degree*C),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -5.5, y = 9.8, label = bquote(Warm:~italic(b)[MMR] == .(MMR_slope_w)),size=5, hjust=0, family="Arial", color = cols.amr[1])+
  annotate("text",  x = 3.0, y = -2, label = bquote(~0*degree*C:~italic(b)[MMR] == .(round(AMR.slopes$lnBWg.trend[1],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.0, y = -3, label = bquote(10*degree*C:~italic(b)[MMR] == .(round(AMR.slopes$lnBWg.trend[2],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.0, y = -4, label = bquote(20*degree*C:~italic(b)[MMR] == .(round(AMR.slopes$lnBWg.trend[3],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.0, y = -5, label = bquote(30*degree*C:~italic(b)[MMR] == .(round(AMR.slopes$lnBWg.trend[4],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.0, y = -6, label = bquote(40*degree*C:~italic(b)[MMR] == .(round(AMR.slopes$lnBWg.trend[5],3))),size=4, hjust=0, family="Arial", color = "black")+
  scale_fill_gradient( low = cols.amr[3], high = cols.amr[1])+
  scale_color_gradient( low = "grey80", high = "grey0")+
  ylim(x = -6.5, 12)+
  xlim(x = -6.5, 12)+
  annotate("text", label = paste("n = ", nrow(data.amrER), sep=""),   x = -5.5, y = 8.5, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.amr.test), sep=""),  x = -5.5, y = 7.6, size=3, hjust=0, family="Arial", color = cols.amr[1])
ggformat(AMRmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), print = F)

RMRmodel_plot1<-ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,2)==3.39,], mapping=aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", linewidth=1, lty=1, show.legend=FALSE) +
  geom_point(data=data.rmr.test, aes(x=lnBWg, y=lnRMR, fill=tempTest), alpha=0.9,  size=2, pch=21, show.legend = FALSE)+
  geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTestK1000,2)==3.39,],
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="#7C0003", linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.5, y = 11.5, label = bquote(Optimal:~italic(b)[RMR] == .(RMR_slope)),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -5.5, y = 9.8, label = bquote(Warm:~italic(b)[RMR] == .(RMR_slope_w)),size=5, hjust=0, family="Arial", color = cols.rmr[1])+
  scale_fill_gradient( low = cols.rmr[3], high = cols.rmr[1])+
  scale_color_gradient( low = cols.rmr[3], high = cols.rmr[1])+
  ylim(x = -6.5, 12)+
  xlim(x = -6.5, 12)+
  annotate("text", label = paste("n = ", nrow(data.rmrER), sep=""),   x = -5.5, y = 8.5, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.rmr.test), sep=""),  x = -5.5, y = 7.6, size=3, hjust=0, family="Arial", color = cols.rmr[1])
ggformat(RMRmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = FALSE)

# FAS! 
FASmodel_plot1<-ggplot(data=data.fasER, aes(x=lnBWg, y=log(FAS))) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_line(data=data.plotFAS_ER, aes(y = model_predFE, x=lnBWg,  group=tempTest), color="black", linewidth=1, lty=1, show.legend=FALSE) +
  geom_point(data=data.fas.test, aes(x=lnBWg, y=log(FAS), fill=tempTest), alpha=0.9,  size= 2, pch=21, show.legend=FALSE)+
  geom_line(data=data.plotFAS_warm, aes(y = model_predFE, x=lnBWg,  group=tempTest), color="#475500", linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.5, y = 3.9, label = bquote(Optimal:~italic(b)[FAS] == .(FAS_slope)),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -5.5, y = 3.5, label = bquote(Warm:~italic(b)[FAS] == .(FAS_slope_w)),size=5, hjust=0, family="Arial", color = "#475500")+
  annotate("text", label = paste("n = ", nrow(data.fasER), sep=""),   x = -5.5, y = 3.1, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.fas.test), sep=""),  x = -5.5, y = 2.85, size=3, hjust=0, family="Arial", color = "#475500")+
  scale_fill_gradient(high = "yellow", low = "black")+
  ylim(0,4)
ggformat(FASmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*FAS), print = FALSE)

# AS 
ASmodel_plot1<-ggplot(data=data.asER, aes(x=lnBWg, y=lnAS)) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_line(data=data.plotAS_ER[round(data.plotAS_ER$tempTestK1000,2)==3.39,], mapping=aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", linewidth=1, lty=1, show.legend=FALSE) +
  geom_point(data=data.as.test, aes(x=lnBWg, y=lnAS, fill=tempTest), alpha=0.9,  size=2, pch=21, show.legend = FALSE)+
  geom_line(data=data.plotAS_warm[round(data.plotAS_warm$tempTestK1000,1)==round(C20inTempTestK1000,1),],
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="#00785B", linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.5, y = 11.5, label = bquote(Optimal:~italic(b)[AS] == .(AS_slope)),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -5.5, y = 9.8, label = bquote(Warm:~italic(b)[AS] == .(AS_slope_w)),size=5, hjust=0, family="Arial", color = cols.as[1])+
  scale_fill_gradient( low = cols.as[3], high = cols.as[1])+
  scale_color_gradient( low = cols.as[3], high = cols.as[1])+
  ylim(x = -6.5, 12)+
  xlim(x = -6.5, 12)+
  annotate("text", label = paste("n = ", nrow(data.asER), sep=""),   x = -5.5, y = 8.5, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.as.test), sep=""),  x = -5.5, y = 7.6, size=3, hjust=0, family="Arial", color = cols.as[1])
ggformat(ASmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1)), print = F)


scaling<-cowplot:::plot_grid(AMRmodel_plot1, RMRmodel_plot1,
                             ASmodel_plot1, FASmodel_plot1,
                              align = "hv",
                              axis = "l",
                              nrow = 2,
                              ncol = 2,
                              labels = "AUTO",
                              label_size = 17)
ggsave(filename = paste("./Figures/Fig1_Scaling_aLL_", Sys.Date(), ".png", sep=""),
       plot=scaling, width = 8.5, height = 8.5, units = "in")


# Both Activation energies together: -------
# model_predFE << is in lnMR units 
# MMR
AMRmodel_plot3.E_ALL<-ggplot(data.amrER[data.amrER$lnBWg==0,]) +
  geom_point(data.amrER, mapping=aes(x=tempTestK1000, y=log(mass_specamr), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21, show.legend = FALSE)+
  geom_point(data.amr.test, mapping=aes(x=tempTestK1000, y=log(mass_specamr), fill= tempTest, size = lnBWg), color="black",  alpha=1, pch=21,  show.legend = FALSE)+
  scale_fill_gradient( low = cols.amr[3], high = cols.amr[1])+
  scale_color_gradient( low = "grey80", high = "grey0")+
  ylim(-5.7,4)+
  annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[MMR] == change~with~mass),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[MMR] == .(round(MMR_E_W_eV,3))),size=5, hjust=0, family="Arial", color = cols.amr[1])+
  annotate("text",  x = 3.27, y = -3.2, label = bquote(1*g~italic(E)== .(round(MMR_E_ER$MMR_E_ER_eV[1],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = -4, label = bquote(10*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[2],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = -4.8, label = bquote(100*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[3],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = -5.6, label = bquote(1000*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[4],3))),size=4, hjust=0, family="Arial", color = "black")+
  scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
  geom_line(data = data.plotAMR_warm[which(round(data.plotAMR_warm$lnBWg, 1) == round(log(1.4), 1)),],
            aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg) ,color=cols.amr[2], linewidth=1, lty=1, show.legend=FALSE)+
  geom_line(data = data.plotAMRint_ER[exp(data.plotAMRint_ER$lnBWg) <= 1000, ], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg, color = lnBWg), linewidth=0.5, lty=1, show.legend=FALSE)
ggformat(AMRmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/K)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1~g^-1)), print = F)

RMRmodel_plot3.E_ALL<-ggplot(data.rmrER[data.rmrER$lnBWg==0,]) +
  geom_point(data.rmrER, mapping=aes(x=tempTestK1000, y=log(mass_specrmr), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21)+
  geom_point(data.rmr.test, mapping=aes(x=tempTestK1000, y=log(mass_specrmr), fill= tempTest, size = lnBWg), color = "black", alpha=1, pch=21, show.legend = FALSE)+
  scale_fill_gradient( low = cols.rmr[3], high = cols.rmr[1])+
  annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[RMR] == .(round(RMR_E_ER_eV,3))),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[RMR] == .(round(RMR_E_W_eV,3))),size=5, hjust=0, family="Arial", color = cols.rmr[1])+
  ylim(-5.7,4)+
  scale_size_continuous(name="Mass (g)",
                      breaks=c(log(0.1), log(10), log(1000)),
                      labels=c("0.1", "10", "1000"))+
  scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
  geom_line(data = data.plotRMR_ER[which(round(data.plotRMR_ER$lnBWg, 1) == round(log(1.35), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg ) ,color="black", linewidth=1, lty=1, show.legend=FALSE)+
  geom_line(data = data.plotRMR_warm[which(round(data.plotRMR_warm$lnBWg, 1) == round(log(1.35), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg), color = cols.rmr[2], linewidth=1, lty=1, show.legend=FALSE)
ggformat(RMRmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/K)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1~g^-1)), print = FALSE)
RMRmodel_plot3.E_ALL <- RMRmodel_plot3.E_ALL + theme(legend.position = c(0.87, 0.82))

ASmodel_plot3.E_ALL<-ggplot(data.asER[data.asER$lnBWg==0,]) +
  geom_point(data.asER, mapping=aes(x=tempTestK1000, y=log(mass_specas), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21, show.legend = FALSE)+
  geom_point(data.as[!c(data.as$test_category=="ecol_relev"),], mapping=aes(x=tempTestK1000, y=log(mass_specas), fill= tempTest, size = lnBWg), color = "black", alpha=1, pch=21, show.legend = FALSE)+
  scale_fill_gradient( low = cols.as[3], high = cols.as[1])+
  annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[AS] == .(round(AS_E_ER_eV,3))),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[AS] == .(round(AS_E_W_eV,3))),size=5, hjust=0, family="Arial", color = cols.as[1])+
  ylim(-5.7,4)+
  scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
  geom_line(data = data.plotAS_ER[which(round(data.plotAS_ER$lnBWg, 1) == round(log(1.35), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg ) ,color="black", linewidth=1, lty=1, show.legend=FALSE)+
  geom_line(data = data.plotAS_warm[which(round(data.plotAS_warm$lnBWg, 1) == round(log(1.4), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg), color = cols.as[2], linewidth=1, lty=1, show.legend=FALSE)
ggformat(ASmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/K)), y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1~g^-1)), print = T)

Arh.plot<-cowplot:::plot_grid(AMRmodel_plot3.E_ALL, RMRmodel_plot3.E_ALL, ASmodel_plot3.E_ALL, 
          align = "hv",
          axis = "l",
          nrow = 1,
          ncol = 3,
          label_size = 17)
ggsave(filename = paste("./Figures/Fig_ArrheniusFigMMR-RMR-AS_", Sys.Date(), ".png", sep=""),
       plot=Arh.plot, width = 13, height = 5, units = "in")


# Both MMR and RMR together - AS punchline plots ---------
# Overall all fish together:
MRmodel_plot1<-ggplot() +
  geom_ribbon(data=data.plotRMR_ER, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = "grey30")+
  geom_ribbon(data.plotAMRint_ER, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = "grey30")+
  geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTestK1000,2)==3.49,], mapping = aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="black", size=0.7, lty=1, show.legend=FALSE) +
  geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,2)==3.49,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=0.7, lty=1, show.legend=FALSE) +
  scale_y_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))+
  scale_x_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))
ggformat(MRmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), print = T)

MRmodel_plotW<-ggplot() +
  geom_ribbon(data.plotAMR_warm, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = cols.amr[2])+
  geom_ribbon(data=data.plotRMR_warm, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = cols.rmr[2])+
  geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTestK1000,2)==3.49,], mapping = aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color=cols.amr[1], size=0.7, lty=1, show.legend=FALSE) +
  geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTestK1000,2)==3.49,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color=cols.rmr[1], size=0.7, lty=1, show.legend=FALSE) +
  scale_y_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))+
  scale_x_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))
ggformat(MRmodel_plotW, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), pr = T)


FinalMR<-cowplot:::plot_grid(MRmodel_plot1, MRmodel_plotW,
          align = "hv",
          axis = "l",
          nrow = 1,
          ncol = 2,
          label_size = 17)
ggsave(filename = paste("./Figures/Fig_Final_", Sys.Date(), ".png", sep=""),
       plot=FinalMR, width = 7.5, height = 3.5, units = "in")


# Overall all fish together:
MRmodel_plot2<-ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
  geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=0.7, lty=1, show.legend=FALSE) +
  geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=0.7, lty=1, show.legend=FALSE) +
  geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color=cols.rmr[2], size=1, lty=1, show.legend=FALSE) +
  geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color=cols.amr[2], size=1, lty=1, show.legend=FALSE) +
   ylim(x = -5, 12)+
   xlim(x = -5, 12)
ggformat(MRmodel_plot2, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), print = F)

ggsave(filename = paste("./Figures/Fig5-final_MMR-RMR_", Sys.Date(), ".png", sep=""),
       plot=MRmodel_plot2, width = 4, height = 4, units = "in")

# Ecology specific groupings (sufficient for scaling ):  -------
ecol_rmr_sum1<-ggplot(ecology_data.RMRd.g, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)), show.legend = FALSE )+
  ylim(x = -7,12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.rmr[2]))
  # facet_grid(.~test_category3)
ggformat(ecol_rmr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = F)
ggsave(filename = paste("./Figures/FigVAR_ecolRMR_v2_", Sys.Date(), ".png", sep=""),
       plot=ecol_rmr_sum1, width = 4, height = 4, units = "in")

ecol_amr_sum1<-ggplot(ecology_data.AMRd.g, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)) , show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.amr[2]))
  # facet_grid(.~test_category3)
ggformat(ecol_amr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), print = F)
ggsave(filename = paste("./Figures/FigVAR_ecolAMR_v2_", Sys.Date(), ".png", sep=""),
       plot=ecol_amr_sum1, width = 4, height = 4, units = "in")

# Species specific groupings: 
species_amr_sum1<-ggplot(species.amr, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)) , show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.amr[2]))
  # facet_grid(.~test_category3)
ggformat(species_amr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), print = F)
ggsave(filename = paste("./Figures/FigVAR_speciesAMR_v2_", Sys.Date(), ".png", sep=""),
       plot=species_amr_sum1, width = 4, height = 4, units = "in")

species_rmr_sum1<-ggplot(species.rmr, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)) , show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.rmr[2]))
  # facet_grid(.~test_category3)
ggformat(species_rmr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = F )
ggsave(filename = paste("./Figures/FigVAR_speciesRMR_v2_", Sys.Date(), ".png", sep=""),
       plot=species_rmr_sum1, width = 4, height = 4, units = "in")





# SUPPLEMENTAL-----------
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

# Boxplots: FAS, AS, MR, mass-independent models ----------

fas_boxplot<-ggplot(data=data.fas, aes(y=FAS, x = test_category, label=species, fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.fas)+
  scale_x_discrete(labels=c("acclim" = "Acclimated \n warm", "acute" = "Acute \n warm", "ecol_relev" = "Optimal"))
ggformat(fas_boxplot, y_title = expression(FAS~(MMR/RMR)), x_title = "", print = F)
fas_boxplot<-fas_boxplot+theme(legend.position = "none")

rmr_boxplot<-ggplot(data=data.rmr, aes(y=mass_specrmr, x = test_category, label=species, fill=test_category))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.rmr)+
  scale_x_discrete(labels=c("acclim" = "Acclimated \n warm", "acute" = "Acute \n warm", "ecol_relev" = "Optimal"))
ggformat(rmr_boxplot, y_title = expression(RMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
rmr_boxplot<-rmr_boxplot+theme(legend.position = "none")

amr_boxplot<-ggplot(data=data.amr, aes(y=mass_specamr, x = test_category, label=species, fill=test_category))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.amr)+
  scale_x_discrete(labels=c("acclim" = "Acclimated \n warm", "acute" = "Acute \n warm", "ecol_relev" = "Optimal"))
ggformat(amr_boxplot, y_title = expression(MMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
amr_boxplot<-amr_boxplot+theme(legend.position = "none")

cowplot:::plot_grid(amr_boxplot,rmr_boxplot,fas_boxplot,
          align = "hv",
          axis = "l",
          nrow = 1,
          ncol = 3) %>%
ggsave(filename = "./Figures/FigSUP_boxplots_mar2022.png", width = 12.5, height =4)



# COMBO warm fish and ecol relev ********************
rmr_boxplot<-ggplot(data=data.rmr, aes(y=mass_specrmr, x = test_category3, label=species, fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=c("grey50", cols.rmr[3]))+
  scale_x_discrete(labels=c("warm" = "Warm", "ecol_relev" = "Optimal"))
ggformat(rmr_boxplot, y_title = expression(RMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
rmr_boxplot<-rmr_boxplot+theme(legend.position = "none")
ggsave(filename = paste("./Figures/Fig_box_rmr1_", Sys.Date(), ".png", sep=""),
       plot=rmr_boxplot, width = 4, height = 4, units = "in")

amr_boxplot<-ggplot(data=data.amr, aes(y=mass_specamr, x = test_category3, label=species, fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=c("grey50", cols.amr[3]))+
  scale_x_discrete(labels=c("warm" = "Warm", "ecol_relev" = "Optimal"))
ggformat(amr_boxplot, y_title = expression(MMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
amr_boxplot<-amr_boxplot+theme(legend.position = "none")
ggsave(filename = paste("./Figures/Fig_box_amr1_", Sys.Date(), ".png", sep=""),
       plot=amr_boxplot, width = 4, height = 4, units = "in")

fas_boxplot<-ggplot(data=data.fas, aes(y=FAS, x = test_category3, label=species, fill=test_category3))+
  geom_boxplot(show.legend = FALSE)+
  scale_fill_manual(values=c("grey50", cols.fas[3]))+
  scale_color_manual(values=cols.fas)+
  scale_x_discrete(labels=c("warm" = "Warm", "ecol_relev" = "Optimal"))
ggformat(fas_boxplot, y_title = "FAS (MMR / RMR)", x_title = "", print = F)
fas_boxplot<-fas_boxplot+theme(legend.position = "none")
ggsave(filename = paste("./Figures/Fig_box_fas1_", Sys.Date(), ".png", sep=""),
       plot=fas_boxplot, width = 4, height = 4, units = "in")


# Plots: FAS, AS with temp and size ---------------

# All scopes by temperature
AS_temp<-ggplot(data=data.as, aes(y=AS/BW_g^0.81, x=tempTest, size=BW_g, fill=tempTest,  color=tempTest))+
  geom_point(position=position_jitter(), pch=21, alpha=1, show.legend = FALSE)+
  # geom_boxplot(alpha=0, colour="black", show.legend = FALSE)+
  scale_fill_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_color_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  facet_grid(.~test_category)+
  ggtitle("AS scaled to b=0.82")
ggformat(AS_temp, y_title = expression(Aerobic~scope~(mgO[2]~g^-1~h^-1)), x_title = expression(Temperature~(degree*C)) , print=F)

# AS_temp <- AS_temp+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

FAS_temp<-ggplot(data=data.amr, aes(y=lnFAS, x=tempTest, size=BW_g, fill=tempTest, color=tempTest))+
  geom_point(position=position_jitter(), pch=21, alpha=0.8, show.legend = FALSE)+
  scale_fill_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_color_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  facet_wrap(DemersPelag~test_category3)+
  geom_hline(yintercept = 3, color="black", lty=2)+
  scale_size(breaks=c(10, 1000,75000), range=c(1,8))
  # geom_smooth(formula = y ~ poly(x, 2), color="#845BB7", show.legend = FALSE)
ggformat(FAS_temp, y_title = expression(Factorial~aerobic~scope), x_title = expression(Temperature~(degree*C)) , print=F)
# ecolFAS1 <- ecolFAS1+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

FAS_size<-ggplot(data=data.amr, aes(y=FAS, x=lnBWg, fill=tempTest, color = tempTest))+
  geom_point( pch=21, alpha=0.8, show.legend = FALSE)+
  # geom_boxplot(alpha=0, colour="black", show.legend = FALSE)+
  facet_wrap(DemersPelag~test_category3)+
  ylim(0,12)+
  scale_fill_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  scale_color_gradient2(low="blue", high="black", mid = "red", midpoint = 20)+
  geom_hline(yintercept = 3, color="black", lty=2)
ggformat(FAS_size, y_title = expression(Factorial~aerobic~scope), x_title = expression(italic(ln)~BW~(g)) , print=F)
# ecolFAS1 <- ecolFAS1+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

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


# Ecologies -------
plot_hist2_amrlm<-ggplot(ecology_data.AMRd, aes(x=slope)) +
  geom_vline(xintercept = mean(ecology_data.AMRd_er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = mean(ecology_data.AMRd_w$slope), color="red", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= 1.3, y = 4.6), label = paste(round(mean(ecology_data.AMRd_er$slope),2), " (", round(sd(ecology_data.AMRd_er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= 1.3, y = 4), label = paste(round(mean(ecology_data.AMRd_w$slope),2), " (", round(sd(ecology_data.AMRd_w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_histogram(show.legend = FALSE, bins=30, fill = cols.amr[1] , alpha = 1)+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,5)+
  scale_x_continuous(limits=c(0.5, 1.5))+
  facet_wrap(.~test_category3)
ggformat(plot_hist2_amrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist2_amrlm <- plot_hist2_amrlm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))

plot_hist2_rmrlm<-ggplot(ecology_data.RMRd, aes(x=slope)) +
  geom_vline(xintercept = mean(ecology_data.RMRd_er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = mean(ecology_data.RMRd_w$slope), color="red", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= 1.3, y = 4.6), label = paste(round(mean(ecology_data.RMRd_er$slope),2), " (", round(sd(ecology_data.RMRd_er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= 1.3, y = 4), label = paste(round(mean(ecology_data.RMRd_w$slope),2), " (", round(sd(ecology_data.RMRd_w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_histogram(show.legend = FALSE, bins=30, fill = cols.rmr[1] , alpha = 1)+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,5)+
  scale_x_continuous(limits=c(0.4, 1.5))+
  facet_wrap(.~test_category3)
ggformat(plot_hist2_rmrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist2_rmrlm <- plot_hist2_rmrlm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))


plot_hist2_faslm<-ggplot(ecology_data.FASd, aes(x=slope)) +
  geom_vline(xintercept = mean(ecology_data.FASd_er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = 0, color="grey30", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = mean(ecology_data.FASd_w$slope), color="red", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= -0.5, y = 5.5), label = paste(round(mean(ecology_data.FASd_er$slope),2), " (", round(sd(ecology_data.FASd_er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= -0.5, y = 4.5), label = paste(round(mean(ecology_data.FASd_w$slope),2), " (", round(sd(ecology_data.FASd_w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_histogram(show.legend = FALSE, bins=40, fill = cols.fas[1], alpha = 1)+
  scale_color_manual(values = c("black", "red"))+
  # ylim(0,12)+
  # scale_x_continuous(limits=c(-0.1, 1.7))+
  facet_wrap(.~test_category3)
ggformat(plot_hist2_faslm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist2_faslm <- plot_hist2_faslm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))


ecolhist.plot<-cowplot:::plot_grid(plot_hist2_amrlm, plot_hist2_rmrlm,plot_hist2_faslm,
                                   align = "hv",
                                   axis = "l",
                                   nrow = 3,
                                   ncol = 1) 
save_plot(filename = paste("./Figures/Suppl_ecol_hist", Sys.Date(), ".png", sep=""),
          ecolhist.plot, base_width = 8, base_height = 12, units = "in")
# 
# 
# # good for scaling only & pairs:
# scopeWARM.ER.good<-ggplot(ecology_data.AMRd.g,
#                           aes(x=slope, y=ecol_temp_cat),
#                           color=test_category3, fill=MR,
#                           group = ecology_subgroup)+
#   # geom_text(aes(label = n_data_n_species,  x = 0.61),   family = "Arial", size=4, fontface="bold")+
#   # geom_text(data = ecology_data.AMRd.g, aes(label = n_data_n_species,  x = 0.5),  family = "Arial", size=4, fontface="bold")+
#   geom_linerange(data = ecology_data.RMRd.g,
#                  aes(xmin = as.numeric(slope.ciL) , xmax = as.numeric(slope.ciH)))+
#   geom_point(pch=21, size=1, stroke=1, fill="grey", alpha=1)+
#   geom_linerange(data = ecology_data.AMRd.g, mapping=aes(xmin = slope.ciL, xmax = slope.ciH), alpha=1)+
#   geom_point(data=ecology_data.AMRd.g, mapping = aes(x=slope, y=ecol_temp_cat, color=test_category3),
#              pch=23, size=1, stroke=1, fill="black", alpha=1)+
#   geom_line(arrow = arrow(length=unit(0.30,"cm"), ends="first", type = "closed"), size = 1)+
#   scale_color_manual(values=c("black", "red"))+
#   xlab(expression(Slope~value~(italic(b))))+
#   scale_x_continuous(limits = c(0.46,1.1), breaks = c(0.66, 0.75, 0.82, 0.9, 1))+
#   theme_classic()+
#   geom_vline(xintercept = c(0.66, 0.75,1), lty = "dashed", color = "grey", size = 0.6)+
#   theme(axis.text.y = element_text(face = "italic", color = "black", size = 15),
#         axis.text.x = element_text( color = "black", size = 15),
#         axis.title.y = element_blank(),
#         legend.position = "none",
#         axis.line.y=element_line(colour = 'black',size=0.5),
#         axis.line.x=element_line(colour = 'black',size=0.5),
#         axis.ticks.y=element_line(size=0.5),
#         axis.ticks.x=element_line(size=0),
#         text=element_text(size=20,  family="Arial"))
# scopeWARM.ER.good
# ggsave(filename = paste("./Figures/Figxx_ecology1_ScalingSuited", Sys.Date(), ".png", sep=""),
       # plot=scopeWARM.ER.good, width = 8, height = 8, units = "in")


rmr.WARM.ER.good<-ggplot(data=ecology_data.RMRd.g,
                         aes(x=slope,
                             y=factor(ecol_temp_cat, level = c(as.character(ecol_temp_cat))),
                             color=interaction(test_category3, MR),
                             fill=interaction(test_category3, MR), 
                             group= ecology_subgroup))+
  scale_y_discrete(limits=rev,
    labels = c("",
               "Any salinity", 
               "",
               "Marine",
               "",
               "Reef-associated",
              "",
               "Demersal",
              "",
               "Benthopelagic", 
               "",
               "Tropical",
               "",
               "Temperate",
               "",
               "Suptropical",
               "",
               "Fusiform"))+
  geom_text(aes(label = slope,  x = 0.25),
            family = "Arial", fontface ="bold", size=3.5, color = cols.rmr[1], hjust = 1)+
  geom_text(data = ecology_data.AMRd.g, aes(label = slope,  x = 0.3),
            family = "Arial", fontface ="bold", size=3.5, color = cols.amr[1], hjust = 0)+
  
  geom_text(data = ecology_data.RMRd.g,
            mapping = aes(label =  n_data_n_species, x = 1.5),
            family = "Arial", size=3.5, color = cols.rmr[1], hjust = 1)+
  geom_text(data = ecology_data.AMRd.g,
            mapping = aes(label =  n_data_n_species, x = 1.55),
            family = "Arial", size=3.5, color = cols.amr[1], hjust = 0)+
  
  geom_linerange(data = ecology_data.RMRd.g,
                 aes(xmin = as.numeric(slope.ciL),
                     xmax = as.numeric(slope.ciH)), alpha = 0.4)+
  geom_vline(xintercept = 1, lty = "dotted", color="grey")+
  geom_linerange(aes(xmin = slope.ciL, xmax = slope.ciH),
                 alpha = 0.4)+
  geom_linerange(data = ecology_data.AMRd.g,
                 mapping = aes(xmin = as.numeric(slope.ciL), xmax = as.numeric(slope.ciH)),
                 alpha = 0.4)+
  geom_point(pch=21, size=4, stroke=0.6, alpha=1)+
  geom_point(data = ecology_data.AMRd.g,
             pch=21, size=4, stroke=0.6, alpha=1)+
  geom_line(arrow = arrow(length=unit(0.20,"cm"), ends="last", type = "closed"),
            size = 0.7, color = cols.rmr[1])+
  geom_line(data = ecology_data.AMRd.g,
            arrow = arrow(length=unit(0.20,"cm"), ends="last", type = "closed"),
            size = 0.7, color = cols.amr[1])+
  scale_color_manual(values=c("black", "black",cols.amr[1], cols.rmr[1]))+
  scale_fill_manual(values=c(cols.amr[1],cols.rmr[1],cols.amr[2], cols.rmr[2]))+
  xlab(expression(Slope~value~(italic(b))))+
  scale_x_continuous(limits = c(0.1,1.8))+
  # scale_y_discrete(limits=rev)+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 15),
        axis.text.x = element_text( color = "black", size = 15),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.line.y=element_line(colour = 'black',size=0.5),
        axis.line.x=element_line(colour = 'black',size=0.5),
        axis.ticks.y=element_line(size=0.5),
        axis.ticks.x=element_line(size=0),
        text=element_text(size=20,  family="Arial"))
rmr.WARM.ER.good


fas.WARM.ER.good<-ggplot(data=ecology_data.FASd.g,
                         aes(x=slope,
                             y=factor(ecol_temp_cat, level = c(as.character(ecol_temp_cat))),
                             color=test_category3,
                             fill=test_category3, 
                             group= ecology_subgroup))+
  scale_y_discrete(limits=rev,
    labels = c("",
               "Any salinity", 
               "",
               "Marine",
               "",
               "Reef-associated",
              "",
               "Demersal",
              "",
               "Benthopelagic", 
               "",
               "Tropical",
               "",
               "Temperate",
               "",
               "Suptropical",
               "",
               "Fusiform"))+
  geom_text(aes(label = n_data_n_species,  x = 0.2), family = "Arial", size=3.5, hjust=0)+
  geom_text(aes(label = slope, x = -0.2, family = "Arial"), size=3.5, hjust=1, fontface="bold")+
  geom_linerange(data = ecology_data.FASd.g,
                 aes(xmin = as.numeric(slope.ciL) , xmax = as.numeric(slope.ciH)), alpha = 0.4)+
  geom_vline(xintercept = 0, lty = "dotted", color="grey")+
  geom_linerange(aes(xmin = slope.ciL, xmax = slope.ciH), alpha = 0.4)+
  geom_point(pch=21, size=4, stroke=0.6, alpha=1)+
  geom_line(arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "closed"), size = 0.7)+
  scale_color_manual(values=c("black", cols.fas[1]))+
  scale_fill_manual(values=c("black", cols.fas[2]))+
  xlab(expression(Slope~value~(italic(b))))+
  scale_x_continuous(limits = c(-0.3,0.4))+
  # scale_y_discrete(limits=rev)+
  theme_classic()+
  theme(
    # axis.text.y = element_text(face = "italic", color = "black", size = 15),
        axis.text.x = element_text( color = "black", size = 15),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.line.y=element_line(colour = 'black',size=0.5),
        axis.line.x=element_line(colour = 'black',size=0.5),
        axis.ticks.y=element_line(size=0.5),
        axis.text.y = element_blank(),
        axis.ticks.x=element_line(size=0),
        text=element_text(size=20,  family="Arial"))
fas.WARM.ER.good

cowplot::plot_grid(rmr.WARM.ER.good,
                   fas.WARM.ER.good,
                  nrow = 1, 
                  labels = "AUTO", 
                  rel_widths = c(1, 0.5),
                  label_x = c(0.435, 0.25),
                  label_y = c(0.98, 0.98)) %>% 
  ggsave(filename = paste("./Figures/Fig_ecology1_ScalingSuited_FAS", Sys.Date(), ".png", sep=""),
         width = 9.5, height = 6)


# Violin plots: AS and FAS ecology size independent -------------

# Demersal Pelagic
ecolAS1<-ggplot(data=data.as, aes(y=mass_specas, x=DemersPelag, size=BW_g, fill=test_category3, group = interaction(test_category3, DemersPelag), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  scale_size(breaks=c(10, 1000,75000), range=c(1,8))+
  scale_x_discrete(labels=c("Bentho- \n pelagic", "Demersal", "Pelagic", "Reef- \n associated"))+
  ylim(0,5)+
  labs(fill = expression(degree*C), color=  expression(degree*C), size = "g")+
  annotate("text", label = "Lifestyle", x= 3.25, y=4.8, hjust=-0.5, size=5)+
  annotate("text", label = "*** (ecol subgr)",  x= 0.6, y=4.8, hjust=-0.5, size=3.5)
ggformat(ecolAS1, y_title = bquote("AS" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=FALSE)
ecolAS1 <- ecolAS1 + theme(
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"), 
  axis.text.y = element_text( color = "black", size = 15), 
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "right",
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5))
# ecolAS1

ecolFAS1<-ggplot(data=data.fas, aes(y=FAS, x=DemersPelag, size=BW_g, fill=test_category3, group = interaction(test_category3, DemersPelag), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  scale_size(breaks=c(10, 1000,75000), range=c(1,8))+
  scale_x_discrete(labels=c("Bentho- \n pelagic", "Demersal", "Pelagic", "Reef- \n associated"))+ 
  annotate("text", x = 3.5, y = 16, label = "Lifestyle", hjust=0, family="Arial", size=5)+
  ylim(0,16.5)+
  labs(fill = expression(degree*C), color=  expression(degree*C), size = "g")
ggformat(ecolFAS1, y_title = expression(FAS~(MMR/RMR)), x_title = element_blank() , print=FALSE)
ecolFAS1 <- ecolFAS1 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"), 
                             axis.text.y = element_text( color = "black", size = 15), 
                             text = element_text( color = "black", size = 15),
                             legend.position = "right",
                             legend.direction = "vertical",
                             legend.justification='center',
                             legend.margin=margin(0,0,0,0),
                             legend.box.margin=margin(0,0,6,0),
                             plot.margin = margin(5.5, 5.5, 5.35, 5.5))
# ecolFAS1

# Body shape
ecolAS2<-ggplot(data=data.as, aes(y=mass_specas, x=BodyShapeI, size=BW_g, fill=test_category3, group = interaction(test_category3, BodyShapeI), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = TRUE)+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  scale_fill_manual(values = c("black", "red"), guide = "none" )+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+# , guide = "none"
  scale_size(breaks=c(10, 1000,75000), range=c(1,8))+
  scale_x_discrete(labels=c("Elongated", "Fusiform", "Short- \n deep"))+
  labs(fill = expression(degree*C), color=  expression(degree*C), size = "g")+
  ylim(0,5)+
  annotate("text", label = "Morphology", x= 2.8, y=4.8, hjust=0.5, size=5)
ggformat(ecolAS2, y_title =  bquote("Aerobic scope" ~ (mgO[2] ~ g^- 1 ~ h^-1)), x_title = element_blank() , print=FALSE)
ecolAS2 <- ecolAS2 + theme(
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(), 
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"),
  axis.text.y = element_text( color = "black", size = 15, family="Arial"),
  legend.position = "none",
  plot.margin = margin(unit(c(5.5, 100,5.5,5.5), "points")))
# ecolAS2


# ecolAS2<- ecolAS2<-theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, family="Arial"), 
#                          axis.text.y = element_text(face = "italic", color = "black", size = 10, family="Arial"),
#                           text=element_text(size=15,  family="Arial"), 
#                          legend.position = "none")

ecolFAS2 <- ggplot(data=data.fas, aes(y=FAS, x=BodyShapeI, size=BW_g, fill=test_category3, group = interaction(test_category3, BodyShapeI), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  scale_x_discrete(labels=c("Elongated", "Fusiform", "Short-\n deep"))+
  ylim(0,16.5)+
  annotate("text", label = "Morphology", x= 2.8, y=16, hjust=0.5, size=5)
ggformat(ecolFAS2, y_title = expression(FAS~(MMR/RMR)), x_title = element_blank() , print=FALSE)
ecolFAS2 <- ecolFAS2 + theme( axis.title.x = element_blank(), 
                              axis.title.y = element_blank(), 
                              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"),
                              axis.text.y = element_text( color = "black", size = 15, family="Arial"),
                              text = element_text( color = "black", size = 15),
                              legend.position = "top",
                              legend.direction = "horizontal",
                              legend.justification='left',
                              plot.margin = margin(unit(c(5.5, 100,5.5, 5.5), "points")))
# ecolFAS2

# Climate
ecolAS3 <-ggplot(data=data.as, aes(y=mass_specas, x=Climate, size=BW_g, fill=test_category3, group = interaction(test_category3, Climate), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  geom_violin(show.legend = FALSE, alpha=0.4, color = "black")+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15, guide="none")+
  labs(fill = expression(degree*C), color=  expression(degree*C), size = "g")+
  ylim(0,5)+
  annotate("text", label = "Climate", x= 3.8, y=4.8, hjust=0.3, size=5)
  # annotate("text", label = "***",  x= 2.5,y=4.8, hjust=-0.5, size=5)
ggformat(ecolAS3, y_title =  bquote("AS" ~ (mgO[2] ~ g^- 1 ~ h^-1)), x_title = element_blank() , print=FALSE)
ecolAS3 <- ecolAS3 + theme(
  axis.title.x = element_blank(), 
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"),
  axis.text.y = element_text( color = "black", size = 15, family="Arial"),
  text = element_text( color = "black", size = 15),
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.justification='left')
# ecolAS3


ecolFAS3<- ggplot(data=data.fas, aes(y=FAS, x=Climate, size=BW_g, fill=test_category3, group = interaction(test_category3, Climate), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  ylim(0,16.5)+
  geom_violin(show.legend = FALSE, alpha=0.4, color = "black")+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  annotate("text", label = "Climate", x= 3.8, y=16, hjust=0.3, size=5)+
  annotate("text", label = "*** (ecol subgr)",  x= 0.6, y=16.1, hjust=-0.5, size=3.5)+
  annotate("text", label = "*** (temp categ)",  x= 0.57, y=14.9, hjust=-0.5, size=3.5)
ggformat(ecolFAS3, y_title = expression(FAS~(MMR/RMR)), x_title = element_blank() , print=FALSE)
ecolFAS3 <- ecolFAS3 + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"), 
  axis.text.y = element_text( color = "black", size = 15, family="Arial"),
  text = element_text( color = "black", size = 15),
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.justification='left')
# ecolFAS3

# Salinity
ecolAS4 <- ggplot(data=data.as, aes(y=mass_specas, x=salintyComb, size=BW_g, fill=test_category3, group = interaction(test_category3, salintyComb), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  scale_size(breaks=c(10, 1000,75000), range=c(1,8), guide="none")+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  labs(fill = expression(degree*C), color=  expression(degree*C), size = "g")+
  scale_x_discrete(labels=c("Freshwater" , "Freshwater-\n brackish", "Marine", "Marine-\n brackish" ,"All salinities"))+
  ylim(0,5)+
  annotate("text", label = "Salinity",  x= 4.2 ,y=4.8, hjust=-0.5, size=5)
ggformat(ecolAS4, y_title =  bquote("Aerobic scope" ~ (mgO[2] ~ g^- 1 ~ h^-1)), x_title = element_blank() , print=FALSE)
ecolAS4 <- ecolAS4 + theme(
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(), 
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"),
  axis.text.y = element_text( color = "black", size = 15, family="Arial"),
  text = element_text( color = "black", size = 15), 
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.justification='left')
# ecolAS4


ecolFAS4 <- ggplot(data=data.fas, aes(y=FAS, x=salintyComb, size=BW_g, color=tempTest, fill=test_category3, group = interaction(test_category3, salintyComb)))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  scale_x_discrete(labels=c("Freshwater" , "Freshwater-\n brackish", "Marine", "Marine-\n brackish" ,"All salinities"))+
  # annotate("text", label = "***",  x= 2.5,y=16, hjust=-0.5, size=5)+
  annotate("text", label = "Salinity",  x= 4.3 ,y=16, hjust=-0.5, size=5)+
  ylim(0,16.5)
ggformat(ecolFAS4, y_title = expression(FAS~(MMR/RMR)), x_title = element_blank() , print=FALSE)
ecolFAS4 <- ecolFAS4 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"), 
                             axis.title.x = element_blank(), 
                             axis.title.y = element_blank(),
                             text = element_text( color = "black", size = 15),
                             plot.title = element_text(face = "bold", size=15, hjust = 0.5),
                             legend.position = "none",
                             legend.direction = "horizontal",
                             legend.justification='left')
# ecolFAS4

plot_grid(ecolAS1, ecolAS2,
          ecolAS3, ecolAS4, 
          align = "h",
          labels = c('A', 'B', 'C', "D"),
          label_size = 15,
          nrow = 2,
          ncol=2,
          label_x = c(0.145, 0.08),
          label_y = c(0.89, 0.89))%>% 
  ggsave(filename = "./Figures/Figure_AS_mar2022-NOlegend.png", width = 10, height = 8)


plot_grid(ecolFAS1, ecolFAS2,
          ecolFAS3, ecolFAS4, 
          align = "h",
          labels = c('A', "B", 'C', "D"),
          label_size = 15,
          nrow = 2,
          ncol=2,
          label_x = c(0.16, 0.08),
          label_y = c(0.89, 0.89)) %>% 
  ggsave(filename = "./Figures/Figure_FAS_mar2022-NOlegend.png", width = 10, height = 8)


# SPECIES ------

# PLOTS:
plot_hist3_amrlm<-ggplot(species.amr, aes(x=slope)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  geom_vline(xintercept = mean(species.amr.er$slope), color="black", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= 2, y = 11.5), label = paste(round(mean(species.amr.er$slope),2), " (", round(sd(species.amr.er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= 2, y = 10.1), label = paste(round(mean(species.amr.w$slope),2), " (", round(sd(species.amr.w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_vline(xintercept = mean(species.amr.w$slope), color="red", lty="solid", lwd=1)+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.amr[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,12)+
  scale_x_continuous(limits=c(-1, 3))+
  facet_wrap(.~test_category3)
ggformat(plot_hist3_amrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist3_amrlm <- plot_hist3_amrlm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))


plot_hist3_rmrlm<-ggplot(species.rmr, aes(x=slope)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  geom_vline(xintercept = mean(species.rmr.er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = mean(species.rmr.w$slope), color="red", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= 1.7, y = 11.5), label = paste(round(mean(species.rmr.er$slope),2), " (", round(sd(species.rmr.er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= 1.7, y = 10.1), label = paste(round(mean(species.rmr.w$slope),2), " (", round(sd(species.rmr.w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.rmr[1] , alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,12)+
  scale_x_continuous(limits=c(-0.5, 2.5))+
  facet_wrap(.~test_category3)
ggformat(plot_hist3_rmrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist3_rmrlm <- plot_hist3_rmrlm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))


plot_hist3_faslm<-ggplot(species.fas, aes(x=slope)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  geom_vline(xintercept = mean(species.fas.er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = 0, color="grey30", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = mean(species.fas.w$slope), color="red", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= 2, y = 11.5), label = paste(round(mean(species.fas.er$slope),2), " (", round(sd(species.fas.er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= 2, y = 10.1), label = paste(round(mean(species.fas.w$slope),2), " (", round(sd(species.fas.w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_histogram(show.legend = FALSE, bins=40, fill = cols.fas[1], alpha = 1)+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,12)+
  scale_x_continuous(limits=c(-2, 3.5))+
  facet_wrap(.~test_category3)
ggformat(plot_hist3_faslm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist3_faslm <- plot_hist3_faslm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))

spechist.plot<-cowplot:::plot_grid(plot_hist3_amrlm, plot_hist3_rmrlm,plot_hist3_faslm,
                              align = "hv",
                              axis = "l",
                              nrow = 3,
                              ncol = 1) 
save_plot(filename = paste("./Figures/Suppl_species_hist", Sys.Date(), ".png", sep=""),
      spechist.plot, base_width = 8, base_height = 12, units = "in")


sp1_amr<-ggplot(data=species.amr, aes(y=slope, x=fct_reorder(species, (slope.ciL-slope.ciH)), color = test_category3))+
  geom_hline(yintercept = 0.81, colour="black", lty="solid")+
  ylim(min(species.amr$slope.ciL)-0.1, max(species.amr$slope.ciH)+0.1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciH, color = test_category3), size=1.2, alpha=1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciL, color = test_category3), size=1.2, alpha=1)+
  geom_point(pch=1, color="black", size=1)+
  geom_text(aes(label = n, y = -3.4, color = test_category3), hjust = "inward",  family = "Arial", size=3)+
  ylab(expression(Slope~value~(italic(b))))+
  scale_color_manual(values = c("black", "red"))+
  facet_grid(.~test_category3)+
  coord_flip()+
  theme_pubr()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 10),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.8, 5.5, 5.5), "pt"),
        text=element_text(size=20,  family="Arial"))
ggsave(filename = paste("./Figures/FigSUP_sp_slopes_amr_", Sys.Date(), ".png", sep=""),
       plot=sp1_amr, width = 6, height = 8.8, units = "in")


sp1_rmr<-ggplot(data=species.rmr, aes(y=slope, x=fct_reorder(species, (slope.ciL-slope.ciH)), color = test_category3))+
  geom_hline(yintercept = 0.81, colour="black", lty="solid")+
  ylim(min(species.rmr$slope.ciL)-0.1, max(species.rmr$slope.ciH)+0.1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciH, color = test_category3), size=1.2, alpha=1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciL, color = test_category3), size=1.2, alpha=1)+
  geom_point(pch=1, color="black", size=1)+
  geom_text(aes(label = n, y = -4, color = test_category3), hjust = "inward",  family = "Arial", size=3)+
  ylab(expression(Slope~value~(italic(b))))+
  scale_color_manual(values = c("black", "red"))+
  facet_grid(.~test_category3)+
  coord_flip()+
  theme_pubr()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 10),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.8, 5.5, 5.5), "pt"),
        text=element_text(size=20,  family="Arial"))
ggsave(filename = paste("./Figures/FigSUP_sp_slopes_rmr_", Sys.Date(), ".png", sep=""),
       plot=sp1_rmr, width = 8, height = 10, units = "in")


sp1_fas<-ggplot(data=species.fas, aes(y=slope, x=fct_reorder(species, (slope.ciL-slope.ciH)), color = test_category3))+
  geom_hline(yintercept = 0, colour="black", lty="solid")+
  ylim(min(species.fas$slope.ciL)-0.1, max(species.fas$slope.ciH)+0.1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciH, color = test_category3), size=1.2, alpha=1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciL, color = test_category3), size=1.2, alpha=1)+
  geom_point(pch=1, color="black", size=1)+
  geom_text(aes(label = n, y = -5, color = test_category3), hjust = "inward",  family = "Arial", size=3)+
  ylab(expression(Slope~value~(italic(b))))+
  scale_color_manual(values = c("black", "red"))+
  facet_grid(.~test_category3)+
  coord_flip()+
  theme_pubr()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 10),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.8, 5.5, 5.5), "pt"),
        text=element_text(size=20,  family="Arial"))
ggsave(filename = paste("./Figures/FigSUP_sp_slopes_fas_", Sys.Date(), ".png", sep=""),
       plot=sp1_fas, width = 6, height = 8, units = "in")

