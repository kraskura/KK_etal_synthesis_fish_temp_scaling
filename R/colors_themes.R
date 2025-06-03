# colors 
# if neccessary : 
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("ggtree")

# install.packages("lme4", type = "source")
# install.packages("evolvability", type = "source")
library(ggh4x)
library(tidyverse)
library(dplyr)

library(colorBlindness)
library(patchwork)
library(weathermetrics)

library(lme4)
library(emmeans)
library(car)

library(evolvability) # Almer function
library(ape)
library(rotl)
library(ggtree)
library(Matrix)

# library(kableExtra)

library(pryr)
library(viridisLite)
library(TDbook)
library(ggimage)

library(ggformat2) # from github
library(weathermetrics)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(forcats)
library(patchwork)
library(reshape2)
library(emmeans)
library(systemfonts)
library(flextable)


library(here)

# displayAvailablePalette(color="white")
# colorBlindness::SteppedSequential5Steps

cols.amr<-c(colorBlindness::SteppedSequential5Steps[c(21, 22, 23, 24,  25)]) 
cols.rmr<-c(colorBlindness::SteppedSequential5Steps[c(11, 12, 13,14,  15)]) # "#6B990F" "#A3CC51" "#E5FFB2"
cols.fas<-c(colorBlindness::SteppedSequential5Steps[c(1,2, 3,4, 5)])  # "#990F0F" "#CC5151" "#FFB2B2"
cols.as<-c(colorBlindness::SteppedSequential5Steps[c(16,17,  18, 19, 20)]) 


cols<-c(cols.amr[1],cols.rmr[1],cols.amr[3],cols.rmr[3], cols.amr[5],cols.rmr[5], cols.as[2], cols.fas[2])
        # AMR -rmr- AMR dark - rmr dark - light - as - fas