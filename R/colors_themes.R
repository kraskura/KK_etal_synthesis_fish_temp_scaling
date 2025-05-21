# colors 
# 
library(colorBlindness)
library(patchwork)
# displayAvailablePalette(color="white")
# colorBlindness::SteppedSequential5Steps

cols.amr<-c(colorBlindness::SteppedSequential5Steps[c(21, 22, 23, 24,  25)]) 
cols.rmr<-c(colorBlindness::SteppedSequential5Steps[c(11, 12, 13,14,  15)]) # "#6B990F" "#A3CC51" "#E5FFB2"
cols.fas<-c(colorBlindness::SteppedSequential5Steps[c(1,2, 3,4, 5)])  # "#990F0F" "#CC5151" "#FFB2B2"
cols.as<-c(colorBlindness::SteppedSequential5Steps[c(16,17,  18, 19, 20)]) 


cols<-c(cols.amr[1],cols.rmr[1],cols.amr[3],cols.rmr[3], cols.amr[5],cols.rmr[5], cols.as[2], cols.fas[2])
        # AMR -rmr- AMR dark - rmr dark - light - as - fas