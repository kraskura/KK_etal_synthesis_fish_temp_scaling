
library(here)
source(here("R", "colors_themes.R"))

data<-read.csv("./Data/Summary_scaling_b_mmr_rmr_jul192022.csv")
data$Running_ID<-factor(data$Running_ID)
data$star<-NA
for(i in 1:nrow(data)){
  if(is.na(data$Order_MMR[i])){
    data$Order_MMR[i]<- data[which(data$Running_ID[i] == data$Running_ID)[1],"Order_MMR"]
  }
  if(!data$notes.mechan[i] == "REPEAT" & !data$notes.mechan[i] == ""){
    data$star[i]<-"•"
  }
}

# data

plot <- ggplot(data = data, mapping = aes(x = b, y=factor(Order_MMR), size = MR,
                                    Color = Test_C, fill = Test_C,
                                    shape = MR, group = Running_ID))+
  geom_vline(xintercept = c(1, 0.75), lty = "dashed", color="black", linewidth= 0.3)+
  geom_path(color = "grey", arrow = arrow(length = unit(0.1, "cm"), angle = 1), 
            show.legend = F)+
  geom_point(alpha=1)+
  xlim(0, 1.55)+
  geom_text(aes(label = Study, x = 0.05), family = "Helvetica", size=3, hjust = 0)+
  geom_text(aes(label = star, x = 0.0), family = "Helvetica", size=7, color = "red",
            hjust = 0)+
  scale_color_viridis_c(option = "A")+
  scale_fill_viridis_c(option = "A", name = "T, ºC")+
  xlab(expression(Slope~value~(italic(b))))+
  scale_shape_manual(values = c(21, 21), name = "")+
  scale_size_manual(values = c(3, 1), name = "")+
    scale_y_discrete(
      labels = c("Brown trout"  , "Leopard coral grouper", "Mudfish "   ,          
                 "Black carp"   ,         "Crucian carp"    ,      "European perch",       
                 "Leopard coral grouper" ,"Atlantic cod"     ,     "Barramundi"    ,       
                 "European perch"   ,     "European perch"    ,    "Atlantic salmon" ,     
                  "Grass carp"     ,       "Atlantic cod"     ,     "Brown trout"    ,      
                  "Brown trout"   ,        "Brown trout"     ,      "European perch" ,      
                  "Atlantic cod"  ,        "Brown trout"      ,     "Atlantic cod"     ,    
                  "Cunner"       ,         "Atlantic cod"      ,    "Cunner"            ,   
                  "Zebrafish"    ,         "Barramundi"        ,    "Cunner"             ,  
                  "Cunner"      ,          "Sockeye salmon"    ,    "Walleye"             , 
                  "Cunner"       ,        "Northern pike"      ,   "Common killifish "    ,
                  "Rainbow trout"   ,      "Sockeye salmon"    ,    "Sockeye salmon"       ,
                  "Banded kōkopu"    ,     "Barramundi"        ,    "Rainbow trout"        ,
                  "Common killifish "   ,  "Common killifish " ,    "Common killifish "    ,
                  "Common killifish "))+  
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 9.5),
        axis.text.x = element_text( color = "black", size =10),
        axis.title.y = element_blank(),
        legend.position = c(0.78, 0.1),
        axis.line.y=element_line(colour = 'black',size=0.5),
        axis.line.x=element_line(colour = 'black',size=0.5),
        axis.ticks.y=element_line(size=0.5),
        axis.ticks.x=element_line(size=0),
        text=element_text(size=12,  family="Helvetica"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.margin = margin(-0.9,0,0,0, unit="cm"), legend.box = "horizontal")+
  guides(size = "legend", colour = "legend") 

  # plot 
  
ggsave(filename = paste("./Figures/Figure1_main.png", sep=""),
       width = 6, height = 7.5, units = "in")
 


