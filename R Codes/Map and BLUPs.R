setwd("E:\\Fatma\\Stability Yield")
library(ggplot2)
library(maps)
library(ggrepel)

locations <- read.table("E:\\Fatma\\Stability Yield\\envRtype.txt", header = T)
head(locations)

us_map <- map_data("state")

p <- ggplot() +
  geom_polygon(data = us_map, aes(x=long, y=lat, group=group), fill="white", color="gray90") +
  geom_point(data=locations, aes(x=lon, y=lat, color= Tester),size=2) +
  geom_text_repel(data=locations, aes(x=lon, y=lat, label=location),size=3,
                  position=position_jitter(width=0.5,height=0.5))+
  theme_bw() +
  labs(title="") +
  facet_grid(~year)+
  scale_color_brewer(palette = "Set1")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(#legend.background = element_blank(),
    legend.background =  element_rect(fill = alpha("gray80", 0.3) , size=0.1, linetype="solid"),
    legend.key.height = unit(0.1, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.59, 0.15))+
  scale_x_continuous("Longitude")+
  scale_y_continuous("Latitude")
p

jpeg("US.map.jpeg", width = 11,height =4.5,units = "in", res=500)
p
dev.off()







###########################
######## Yield ############
###########################
rm(list=ls())
setwd("E:/Fatma")
library(nlme)
library(sjstats)
library(lme4)
library("RColorBrewer")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(dplyr)

df20 <- read.csv("g2f_2020_phenotypic_clean_data.csv")
df20 <- df20[df20$Field.Location %in% c("DEH1","GAH1","INH1","MNH1","NEH1","NEH2","NEH3","OHH1","TXH1","TXH2","WIH1","WIH2"),]

df21 <- read.csv("g2f_2021_phenotypic_clean_data.csv")
df21 <- df21[df21$Field.Location %in% c( "DEH1", "GAH1", "ILH1", "MNH1", "NCH1", "NEH1", "SCH1", "TXH2", "TXH3", "WIH1", "WIH2", "WIH3",
                                         "IAH1", "IAH2", "IAH3", "IAH4", "MIH1"),]

df <- rbind(df20,df21)
df$Env <- paste(df$Field.Location,df$Year, sep = ".")
df <- df[!df$Replicate %in% c(3,4),]


df$Range <- as.factor(df$Range)
df$Pass <- as.factor(df$Pass)
df$Replicate  <- as.factor(df$Replicate )
df$Pedigree <- as.factor(df$Pedigree)
df$Env <- as.factor(df$Env)
df$Grain.Yield..bu.A. <- as.numeric(df$Grain.Yield..bu.A.)


anova<- lmer(Grain.Yield..bu.A.*.0673 ~
               (1|Pedigree)+
               (1|Env)+
               (1|Pedigree:Env)+
               (1|Range)+
               (1|Pass)+
               (1|Replicate:Env), df )  
cv(anova)
rmse<- rmse(anova)
R<- MuMIn::r.squaredGLMM(anova)[,2]
VC<-as.data.frame(print(VarCorr(anova ), comp=c("Variance")))
VC$Percent<-round(VC$vcov/sum(VC$vcov)*100,2)
VC
heritability <- VC[2,6] / (VC[2,6] +  VC[1,6]/length(unique(df$Env)) +  (VC[7,6]/length(unique(df$Env))*2))
round(heritability,3)  
VC[,7] <-rmse
VC[,8] <-heritability
VC[,9] <-R
names(VC)[7:9] <- c("Rmse","Heritability","R")
VC
write.csv(VC, "VC_yield.csv")


E <- ranef(anova)$'Env'
G <- ranef(anova)$'Pedigree'
GE <- coef(anova)$'Pedigree:Env'
GE[,2] <- rownames(GE)
head(GE)

GE$Pedigree <- lapply(strsplit(as.character(GE$V2), "\\:"), "[", 1)
GE$Env <- lapply(strsplit(as.character(GE$V2), "\\:"), "[", 2)
head(GE)
GE <- as.data.frame(lapply(GE, unlist))
head(GE)
names(GE)[1:2] <- c("Yield","interaction" )
head(GE)

for (k in 1:nrow(GE)){
  GE[,1][k]<-GE[,1][k]+E$`(Intercept)`[row.names(E)==GE$Env[k]]
}

for (k in 1:nrow(GE)){
  GE[,1][k]<-GE[,1][k]+G$`(Intercept)`[row.names(G)==GE$Pedigree[k]]
}

head(GE)

GE <- GE[grepl( 'W10004', GE$Pedigree   ), ]
GE$Tester <- lapply(strsplit(as.character(GE$Pedigree), "\\/"), "[", 2)
GE <- as.data.frame(lapply(GE, unlist))

GE$Location <- lapply(strsplit(as.character(GE$Env), "\\."), "[", 1)
GE$Year <- lapply(strsplit(as.character(GE$Env), "\\."), "[", 2)
GE <- as.data.frame(lapply(GE, unlist))
Yield <- GE
mean <- Yield %>% group_by(Env) %>% summarise(Mean = mean(Yield)) %>% as.data.frame()


library(tidytext)
head(Yield)

p<-ggplot(Yield, aes(x=reorder(Env, -Yield), y=Yield, fill=Tester)) +
  geom_boxplot(position = position_dodge2(preserve = "single"), size=0.2, outlier.size = 0.3)+
  theme_bw()+
  #  facet_grid(Year~Tester)+
  scale_x_discrete("Environment")+
  scale_y_continuous("Grain yield (t/ha)")+
  theme(legend.position = c(0.93,0.83),
        axis.text.x = element_text(size=8, angle = 45, hjust = 1),
        legend.background = element_rect(fill="lightblue",
                                         size=0.5,
                                         linetype="solid"))
#  ggtitle("B) Genotypic effect of hybrid for grain yield (t/ha)")
p

VC$Trait <- "Grain yield (t/ha)"
names(VC)[8] <- "Repeatability"

VC$grp <- gsub("Pedigree", "g", VC$grp)
VC$grp <- gsub("Env", "E", VC$grp)
VC$grp <- gsub("Replicate", "Rep", VC$grp)
VC$grp <- gsub("Pass", "Row", VC$grp)


VC$grp <- factor(VC$grp, levels = c("g", "E", "g:E", "Range", "Row", "Rep:E", "Residual"))




p1 <-  ggplot(VC, aes(x=Trait , y=Percent/100)) +
  geom_col(mapping=aes(fill=grp), width= 0.9) +
  geom_point(aes(y = Repeatability, shape="Repeatability"), fill="white",  color="black",size=4, alpha=1) +
  geom_point(aes(y = R, shape="R-squared"), fill="white",  color="black",size=4, alpha=1) +
  #  scale_x_discrete("")+
  scale_y_continuous("Explained percent variation, repeatability and R-squared") +
  #  scale_fill_manual(values = c("#B15928","#A6CEE3","#6A3D9A","#CAB2D6","#FF7F00","#FDBF6F","#E31A1C","#FB9A99","#33A02C" ,"#B2DF8A","#1F78B4" ))+
  #scale_fill_brewer(palette = "Paired")+
  #  scale_colour_manual(values =c("black","black"))+
  scale_shape_manual(values=c(21,24))+
  #  scale_fill_viridis_d()+
  theme_bw() +
  scale_x_discrete("")+
  theme(legend.position= c(0.53,0.3),
        legend.background =  element_rect(fill = alpha("gray80", 0.7) , linetype="solid"),
        legend.key.width = unit(0.3, 'cm'),
        legend.key.height = unit(0.1, 'cm'))+
  #      legend.key.width = unit(0.3, 'cm'),
  #      legend.key.height = unit(0.3, 'cm'))+
  #      facet_grid(~Type, scales = "free", space = "free")+
  #  labs(title = "A-) Multispectral HTP platform")+
  guides(fill=guide_legend(title="Variance component", nrow=6, byrow=TRUE))
#  ggtitle("A) Variance component")


p1


p <- ggplot(Yield, aes(x= reorder_within(Env, -Yield, Tester), y=Yield,fill=Tester)) +
  geom_boxplot(position = position_dodge2(preserve = "single"), size=0.2, outlier.size = 0.3)+
  # stat_summary(fun.y=mean, geom="point", aes(fill=CRI))+
  # geom_point(aes(x= reorder_within(Env, -Mean, Tester), y=Mean,  group = interaction(Env,CRI,Tester) ))+
  facet_grid(~Tester, scales = "free_x")+
  scale_x_reordered("Environment") +
  scale_y_continuous("Grain yield (t/ha)")+
  theme_bw()+
  theme(legend.position = c(0.5,0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8, angle = 45, hjust = 1),
        legend.background = element_rect(fill="lightblue",
                                         size=0.5,
                                         linetype="solid"),
        #        legend.key.size = unit(0.2, 'cm'), #change legend key size
        #        legend.key.height = unit(0.5, 'cm'), #change legend key height
        #        legend.key.width = unit(0.5, 'cm'), #change legend key width
        #        legend.text = element_text(size=8)
  )+
  ggtitle("B) Genotypic effect of hybrid for grain yield (t/ha)")
p

library(gridExtra)

grid.arrange(p1,p,layout_matrix = rbind(c(1,2,2)))


jpeg("D:\\Fatma\\Stability Yield\\Var. comp and Blup of Yield.jpeg", width = 11,height =4.5,units = "in", res=500)

grid.arrange(p1,p,layout_matrix = rbind(c(1,2,2,2)))

dev.off()
.  