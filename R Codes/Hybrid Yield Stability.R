
#############################################
############## Stability Yield ##############
#############################################

GY.wide <- read.csv("Yield.wide.csv")
GY.phk76 <- GY.wide[,c(2,3,4,5,6,7,8,9,10,11,12,13)]
GY.phk76 <- cbind(GY.phk76[,1], stack(GY.phk76[,2:ncol(GY.phk76)]))
names(GY.phk76) <- c("Inbred", "GY", "env")
GY.phk76 <- GY.phk76 %>% tidyr::separate(env, into = c("Tester", "env"), sep = "_")

head(GY.phk76)
head(env.index.phk76.GY)

GY.phk76.stability <- dplyr::left_join(GY.phk76, env.index.phk76.GY, by="env")
head(GY.phk76.stability)



GY.wide <- read.csv("Yield.wide.csv")
head(GY.wide)
GY.php02 <- GY.wide[,c(2,15:28)]
head(GY.php02)
GY.php02 <- cbind(GY.php02[,1], stack(GY.php02[,2:ncol(GY.php02)]))
names(GY.php02) <- c("Inbred", "GY", "env")
GY.php02 <- GY.php02 %>% tidyr::separate(env, into = c("Tester", "env"), sep = "_")

head(GY.php02)
head(env.index.php02.GY)

GY.php02.stability <- dplyr::left_join(GY.php02, env.index.php02.GY, by="env")
head(GY.php02.stability)



GY.wide <- read.csv("Yield.wide.csv")
head(GY.wide)
GY.phz51 <- GY.wide[,c(2,30:49)]
head(GY.phz51)
GY.phz51 <- cbind(GY.phz51[,1], stack(GY.phz51[,2:ncol(GY.phz51)]))
names(GY.phz51) <- c("Inbred", "GY", "env")
GY.phz51 <- GY.phz51 %>% tidyr::separate(env, into = c("Tester", "env"), sep = "_")

head(GY.phz51)
head(env.index.phz51.GY)

GY.phz51.stability <- dplyr::left_join(GY.phz51, env.index.phz51.GY, by="env")
head(GY.phz51.stability)

GY.stability <- rbind(GY.phk76.stability,GY.php02.stability,GY.phz51.stability)
head(GY.stability)


GY.stability <- GY.stability %>% group_by(env,Tester) %>% mutate(ScaledGY = scale(GY, center = TRUE, scale = TRUE)) %>% as.data.frame()

GY.stability <- GY.stability %>% group_by(Inbred,Tester) %>% dplyr::mutate(Slope = coef(lm(ScaledGY ~ mean_temperature))[2],
                                                                           Intercept = coef(lm(ScaledGY ~ mean_temperature))[1],
                                                                           Predicted.ScaledGY = predict(lm(ScaledGY ~ mean_temperature))) %>% as.data.frame()

GY.stability <- as.data.frame(lapply(GY.stability, unlist))
GY.stability$Slope <- as.numeric(GY.stability$Slope)


result <- GY.stability %>%  group_by(env,Tester) %>% dplyr::mutate(Meann = mean(Predicted.ScaledGY)) %>% as.data.frame()
result <- result %>% distinct(env,Tester, .keep_all = TRUE)
head(result)



p12<- ggplot(GY.stability, aes(x=mean_temperature, y=Predicted.ScaledGY, group=Inbred, color=Slope)) +
  geom_line(alpha=0.4) +
  # geom_smooth(se=F)+
  geom_point(alpha=0.3)+
  # geom_text( result, mapping=aes(label = Env, x=Mean, y=Meann-3), color='black', angle=90, size=1, alpha=0.8)+
  geom_point(result, mapping=aes(x=mean_temperature, y=2.48 ),color="black" ,shape=24,alpha=0.5   )+
  geom_text_repel(result,mapping=aes(x=mean_temperature, y=2.48, label=env ),color="black",size=2,
                  force_pull   = 0, # do not pull toward data points
                  nudge_y      = 0.05,
                  direction    = "x",
                  angle        = 90,
                  hjust        = 0,
                  segment.size = 0.2,
                  max.iter = 1e4, max.time = 1
  )+
  #scale_color_brewer(palette = "Dark2")+
  # geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2,
  #                position=position_dodge(0.05))+
  scale_y_continuous(name="Predicted GY")+
  theme(legend.position="right") +
  theme_bw()+
  facet_grid(~Tester,scales ="free_x")+
  scale_x_continuous(name = "Environmental index"  )+
  theme(#legend.background = element_blank(),
    legend.background =  element_rect(fill = alpha("gray80", 0.3) , size=0.2, linetype="solid"),
    legend.key.height = unit(0.5, 'cm'),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.key.size = unit(12,"pt"),
    axis.text.y = element_text(angle = 90, vjust=0.8,  hjust=0.5),
    legend.position   = c(0.75, 0.25))+
  scale_colour_gradient2("Slope", low = "navy", mid = "cyan", high = "red", midpoint = mean(GY.stability$Slope, na.rm=T) )
#  guides(color = guide_legend(override.aes = list(size = 1,alpha = 1) ) )+
#  ggtitle("A) Slope of the temporal plant height across environments")

p12



jpeg( "D:\\Fatma\\Stability Yield\\Stability of GY by environmental index.phk76.php02,phz51.jpeg",   width = 8,height =4,units = "in", res=500)
p12
dev.off()


head(GY.stability)

p <- ggplot(GY.stability, aes(x=ScaledGY, y=Predicted.ScaledGY)) +
  geom_point()+
  geom_smooth(method = lm, se = T)+
  stat_cor(method = "pearson", color="red")+
  facet_wrap(Tester~env,scales = "free")+
  theme(#axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    legend.key.height = unit(0.5, 'cm'),
    strip.text = element_text(size = 8, margin = margin(t = 1, b = 1, unit = "pt")))+
  scale_x_continuous("GY scaled")+
  scale_y_continuous("GY scaled predicted")


p

jpeg( "D:\\Fatma\\Stability Yield\\Correlation GY.jpeg",   width = 13,height =14,units = "in", res=500)
p
dev.off()


GY.stability.wide <- GY.stability %>% dplyr::distinct(Inbred, Tester, .keep_all = TRUE)
GY.stability.wide <-data.table::dcast(setDT(GY.stability.wide), Inbred ~ Tester, value.var = c("Slope", "Intercept")) %>% as.data.frame()

write.csv(GY.stability.wide, "GY.stability.wide.Slope.Intercept.csv")