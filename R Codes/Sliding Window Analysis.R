rm(list=ls())
setwd("E:\\Fatma_weather\\Weather data")

#devtools::install_github('allogamous/EnvRtype',force=TRUE)
library(dplyr)
library(EnvRtype)
data <- data.table::fread("E:\\Fatma\\Stability Yield\\envRtype.txt") %>% as.data.frame()
head(data)

env   = data$env
lon   = data$lon
lat   = data$lat
start = data$start
end   = data$end125 ## This is 100 days
country = rep('USA1',length(lon))

env.data = EnvRtype::get_weather(env.id = env,
                                 lat = lat,
                                 lon = lon,
                                 start.day = start,
                                 end.day = end,
                                 country = country,
                                 parallel = TRUE)


df.clim = processWTH(env.data = env.data,Tbase1 = 8,Tbase2 = 45,Topt1 = 30,Topt2 = 37)
head(df.clim)
df.clim$PTR <- df.clim$N/df.clim$GDD
df.clim$PTT <- df.clim$GDD*df.clim$N
head(df.clim)
#write.csv(df.clim, "D:\\Fatma\\Stability Yield\\df.clim.csv")

df.clim <- read.csv("E:\\Fatma\\Stability Yield\\df.clim.csv")
#df.clim <- df.clim[df.clim$daysFromStart <= 60,]
#df.clim$PTR <- (df.clim$N)/(df.clim$GDD)

setwd("E:\\Fatma\\Stability Yield")
yield <- read.csv("Yield.wide.csv")
yield <-  subset(yield, select = -c(PHK76.Yield, PHP02.Yield,PHZ51.Yield))
yield <- cbind(yield$Inbred, stack(yield[,3:ncol(yield)]))

head(yield)
tail(yield)
library(dplyr)
yield <- yield %>% tidyr::separate(ind, into = c("Tester", "Env"), sep = "_")

names(yield) <- c("Inbred", "Yield", "Tester", "env")

yield <- yield %>% group_by(env,Tester) %>% mutate(Mean = mean(Yield)) %>% as.data.frame()

unique(yield$Tester)
#write.csv(unique(yield$env), "env.csv")


phk76 <- yield[yield$Tester=="PHK76",]
php02 <- yield[yield$Tester=="PHP02",]
phz51 <- yield[yield$Tester=="PHZ51",]


phk76 <- phk76[!duplicated(phk76$env),][,c(4,5)]
php02 <- php02[!duplicated(php02$env),][,c(4,5)]
phz51 <- phz51[!duplicated(phz51$env),][,c(4,5)]

head(phk76)
head(php02)
head(phz51)

head(df.clim)

# Assuming your data frame is named 'data' and has columns 'Environment', 'DL', 'GDD', 'Precipitation'
normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

df.clim <- df.clim %>%
  group_by(env) %>%
  mutate(
    DL_normalized = normalize(N),
    GDD_normalized = normalize(GDD),
    Precipitation_normalized = normalize(PRECTOT)
  ) %>%
  ungroup() %>% as.data.frame()

# View the normalized data
head(df.clim)


df.clim$new <-df.clim$DL_normalized/df.clim$GDD_normalized
library(ggplot2)
p<-ggplot(df.clim, aes(x=daysFromStart, y=new, group=env)) +
  geom_line(aes(color=env))
#  geom_point(aes(color=supp))
p



temperature_data <- reshape2::dcast(df.clim,  daysFromStart~env, value.var = "PTR")
head(temperature_data)
names(temperature_data)[1] <- "day"


# Load necessary library
library(dplyr)

# Assuming temperature_data is already defined as before
set.seed(123) # Ensure reproducibility for simulated data

# Define a function to calculate mean temperatures for all window sizes for a given environment
calculate_mean_temps_for_all_windows <- function(temperature_vector) {
  results <- list()
  total_days <- length(temperature_vector)
  
  for (window_size in 1:total_days) {
    mean_temps <- sapply(1:(total_days-window_size+1), function(start_day) {
      mean(temperature_vector[start_day:(start_day+window_size-1)])
    })
    results[[as.character(window_size)]] <- mean_temps
  }
  
  return(results)
}

# Apply the function to each environment column (excluding the 'day' column)
all_env_means <- list()

for (env in names(temperature_data)[-1]) { # Skip 'day' column
  env_data <- temperature_data[[env]]
  all_env_means[[env]] <- calculate_mean_temps_for_all_windows(env_data)
}

# Now all_env_means contains the mean temperatures for every window size for each environment


library(data.table)
library(zoo) # For rollapply, if not already loaded

# Initialize an empty data.table to store the results
results_dt <- data.table(environment = character(),
                         window_size = integer(),
                         start_day = integer(),
                         end_day = integer(),
                         mean_temperature = numeric())

# Environment names, assuming they were the column names in your original temperature_data
environment_names <- names(temperature_data)[-1] # Adjust if necessary

# Flatten the results and populate results_dt
for (i in seq_along(all_env_means)) {
  env_results <- all_env_means[[i]] # Results for the current environment
  env_name <- environment_names[i]
  
  for (window_size in seq_along(env_results)) {
    window_data <- env_results[[window_size]]
    if (is.null(window_data)) next # Skip if there's no data for this window size
    
    # Calculate start and end days for each window
    start_days <- seq_len(length(window_data))
    end_days <- start_days + window_size - 1
    
    # Append to results_dt
    results_dt <- rbindlist(list(results_dt, data.table(environment = env_name,
                                                        window_size = window_size,
                                                        start_day = start_days,
                                                        end_day = end_days,
                                                        mean_temperature = window_data)),
                            use.names = TRUE, fill = TRUE)
  }
}

results_dt <- results_dt %>% as.data.frame()

# View the results

head(results_dt)
names(results_dt)[1] <- "env"

######################
####### PHK76 ########
######################


phk76.cor <- dplyr::left_join(results_dt, phk76, by="env")
head(phk76.cor)
phk76.cor <- na.omit(phk76.cor)

phk76.cor <- phk76.cor %>% group_by(start_day,end_day) %>% mutate(cor=cor(mean_temperature,Mean,use = "complete.obs")) %>% as.data.frame()
str(phk76.cor)
library(ggplot2)
head(phk76.cor)

max_cor <- max(phk76.cor$cor)
min_cor <- min(phk76.cor$cor)

max_cor
min_cor

dff <- na.omit(phk76.cor[phk76.cor$cor==max_cor,])
env.index.phk76.GY <- dff

p <- ggplot(phk76.cor,aes(x= start_day ,y=end_day))+
  geom_tile(aes(fill=cor))+
  annotate("segment",x=0,xend=unique(dff$start_day), y=unique(dff$end_day), yend=unique(dff$end_day),colour="black",size=0.5,linetype="dashed")+
  annotate("segment",x=unique(dff$start_day),xend=unique(dff$start_day), y=0, yend=unique(dff$end_day),colour="black",size=0.5,linetype="dashed")+
  #annotate("segment", x = 0, xend = unique(dff$start_day), y = dff$end_day, yend = unique(dff$end_day),colour = "black")+
  geom_point(dff, mapping=aes(start_day, end_day),shape=1, size=4, color="black")+
  xlab("Begining of window (days)")+
  annotate("text",
           label = paste( "(", unique(dff$start_day),",",unique(dff$end_day)  ,")",   sep=""),
           x = dff$start_day+17, y = dff$end_day, colour = "black"
  )+
  #  geom_label(label = paste( "(", unique(dff$start_day),",",unique(dff$end_day)  ,")",   sep=""))+
  ylab("End of window (days)")+
  theme_bw()+
  scale_fill_distiller("cor (R)",palette = "Spectral",breaks = seq(from = -1, to = 1, by = 0.2))+
  #  scale_fill_continuous(colors = "viridis")+
  #  scale_fill_viridis_c()+
  # scale_fill_gradient2(low = "darkgreen", high = "brown",mid="white",midpoint=0,breaks = seq(from = -1, to = 1, by = 0.2))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour="black"))+
  #ggtitle(name_cor_plot[q])+
  theme(plot.title = element_text(hjust = 0.5),
        #  text = element_text(family = "sans"),
        legend.position = c(0.9,0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(12,"pt"))+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120))+
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120))
p

library(ggrepel)
library(ggpubr)

p1 <- ggplot(data= dff, aes(x=mean_temperature, y=Mean, group=1)) +
  geom_smooth(method = lm, se = T)+
  geom_point(color="red", size=3)+
  geom_text_repel(aes(label = env),size=3 )+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_cor(method = "pearson")+
  scale_y_continuous("Environmental mean of yield (t/ha)")+
  scale_x_continuous( paste("PTR"," ", "(", unique(dff$start_day),",",unique(dff$end_day)  ,")",   sep=""))
p1

library(gridExtra)
grid.arrange(p,p1,layout_matrix = rbind(c(1,2)))

jpeg("D:\\Fatma\\Stability Yield\\heatmap.phk76.Yield.jpeg", width = 8,height =4,units = "in", res=500)
grid.arrange(p,p1,layout_matrix = rbind(c(1,2)))
dev.off()



######################
####### PHP02 ########
######################

php02.cor <- dplyr::left_join(results_dt, php02, by="env")
head(php02.cor)
php02.cor <- na.omit(php02.cor)

php02.cor <- php02.cor %>% group_by(start_day,end_day) %>% mutate(cor=cor(mean_temperature,Mean,use = "complete.obs")) %>% as.data.frame()
str(php02.cor)
library(ggplot2)
head(php02.cor)

max_cor <- max(php02.cor$cor)
min_cor <- min(php02.cor$cor)

max_cor
min_cor

dff <- na.omit(php02.cor[php02.cor$cor==max_cor,])
env.index.php02.GY <- dff

p <- ggplot(php02.cor,aes(x= start_day ,y=end_day))+
  geom_tile(aes(fill=cor))+
  annotate("segment",x=0,xend=unique(dff$start_day), y=unique(dff$end_day), yend=unique(dff$end_day),colour="black",size=0.5,linetype="dashed")+
  annotate("segment",x=unique(dff$start_day),xend=unique(dff$start_day), y=0, yend=unique(dff$end_day),colour="black",size=0.5,linetype="dashed")+
  #annotate("segment", x = 0, xend = unique(dff$start_day), y = dff$end_day, yend = unique(dff$end_day),colour = "black")+
  geom_point(dff, mapping=aes(start_day, end_day),shape=1, size=4, color="black")+
  xlab("Begining of window (days)")+
  annotate("text",
           label = paste( "(", unique(dff$start_day),",",unique(dff$end_day)  ,")",   sep=""),
           x = dff$start_day+17, y = dff$end_day, colour = "black"
  )+
  #  geom_label(label = paste( "(", unique(dff$start_day),",",unique(dff$end_day)  ,")",   sep=""))+
  ylab("End of window (days)")+
  theme_bw()+
  scale_fill_distiller("cor (R)",palette = "Spectral",breaks = seq(from = -1, to = 1, by = 0.2))+
  #  scale_fill_continuous(colors = "viridis")+
  #  scale_fill_viridis_c()+
  # scale_fill_gradient2(low = "darkgreen", high = "brown",mid="white",midpoint=0,breaks = seq(from = -1, to = 1, by = 0.2))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour="black"))+
  #ggtitle(name_cor_plot[q])+
  theme(plot.title = element_text(hjust = 0.5),
        #  text = element_text(family = "sans"),
        legend.position = c(0.9,0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(12,"pt"))+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120))+
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120))
p

library(ggrepel)
library(ggpubr)

p1 <- ggplot(data= dff, aes(x=mean_temperature, y=Mean, group=1)) +
  geom_smooth(method = lm, se = T)+
  geom_point(color="red", size=3)+
  geom_text_repel(aes(label = env),size=3 )+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_cor(method = "pearson")+
  scale_y_continuous("Environmental mean of yield (t/ha)")+
  scale_x_continuous( paste("PTR"," ", "(", unique(dff$start_day),",",unique(dff$end_day)  ,")",   sep=""))
p1
library(gridExtra)
grid.arrange(p,p1,layout_matrix = rbind(c(1,2)))

jpeg("D:\\Fatma\\Stability Yield\\heatmap.php02.Yield.jpeg", width = 8,height =4,units = "in", res=500)
grid.arrange(p,p1,layout_matrix = rbind(c(1,2)))
dev.off()


######################
####### PHK51 ########
######################

phz51.cor <- dplyr::left_join(results_dt, phz51, by="env")
head(phz51.cor)
phz51.cor <- na.omit(phz51.cor)

phz51.cor <- phz51.cor %>% group_by(start_day,end_day) %>% mutate(cor=cor(mean_temperature,Mean,use = "complete.obs")) %>% as.data.frame()
str(phz51.cor)
library(ggplot2)
head(phz51.cor)

max_cor <- max(phz51.cor$cor)
min_cor <- min(phz51.cor$cor)

max_cor
min_cor

dff <- na.omit(phz51.cor[phz51.cor$cor==max_cor,])
env.index.phz51.GY <- dff

p <- ggplot(phz51.cor,aes(x= start_day ,y=end_day))+
  geom_tile(aes(fill=cor))+
  annotate("segment",x=0,xend=unique(dff$start_day), y=unique(dff$end_day), yend=unique(dff$end_day),colour="black",size=0.5,linetype="dashed")+
  annotate("segment",x=unique(dff$start_day),xend=unique(dff$start_day), y=0, yend=unique(dff$end_day),colour="black",size=0.5,linetype="dashed")+
  #annotate("segment", x = 0, xend = unique(dff$start_day), y = dff$end_day, yend = unique(dff$end_day),colour = "black")+
  geom_point(dff, mapping=aes(start_day, end_day),shape=1, size=4, color="black")+
  xlab("Begining of window (days)")+
  annotate("text",
           label = paste( "(", unique(dff$start_day),",",unique(dff$end_day)  ,")",   sep=""),
           x = dff$start_day+17, y = dff$end_day, colour = "black"
  )+
  #  geom_label(label = paste( "(", unique(dff$start_day),",",unique(dff$end_day)  ,")",   sep=""))+
  ylab("End of window (days)")+
  theme_bw()+
  scale_fill_distiller("cor (R)",palette = "Spectral",breaks = seq(from = -1, to = 1, by = 0.2))+
  #  scale_fill_continuous(colors = "viridis")+
  #  scale_fill_viridis_c()+
  # scale_fill_gradient2(low = "darkgreen", high = "brown",mid="white",midpoint=0,breaks = seq(from = -1, to = 1, by = 0.2))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour="black"))+
  #ggtitle(name_cor_plot[q])+
  theme(plot.title = element_text(hjust = 0.5),
        #  text = element_text(family = "sans"),
        legend.position = c(0.9,0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(12,"pt"))+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120))+
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100,110,120))
p

library(ggrepel)
library(ggpubr)

p1 <- ggplot(data= dff, aes(x=mean_temperature, y=Mean, group=1)) +
  geom_smooth(method = lm, se = T)+
  geom_point(color="red", size=3)+
  geom_text_repel(aes(label = env),size=3 )+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_cor(method = "pearson")+
  scale_y_continuous("Environmental mean of yield (t/ha)")+
  scale_x_continuous( paste("PTR"," ", "(", unique(dff$start_day),",",unique(dff$end_day)  ,")",   sep=""))
p1

library(gridExtra)
grid.arrange(p,p1,layout_matrix = rbind(c(1,2)))

jpeg("D:\\Fatma\\Stability Yield\\heatmap.phz51.Yield.jpeg", width = 8,height =4,units = "in", res=500)
grid.arrange(p,p1,layout_matrix = rbind(c(1,2)))
dev.off()
