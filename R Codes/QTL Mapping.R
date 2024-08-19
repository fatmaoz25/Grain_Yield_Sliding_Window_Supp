setwd("D:\\control\\QTL results")

library(readr)
list_csv_files <- list.files(path = "D:\\control\\QTL results")
df <- readr::read_csv(list_csv_files, id = "file_name") %>% as.data.frame()
head(df)

df$file_name <- gsub("QTL.", "", df$file_name)
df$file_name <- gsub(".csv", "", df$file_name)
df <- tidyr::separate(df, file_name, into = c("Tester", "Env"), sep = "_")

head(df)

library(readxl)
library(dplyr)
df$Position  <- as.numeric(df$Position)

chromosome_data$Chr <- as.factor(chromosome_data$Chr)

chromosome_data <- df %>% group_by(chr,Tester) %>%  dplyr::summarise(start= min(Position),
                                                                     end=max(Position)) %>% as.data.frame()

chromosome_data$chr <- factor(chromosome_data$chr, levels = c(1,2,3,4,5,6,7,8,9,10))

qtl_data <- df[df$LOD>5,]
head(qtl_data)

qtl_data <- qtl_data %>%
  group_by(Tester, Env, chr) %>%
  mutate(Peak = max(LOD)) %>% as.data.frame()

qtl_data <- qtl_data[!duplicated(qtl_data$Peak),]

#write.csv(qtl_data, "qtl_data.csv")
qtl_data <- read.csv("qtl_data.csv")
qtl_data$Year <- as.factor(qtl_data$Year)
chromosome_data$chr <- as.factor(chromosome_data$chr)
head(qtl_data)
p <- ggplot() +
  geom_segment(data = chromosome_data, aes(x = start/1000000, xend = end/1000000, y = chr, yend = chr),
               size = 12, lineend = "round") +
  # geom_vline(data = df[df$TraitName=="Sen.85.Drought",],  aes(xintercept = Position), linetype = "solid", color = "white")+
  # geom_segment(data = df, aes(x = Position, xend = Position, y = Chromosome, yend = Chromosome), size = 2)+
  # geom_segment(data = qtl_data[!qtl_data$TraitName=="Slope", ], aes(x = LeftCI, xend = RightCI, y = Chromosome, yend = Chromosome, color=`PVE(%)`),alpha=1, size = 4) +
  # geom_text_repel(data = qtl_data, aes(x = (LeftCI + RightCI) / 2, y = Chromosome, label = label  ), vjust = -1,size=2,box.padding = 0.8) +
  geom_point(data = qtl_data, aes(x =  Position/1000000, y = chr, color = Env, shape=Year), alpha=0.6, size = 2,
             position = position_jitter(height = 0.3, width = 0)) +
  theme_classic() +
  #  scale_color_gradient("PVE(%)",  low="navy", high="red") +
  #  scale_color_brewer(palette = "Dark2")+
  facet_grid(~Tester,scales = "free_x")+
  labs(x = "Genomic position (Mbp)", y = "Chromosome")+
  guides( color = guide_legend(ncol = 2))
p

jpeg("D:\\control\\QTL results\\QTL.result.jpeg", width = 15,height =4.5,units = "in", res=500)
p
#grid.arrange(p, p1, p2, p3, ncol=2)
dev.off()