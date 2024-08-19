#####################################################
################ Pheno for PHK76 ####################
#####################################################
rm(list = ls())
pheno <- read.csv("D:\\FATMA PREDICTION\\Yield.wide.csv")
pheno <- pheno[,-c(1,14,29,50)]
pheno <- reshape2::melt(pheno, id.var = c('Inbred'), variable.name = 'Tester_Env')
pheno <- tidyr::separate(pheno, col=Tester_Env, into=c('Tester', 'Env'), sep='_')


pheno <- pheno[pheno$Tester=="PHK76",]
pheno <- pheno[,-2]
names(pheno)[1:3] <- c("inbred", "env", "yield")
head(pheno)

#########################
####### Env data ########
#########################

library(dplyr)
wht <- data.table::fread("D:\\FATMA PREDICTION\\df.clim.csv") %>% as.data.frame()
wht <- wht[wht$env %in% unique(pheno$env), ] #### subset teh envs related to PHK76
wht <- wht[wht$daysFromStart %in% c(100,101), ] #### Subset the days of environmental index
wht <- wht %>% group_by(env) %>% dplyr::summarize(PTR = mean(PTR) ) %>% as.data.frame()
pheno <- left_join(pheno, wht, by="env")
head(pheno)

#######################
####### GRM ###########
#######################

geno <- data.table::fread("D:\\FATMA PREDICTION\\cm.mbp.ss.100k.sites.SSpopulation_geno.csv") %>% as.data.frame()
row.names(geno) <- geno$lines
geno <- geno[,-1]

geno <- rrBLUP::A.mat(geno, min.MAF = 0.05, impute.method = "EM",  tol = 0.02,n.core = 10,shrink = FALSE,  return.imputed = T)
ZG <- geno$A

CV1 <- list()
CV2 <- list()
CV0 <- list()
CV00 <- list()

library(BGLR)
library(parallel)
library(doParallel)
library(MASS)

env <- unique(pheno$env)
index =1

numCores <- detectCores()
numCores
registerDoParallel(numCores-3)

library(BGLR)

env <- unique(pheno$env)

#foreach(j = 1:200, .packages = c("BGLR", "dplyr")) %dopar%

for (j in 1:100)
{
  for (i in 1:length(env) )
  {
    untested.env <- env[i]
    
    df <- pheno
    df$Y2 <- df$yield
    hybrids <- unique(df$inbred)
    
    train_hybrids <- sample(hybrids, size = length(hybrids) * 0.7)
    test_hybrids <- hybrids[!hybrids %in% train_hybrids]
    
    df$Y2[df$inbred%in%test_hybrids] <- NA
    df$Y2[pheno$env == untested.env] <- NA
    
    train <- df[df$inbred %in% train_hybrids & !df$env == untested.env,]
    test <- df[df$inbred %in% test_hybrids & df$env == untested.env,]
    
    out.linear <- train %>% group_by(inbred) %>% dplyr::summarise(Slope = coef(lm(Y2 ~ PTR))[2],
                                                                  Intercept = coef(lm(Y2 ~ PTR))[1]
                                                                  #Predicted.yield = predict(lm(Y2 ~ PTR))
    ) %>% as.data.frame()
    
    df2 <- dplyr::left_join( df,out.linear, by="inbred")
    df2 <- df2[!duplicated(df2$inbred),]
    
    y_slope <- as.numeric(df2$Slope)
    fit.slope <-  BGLR(y = y_slope, ETA = list(list(K = ZG, model = "RKHS")), nIter = 10000, burnIn = 1000, verbose = F)
    
    y_Intercept <- as.numeric(df2$Intercept)
    fit.intercept <-  BGLR(y = y_Intercept, ETA = list(list(K = ZG, model = "RKHS")), nIter = 10000, burnIn = 1000, verbose = F)
    
    df2$predicted.slope <- fit.slope$yHat
    df2$predicted.intercept <- fit.intercept$yHat
    
    df3 <- dplyr::left_join(df,df2[,c('inbred','predicted.slope','predicted.intercept')],by = "inbred")
    
    df3$predicted.yield <- (df3$predicted.slope*df3$PTR) + df3$predicted.intercept
    head(df3)
    
    df3.CV1 <- df3[df3$inbred %in% train_hybrids & !df3$env == untested.env,]
    df3.CV2 <- df3[df3$inbred %in% test_hybrids & !df3$env == untested.env,]
    df3.CV0 <- df3[df3$inbred %in% train_hybrids & df3$env == untested.env,]
    df3.CV00 <- df3[df3$inbred %in% test_hybrids & df3$env == untested.env,]
    
    
    CV1[[index]] <- df3.CV1 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    CV2[[index]] <- df3.CV2 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    CV0[[index]] <- df3.CV0 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    CV00[[index]]<- df3.CV00 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    
    index = index + 1
    
  }
  
}

CV1.phk76 <- plyr::ldply(CV1, data.frame)
CV2.phk76 <- plyr::ldply(CV2, data.frame)
CV0.phk76 <- plyr::ldply(CV0, data.frame)
CV00.phk76 <- plyr::ldply(CV00, data.frame)

write.csv(CV1.phk76, "D:\\FATMA PREDICTION\\CV1.phk76.csv")
write.csv(CV2.phk76, "D:\\FATMA PREDICTION\\CV2.phk76.csv")
write.csv(CV0.phk76, "D:\\FATMA PREDICTION\\CV0.phk76.csv")
write.csv(CV00.phk76, "D:\\FATMA PREDICTION\\CV00.phk76.csv")





#####################################################
################ Pheno for PHP02 ####################
#####################################################
rm(list = ls())
pheno <- read.csv("D:\\FATMA PREDICTION\\Yield.wide.csv")
pheno <- pheno[,-c(1,14,29,50)]
pheno <- reshape2::melt(pheno, id.var = c('Inbred'), variable.name = 'Tester_Env')
pheno <- tidyr::separate(pheno, col=Tester_Env, into=c('Tester', 'Env'), sep='_')


pheno <- pheno[pheno$Tester=="PHP02",]
pheno <- pheno[,-2]
names(pheno)[1:3] <- c("inbred", "env", "yield")
head(pheno)

#########################
####### Env data ########
#########################

library(dplyr)
wht <- data.table::fread("D:\\FATMA PREDICTION\\df.clim.csv") %>% as.data.frame()
wht <- wht[wht$env %in% unique(pheno$env), ] #### subset teh envs related to PHK76
wht <- wht[wht$daysFromStart %in% c(67,68), ] #### Subset the days of environmental index
wht <- wht %>% group_by(env) %>% dplyr::summarize(PTR = mean(PTR) ) %>% as.data.frame()
pheno <- left_join(pheno, wht, by="env")
head(pheno)

#######################
####### GRM ###########
#######################

geno <- data.table::fread("D:\\FATMA PREDICTION\\cm.mbp.ss.100k.sites.SSpopulation_geno.csv") %>% as.data.frame()
row.names(geno) <- geno$lines
geno <- geno[,-1]

geno <- rrBLUP::A.mat(geno, min.MAF = 0.05, impute.method = "EM",  tol = 0.02,n.core = 10,shrink = FALSE,  return.imputed = T)
ZG <- geno$A

CV1 <- list()
CV2 <- list()
CV0 <- list()
CV00 <- list()

library(BGLR)
library(parallel)
library(doParallel)
library(MASS)

env <- unique(pheno$env)
index =1

#numCores <- detectCores()
#numCores
#registerDoParallel(numCores-3)

library(BGLR)

env <- unique(pheno$env)

#foreach(j = 1:200, .packages = c("BGLR", "dplyr")) %dopar%

for (j in 1:100)
{
  for (i in 1:length(env) )
  {
    untested.env <- env[i]
    
    df <- pheno
    df$Y2 <- df$yield
    hybrids <- unique(df$inbred)
    
    train_hybrids <- sample(hybrids, size = length(hybrids) * 0.7)
    test_hybrids <- hybrids[!hybrids %in% train_hybrids]
    
    df$Y2[df$inbred%in%test_hybrids] <- NA
    df$Y2[pheno$env == untested.env] <- NA
    
    train <- df[df$inbred %in% train_hybrids & !df$env == untested.env,]
    test <- df[df$inbred %in% test_hybrids & df$env == untested.env,]
    
    out.linear <- train %>% group_by(inbred) %>% dplyr::summarise(Slope = coef(lm(Y2 ~ PTR))[2],
                                                                  Intercept = coef(lm(Y2 ~ PTR))[1]
                                                                  #Predicted.yield = predict(lm(Y2 ~ PTR))
    ) %>% as.data.frame()
    
    df2 <- dplyr::left_join( df,out.linear, by="inbred")
    df2 <- df2[!duplicated(df2$inbred),]
    
    y_slope <- as.numeric(df2$Slope)
    fit.slope <-  BGLR(y = y_slope, ETA = list(list(K = ZG, model = "RKHS")), nIter = 10000, burnIn = 1000, verbose = F)
    
    y_Intercept <- as.numeric(df2$Intercept)
    fit.intercept <-  BGLR(y = y_Intercept, ETA = list(list(K = ZG, model = "RKHS")), nIter = 10000, burnIn = 1000, verbose = F)
    
    df2$predicted.slope <- fit.slope$yHat
    df2$predicted.intercept <- fit.intercept$yHat
    
    df3 <- dplyr::left_join(df,df2[,c('inbred','predicted.slope','predicted.intercept')],by = "inbred")
    
    df3$predicted.yield <- (df3$predicted.slope*df3$PTR) + df3$predicted.intercept
    head(df3)
    
    df3.CV1 <- df3[df3$inbred %in% train_hybrids & !df3$env == untested.env,]
    df3.CV2 <- df3[df3$inbred %in% test_hybrids & !df3$env == untested.env,]
    df3.CV0 <- df3[df3$inbred %in% train_hybrids & df3$env == untested.env,]
    df3.CV00 <- df3[df3$inbred %in% test_hybrids & df3$env == untested.env,]
    
    
    CV1[[index]] <- df3.CV1 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    CV2[[index]] <- df3.CV2 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    CV0[[index]] <- df3.CV0 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    CV00[[index]]<- df3.CV00 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    
    index = index + 1
    
  }
  
}

CV1.php02 <- plyr::ldply(CV1, data.frame)
CV2.php02 <- plyr::ldply(CV2, data.frame)
CV0.php02 <- plyr::ldply(CV0, data.frame)
CV00.php02 <- plyr::ldply(CV00, data.frame)

write.csv(CV1.php02, "D:\\FATMA PREDICTION\\CV1.php02.csv")
write.csv(CV2.php02, "D:\\FATMA PREDICTION\\CV2.php02.csv")
write.csv(CV0.php02, "D:\\FATMA PREDICTION\\CV0.php02.csv")
write.csv(CV00.php02, "D:\\FATMA PREDICTION\\CV00.php02.csv")





#####################################################
################ Pheno for PHZ51 ####################
#####################################################
rm(list = ls())
pheno <- read.csv("D:\\FATMA PREDICTION\\Yield.wide.csv")
pheno <- pheno[,-c(1,14,29,50)]
pheno <- reshape2::melt(pheno, id.var = c('Inbred'), variable.name = 'Tester_Env')
pheno <- tidyr::separate(pheno, col=Tester_Env, into=c('Tester', 'Env'), sep='_')


pheno <- pheno[pheno$Tester=="PHZ51",]
pheno <- pheno[,-2]
names(pheno)[1:3] <- c("inbred", "env", "yield")
head(pheno)

#########################
####### Env data ########
#########################

library(dplyr)
wht <- data.table::fread("D:\\FATMA PREDICTION\\df.clim.csv") %>% as.data.frame()
wht <- wht[wht$env %in% unique(pheno$env), ] #### subset teh envs related to PHK76
wht <- wht[wht$daysFromStart %in% c(46,90), ] #### Subset the days of environmental index
wht <- wht %>% group_by(env) %>% dplyr::summarize(PTR = mean(PTR) ) %>% as.data.frame()
pheno <- left_join(pheno, wht, by="env")
head(pheno)

#######################
####### GRM ###########
#######################

geno <- data.table::fread("D:\\FATMA PREDICTION\\cm.mbp.ss.100k.sites.SSpopulation_geno.csv") %>% as.data.frame()
row.names(geno) <- geno$lines
geno <- geno[,-1]

geno <- rrBLUP::A.mat(geno, min.MAF = 0.05, impute.method = "EM",  tol = 0.02,n.core = 10,shrink = FALSE,  return.imputed = T)
ZG <- geno$A

CV1 <- list()
CV2 <- list()
CV0 <- list()
CV00 <- list()

library(BGLR)
library(parallel)
library(doParallel)
library(MASS)

env <- unique(pheno$env)
index =1

#numCores <- detectCores()
#numCores
#registerDoParallel(numCores-3)

library(BGLR)

env <- unique(pheno$env)

#foreach(j = 1:200, .packages = c("BGLR", "dplyr")) %dopar%

for (j in 1:100)
{
  for (i in 1:length(env) )
  {
    untested.env <- env[i]
    
    df <- pheno
    df$Y2 <- df$yield
    hybrids <- unique(df$inbred)
    
    train_hybrids <- sample(hybrids, size = length(hybrids) * 0.7)
    test_hybrids <- hybrids[!hybrids %in% train_hybrids]
    
    df$Y2[df$inbred%in%test_hybrids] <- NA
    df$Y2[pheno$env == untested.env] <- NA
    
    train <- df[df$inbred %in% train_hybrids & !df$env == untested.env,]
    test <- df[df$inbred %in% test_hybrids & df$env == untested.env,]
    
    out.linear <- train %>% group_by(inbred) %>% dplyr::summarise(Slope = coef(lm(Y2 ~ PTR))[2],
                                                                  Intercept = coef(lm(Y2 ~ PTR))[1]
                                                                  #Predicted.yield = predict(lm(Y2 ~ PTR))
    ) %>% as.data.frame()
    
    df2 <- dplyr::left_join( df,out.linear, by="inbred")
    df2 <- df2[!duplicated(df2$inbred),]
    
    y_slope <- as.numeric(df2$Slope)
    fit.slope <-  BGLR(y = y_slope, ETA = list(list(K = ZG, model = "RKHS")), nIter = 10000, burnIn = 1000, verbose = F)
    
    y_Intercept <- as.numeric(df2$Intercept)
    fit.intercept <-  BGLR(y = y_Intercept, ETA = list(list(K = ZG, model = "RKHS")), nIter = 10000, burnIn = 1000, verbose = F)
    
    df2$predicted.slope <- fit.slope$yHat
    df2$predicted.intercept <- fit.intercept$yHat
    
    df3 <- dplyr::left_join(df,df2[,c('inbred','predicted.slope','predicted.intercept')],by = "inbred")
    
    df3$predicted.yield <- (df3$predicted.slope*df3$PTR) + df3$predicted.intercept
    head(df3)
    
    df3.CV1 <- df3[df3$inbred %in% train_hybrids & !df3$env == untested.env,]
    df3.CV2 <- df3[df3$inbred %in% test_hybrids & !df3$env == untested.env,]
    df3.CV0 <- df3[df3$inbred %in% train_hybrids & df3$env == untested.env,]
    df3.CV00 <- df3[df3$inbred %in% test_hybrids & df3$env == untested.env,]
    
    
    CV1[[index]] <- df3.CV1 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    CV2[[index]] <- df3.CV2 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    CV0[[index]] <- df3.CV0 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    CV00[[index]]<- df3.CV00 %>% group_by(env) %>% dplyr::summarize(cor=cor(yield, predicted.yield,use = "complete.obs")) %>% as.data.frame()
    
    index = index + 1
    
  }
  
}

CV1.phz51 <- plyr::ldply(CV1, data.frame)
CV2.phz51 <- plyr::ldply(CV2, data.frame)
CV0.phz51 <- plyr::ldply(CV0, data.frame)
CV00.phz51 <- plyr::ldply(CV00, data.frame)

write.csv(CV1.phz51, "D:\\FATMA PREDICTION\\CV1.phz51.csv")
write.csv(CV2.phz51, "D:\\FATMA PREDICTION\\CV2.phz51.csv")
write.csv(CV0.phz51, "D:\\FATMA PREDICTION\\CV0.phz51.csv")
write.csv(CV00.phz51, "D:\\FATMA PREDICTION\\CV00.phz51.csv")