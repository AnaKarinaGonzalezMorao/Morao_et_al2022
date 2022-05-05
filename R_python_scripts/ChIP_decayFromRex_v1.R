library(dplyr)
library(ggplot2)
library(slider)
library(reshape2)
library(ggpubr)
library(IRanges)
options(scipen = 999)

## for plotting asthetics
addUnits <- function(n) {
  labels <- ifelse(n < 1000, n,  # less than thousands
                   ifelse(n < 1e6, paste0(round(n/1e3), 'kb'),  # in thousands
                          ifelse(n < 1e9, paste0(round(n/1e6), 'Mb'),  # in millions
                                 ifelse(n < 1e12, paste0(round(n/1e9), 'B'), # in billions
                                        ifelse(n < 1e15, paste0(round(n/1e12), 'T'), # in trillions
                                               'too big!'
                                        )))))
  return(labels)
}

## 'TSS' does not need to be TSS, 
# any bed3-column + n number of value columns additional columns
# second 'rex' is bed3, to which each of 'TSS' rows will be assigned to
# ie. each row of 'TSS' is assigned to single closest row of 'rex'
assignTSStoNearestRex<- function(TSS,rex) {
  TSS_mid <- (TSS[,2]+TSS[,3])/2
  rex_mid <- (rex[,2]+rex[,3])/2
  TSS_rex_pair <- list()
  for (i in 1:nrow(TSS)){
    nearest_rex_i <- rex[abs(TSS_mid[i]-rex_mid)==min(abs(TSS_mid[i]-rex_mid)),]
    TSS_rex_pair[[i]] <- cbind(TSS[i,],nearest_rex_i)
  }
  do.call("rbind",TSS_rex_pair) %>% return
}

## generates submatrix, this is used to compute moving average on assignTSStoNearestRex output
Generate_submatrix <-  function(chrInfo,windowSize,stepSize) {
  output <- list()
  for (i in 1:length(chrInfo)) {
    nSteps <- (floor(as.numeric(chrInfo[[i]][2]))-windowSize)/stepSize
    V2=c(stepSize*(0:nSteps))
    V3=c(stepSize*(0:nSteps))+windowSize
    chr_i <- data.frame(V1=rep(chrInfo[[i]][1],nSteps+1),V2,V3)
    output[[chrInfo[[i]][1]]] <- chr_i
  }
  return(output)
}
#ce10 <- list(c("I",15072423),c("II",15279345),c("III",13783700),
#             c("IV",17493793),c("V",20924149),c("X",17718866))

#two column plus: 1st is position, second plus is value to average
MovingAverage <- function(twocolumn,windowSize,stepSize) {
  chrinfo <- list(c("NA",max(twocolumn[,1])))
  submat <- Generate_submatrix(chrinfo,windowSize,stepSize)$"NA"
  MovingAvg <- list()
  for (i in 1:nrow(submat)) {
    MovingAvg[[i]]<- data.frame(twocolumn[(twocolumn[,1]>submat[i,2])&(twocolumn[,1]<submat[i,3]),-1]) %>% colMeans
  }
  data.frame(submat[,c(2:3)],do.call("rbind",MovingAvg)) %>% return
}

# same as moving average instead computes confidence interval 95%
MovingCI <- function(twocolumn,windowSize,stepSize,CI_12) {
  chrinfo <- list(c("NA",max(twocolumn[,1])))
  submat <- Generate_submatrix(chrinfo,windowSize,stepSize)$"NA"
  MovingCI <- list()
  for (i in 1:nrow(submat)) {
    submat_i <- data.frame(twocolumn[(twocolumn[,1]>submat[i,2])&(twocolumn[,1]<submat[i,3]),-1])
    MovingCI[[i]]<-sapply(submat_i,function(x) t.test(x)$"conf.int"[CI_12])
  }
  data.frame(submat[,c(2:3)],do.call("rbind",MovingCI)) %>% return
}

# combined above functions into one big one
ChIPdecay_format <- function(rex,bin_bedgraph,windowSize,stepSize) {
  bin_X <- bin_bedgraph[bin_bedgraph[,1]=="chrX",]
  NumOfChIP <- ncol(bin_X)-3
  paired <- assignTSStoNearestRex(bin_X,rex[,c(1:3)])
  bin_mid <- (paired[,2]+paired[,3])/2
  rex_mid <- (paired[,sum(5,NumOfChIP)]+paired[,sum(6,NumOfChIP)])/2
  rexToBin <- abs(bin_mid-rex_mid)
  ChIPdecay <- data.frame(rexToBin,paired[,(4):(3+NumOfChIP)])
  ChIPdecay_cc <- ChIPdecay[complete.cases(ChIPdecay),]
  ChIPdecay_cc_sort <- ChIPdecay_cc[order(ChIPdecay_cc$rexToBin),]
  MovingAverage(ChIPdecay_cc_sort,windowSize,stepSize) %>% return
}


ChIPdecay_format_CI <- function(rex,bin_bedgraph,windowSize,stepSize,CI_12) {
  bin_X <- bin_bedgraph[bin_bedgraph[,1]=="chrX",]
  NumOfChIP <- ncol(bin_X)-3
  paired <- assignTSStoNearestRex(bin_X,rex[,c(1:3)])
  bin_mid <- (paired[,2]+paired[,3])/2
  rex_mid <- (paired[,sum(5,NumOfChIP)]+paired[,sum(6,NumOfChIP)])/2
  rexToBin <- abs(bin_mid-rex_mid)
  ChIPdecay <- data.frame(rexToBin,paired[,(4):(3+NumOfChIP)])
  ChIPdecay_cc <- ChIPdecay[complete.cases(ChIPdecay),]
  ChIPdecay_cc_sort <- ChIPdecay_cc[order(ChIPdecay_cc$rexToBin),]
  MovingCI(ChIPdecay_cc_sort,windowSize,stepSize,CI_12) %>% return
}

# sum chrX = 1
applyUnityX <- function(multi.tab) {
  M <- multi.tab[complete.cases(multi.tab),]
  M.x <- M[M[,1]=="chrX",]
  M.z <- sapply(M.x[,-c(1:3)],function(x) x/sum(x))
  cbind(M.x[,c(1:3)],M.z) %>% return
}
############################################# end of functions #################################
#rm(list = ls())

rex <- read.csv("C:/Users/kimj5/Desktop/annotation/sarah_strong_annot.bed",header=F,sep="\t",quote="'",check.names=F)
# bed file format: three column, chr, start, end

################################# main: strains with auxin or etoposide  ######################
sample_orders <- c("AKM189_AKM211","AKM104_AK148","AKM333_AKM335_AKM337","AKM171_AKM240","AKM219_AKM221")
sample_labels <- c("No-tag auxin", 
                   "top-2::degron auxin",
                   "top-2::degron etoposide, No-auxin",
                   "top-1::degron auxin",
                   "top-1::degron; top-2::degron auxin")
auxin_colors <- c("#000000","#0033CC","#AA0000","#009999","#9900CC")

auxin_bins <- read.csv("multibins_auxin_100_inputsubt.tab",header=T,sep="\t",quote="'",check.names=T)
# output of deeptools multibigwig summary bins (--outRawCounts, -bs 100)
# applyUnityX divides each data point by the sum of values on X
# ChIPdecay_format params: applyUnityX output, bedfile, windowsize, stepsize
# CI1/2 params: same as ChIPdecay, but additional param 1 or 2 used to compute low or high bounds of 95% confidence interval

auxin_bins_f <- ChIPdecay_format(auxin_bins %>% applyUnityX,rex=rex,windowSize=100000,10000)
auxin_bins_f_CI1 <- ChIPdecay_format_CI(auxin_bins %>% applyUnityX,rex=rex,windowSize=100000,10000,1)
auxin_bins_f_CI2 <- ChIPdecay_format_CI(auxin_bins %>% applyUnityX,rex=rex,windowSize=100000,10000,2)

auxin_bins_f.m <- melt(auxin_bins_f[,-2],id.vars="V2")
auxin_bins_f_CI1.m <- melt(auxin_bins_f_CI1[,-2],id.vars="V2")
auxin_bins_f_CI2.m <- melt(auxin_bins_f_CI2[,-2],id.vars="V2")

auxin_combined_f.m <- auxin_bins_f.m %>%
  cbind(auxin_bins_f_CI1.m$value) %>% 
  cbind(auxin_bins_f_CI2.m$value)
auxin_combined_f.m %>% head
colnames(auxin_combined_f.m) <- c('V2','variable','value','low','high')
auxin_combined_f.m$variable <- factor(auxin_combined_f.m$variable, levels = sample_orders)

gg_auxin_avg <- ggplot(auxin_combined_f.m %>% subset(V2<=500000),
       aes(V2,value,color=variable))+
  geom_line(aes(color=variable),size=1.2)+
  geom_ribbon(aes(ymin=low, 
                  ymax=high,fill=factor(variable)),
              linetype=0,alpha=0.2)+
  scale_x_continuous(labels=addUnits)+theme_bw()+
  scale_fill_manual(values=auxin_colors,names("samples"),labels=sample_labels)+
  scale_color_manual(values=auxin_colors,names("samples"),labels=sample_labels)+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
                     limits=c(3*10^-6, 9.5*10^-6))+
  ylab("DPY-27 ChIP-seq (sum X=1)")+
  xlab("Distance to the nearest strong rex")+
  ggtitle("auxin inputsubt: bin=100bp,window=100kb,step=10kb")
ggsave(plot=gg_auxin_avg,"gg_auxin_avg.pdf",device="pdf",dpi=800,width=10,height=6,unit="in")

####### no auxin avg ####
sample_orders_noAuxin <- c("AKM188_AKM210","AKM102_AKM147","AKM170_AKM239","AKM218_AKM220")
sample_labels_noAuxin <- c("No-tag No-auxin", 
                   "top-2::degron No-auxin",
                   "top-1::degron No-auxin",
                   "top-1::degron; top-2::degron No-auxin")
noAuxin_colors <- c("#000000","#0033CC","#009999","#9900CC")

noAuxin_bins <- read.csv("multibins_noAuxin_100_inputsubt.tab",header=T,sep="\t",quote="'",check.names=T)
noAuxin_bins_f <- ChIPdecay_format(noAuxin_bins %>% applyUnityX,rex=rex,windowSize=100000,10000)
noAuxin_bins_f_CI1 <- ChIPdecay_format_CI(noAuxin_bins %>% applyUnityX,rex=rex,windowSize=100000,10000,1)
noAuxin_bins_f_CI2 <- ChIPdecay_format_CI(noAuxin_bins %>% applyUnityX,rex=rex,windowSize=100000,10000,2)

noAuxin_bins_f.m <- melt(noAuxin_bins_f[,-2],id.vars="V2")
noAuxin_bins_f_CI1.m <- melt(noAuxin_bins_f_CI1[,-2],id.vars="V2")
noAuxin_bins_f_CI2.m <- melt(noAuxin_bins_f_CI2[,-2],id.vars="V2")

noAuxin_combined_f.m <- noAuxin_bins_f.m %>%
  cbind(noAuxin_bins_f_CI1.m$value) %>% 
  cbind(noAuxin_bins_f_CI2.m$value)
noAuxin_combined_f.m %>% head
colnames(noAuxin_combined_f.m) <- c('V2','variable','value','low','high')
noAuxin_combined_f.m$variable <- factor(noAuxin_combined_f.m$variable, levels = sample_orders_noAuxin)

gg_noAuxin_avg <- ggplot(noAuxin_combined_f.m %>% subset(V2<=500000),
       aes(V2,value,color=variable))+
  geom_line(aes(color=variable),size=1.2)+
  geom_ribbon(aes(ymin=low, 
                  ymax=high,fill=factor(variable)),
              linetype=0,alpha=0.2)+
  scale_x_continuous(labels=addUnits)+theme_bw()+
  scale_fill_manual(values=noAuxin_colors,names("samples"),labels=sample_labels_noAuxin)+
  scale_color_manual(values=noAuxin_colors,names("samples"),labels=sample_labels_noAuxin)+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
                     limits=c(3*10^-6, 9.5*10^-6))+
  ylab("DPY-27 ChIP-seq (sum X=1)")+
  xlab("Distance to the nearest strong rex")+
  ggtitle("No-auxin inputsubt: bin=100bp,window=100kb,step=10kb")
ggsave(plot=gg_noAuxin_avg,"gg_noAuxin_avg.pdf",device="pdf",dpi=800,width=10,height=6,unit="in")


##################### etop effect ################
sample_orders_etop <- c("AKM102", "AKM147", "AKM333", "AKM335", "AKM337")
sample_labels_etop <- c("top-2::degron No-auxin, AKM102",
                        "top-2::degron No-auxin, AKM147",
                        "top-2::degron No-auxin,etoposide AKM333",
                        "top-2::degron No-auxin,etoposide AKM335",
                        "top-2::degron No-auxin,etoposide AKM337")
etop_colors <- c("#0033CC","#00A9FF",
                 "#AA0000","#FF68A1","#F8766D")

etop_bins <- read.csv("multibins_etop_100_inputsubt.tab",header=T,sep="\t",quote="'",check.names=T)
head(etop_bins)
etop_bins_f <- ChIPdecay_format(etop_bins %>% applyUnityX,rex=rex,windowSize=100000,10000)
etop_bins_f_CI1 <- ChIPdecay_format_CI(etop_bins %>% applyUnityX,rex=rex,windowSize=100000,10000,1)
etop_bins_f_CI2 <- ChIPdecay_format_CI(etop_bins %>% applyUnityX,rex=rex,windowSize=100000,10000,2)

etop_bins_f.m <- melt(etop_bins_f[,-2],id.vars="V2")
etop_bins_f_CI1.m <- melt(etop_bins_f_CI1[,-2],id.vars="V2")
etop_bins_f_CI2.m <- melt(etop_bins_f_CI2[,-2],id.vars="V2")

etop_combined_f.m <- etop_bins_f.m %>%
  cbind(etop_bins_f_CI1.m$value) %>% 
  cbind(etop_bins_f_CI2.m$value)

colnames(etop_combined_f.m) <- c('V2','variable','value','low','high')
etop_combined_f.m$variable <- factor(etop_combined_f.m$variable, levels = sample_orders_etop)


gg_etoposide_reps<- ggplot(etop_combined_f.m %>% subset(V2<=500000),
       aes(V2,value,color=variable))+
  geom_line(aes(color=variable),size=1.2)+
  geom_ribbon(aes(ymin=low, 
                  ymax=high,fill=factor(variable)),
              linetype=0,alpha=0.2)+
  scale_x_continuous(labels=addUnits)+theme_bw()+
  scale_fill_manual(values=etop_colors,names("samples"),labels=sample_labels_etop)+
  scale_color_manual(values=etop_colors,names("samples"),labels=sample_labels_etop)+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
                     limits=c(3*10^-6, 9.5*10^-6))+
  ylab("DPY-27 ChIP-seq (sum X=1)")+
  xlab("Distance to the nearest strong rex")+
  ggtitle("etoposide inputsubt: bin=100bp,window=100kb,step=10kb")

ggsave(plot=gg_etoposide_reps,"gg_etoposide_reps.pdf",device="pdf",dpi=800,width=10,height=6,unit="in")
