library(ggplot2)
library(gridExtra)

mytheme = theme(
  #	plot.background = element_rect(fill = 'black', color = "black"), 
  axis.ticks.y = element_line(color = "black", linewidth = 1), 
  axis.ticks.x = element_line(color = "black", linewidth = 1), 
  axis.ticks.length = unit(0.5, "cm"), 
  axis.line.x = element_line(color = "black", linewidth = 2), 
  axis.line.y = element_line(color = "black", linewidth = 2), 
  axis.text.x = element_text(color = 'black', size  = 28, margin = margin(rep(4, 4))), 
  axis.text.y = element_text(color = 'black', size = 28, margin = margin(rep(2, 4))), 
  axis.title.x = element_text(color = "black", size = 24, margin = margin(rep(15, 4))), 
  axis.title.y = element_text(color = "black", size = 24,  margin = margin(rep(15, 4))), 
  plot.title = element_text(color = "black", size = 24, margin = margin(rep(15, 4))), 
  panel.background = element_rect(fill = "white"), 
  #	panel.border = element_blank(),
  panel.grid.major.y = element_blank(), 
  panel.grid.minor.y = element_blank(), 
  panel.grid.major.x = element_blank(), 
  panel.grid.minor.x = element_blank(),
  legend.background = element_blank())

tp<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/cor_Comeron_A6_A7_rate_depth_normalized.txt",header=F)#used this version, exclude multiple reads map to the same events doesn't change much. column4 A7, column5 A6
rec_comeron<-data.frame(tp$V1,tp$V2,tp$V3/10)
rec_comeron <- setNames(rec_comeron, c("chr","posi","rec"))
rec_comeron$strain <- "comeron"
chro<-rec_comeron[grep("3R", rec_comeron$chr),]
p1 = ggplot(data = chro, aes(x = posi, y = rec,col=strain))+geom_smooth(method = "loess",span = 0.2, se=F)+ xlab("3R")+ ylab("cM/Mb")+geom_point(size=1)+scale_color_manual(values = c("#818181"))
p1<-p1+mytheme
p1
ggsave(file = "/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/3R_loess_0.2_comeron_cM_Mb.pdf", plot = p1, width = 14, height = 4)

rec_A6<-data.frame(tp$V1,tp$V2,tp$V4*100)
rec_A6 <- setNames(rec_A6, c("chr","posi","rec"))
rec_A6$strain <- "A6"
rec_A7<-data.frame(tp$V1,tp$V2,tp$V5*100)
rec_A7 <- setNames(rec_A7, c("chr","posi","rec"))
rec_A7$strain <- "A7"
chro_sum<-rbind(rec_A6,rec_A7)
chro<-chro_sum[grep("X", chro_sum$chr),]
p1 = ggplot(data = chro, aes(x = posi, y = rec,col=strain))+geom_smooth(method = "loess",span = 0.2, se=F)+ xlab("X")+ ylab("normalized CO")+geom_point(size=1)+scale_color_manual(values = c("#CC79A7", "#56B4E9"))
p1<-p1+mytheme
p1
ggsave(file = "/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/X_loess_0.2_A6_A7_normalized_CO.pdf", plot = p1, width = 14, height = 4)


rec_pool<-data.frame(tp$V1,tp$V2,(tp$V4*100+tp$V5*100)/2)
rec_pool <- setNames(rec_pool, c("chr","posi","rec"))
rec_pool$strain <- "pool"
chro<-rec_pool[grep("3R", rec_pool$chr),]
p1 = ggplot(data = chro, aes(x = posi, y = rec,col=strain))+geom_smooth(method = "loess",span = 0.2, se=F)+ xlab("3R")+ ylab("normalized CO")+geom_point(size=1)+scale_color_manual(values = c("#b0d146"))
p1<-p1+mytheme
p1
ggsave(file = "/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/3R_loess_0.2_average_pool_normalized_CO.pdf", plot = p1, width = 14, height = 4)


chro_sum<-rbind(rec_comeron,rec_pool)
cor.test(tp$V3,(tp$V4+tp$V5)/2,method = "spearman")
cor.test(tp$V3,tp$V5,method = "spearman")

chro<-chro_sum[grep("3R", chro_sum$chr),]
p1 = ggplot(data = chro, aes(x = posi, y = rec,col=strain))+geom_smooth(method = "loess",span = 0.2, se=F)+ xlab("3R")+ ylab("cM/Mb")+geom_point(size=1)+scale_color_manual(values = c("#bfbfbf","#b0d146"))
#p1 = ggplot(data = chro, aes(x = posi, y = rec,col=strain))+geom_smooth(method = "loess",span = 0.2, se=F)+ xlab("X")+ ylab("cM/Mb")+geom_point(size=1)+scale_color_manual(values = c("#CC79A7", "#56B4E9","#bfbfbf"))
p1<-p1+mytheme
p1
ggsave(file = "/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/3R_loess_0.2_comeron_avg_pools_cM_Mb.pdf", plot = p1, width = 14, height = 6)


