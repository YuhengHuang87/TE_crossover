#### Within strains comparisons
####CO numbers 
focal_depth<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7_crossover_num_TE_included_window_exclude_TE_overlapped_average_5000.txt",header=F)
con_depth<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7_crossover_num_control_included_window_exclude_TE_overlapped_average_10000.txt",header=F)
#con_depth<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A6_crossover_num_control_included_window_exclude_TE_overlapped_average_10000_30000_40000",header=F)

focal_depth<-focal_depth[focal_depth$V7 >= 5 & focal_depth$V9 >= 5, ];con_depth<-con_depth[con_depth$V5 >= 5 & con_depth$V7 >= 5, ];
focal_depth$leng=focal_depth$V4-focal_depth$V3+1;
focal_depth<-focal_depth[focal_depth$leng>=200, ];
foo <- data.frame(do.call('rbind', strsplit(as.character(focal_depth$V5),'/',fixed=TRUE)));focal_depth$class <- foo$X1
focal_depth<-focal_depth[focal_depth$class != "Unknown", ];focal_depth$TE_type<-ifelse(focal_depth$class == "DNA","DNA","RNA")

#focal_depth<-focal_depth[focal_depth$V1 != "nearby_ambiguous", ]
#focal_rec<-focal_depth[focal_depth$V2 == "CM010573.1", ];nrow(focal_rec)
focal_depth$co=focal_depth$V10+focal_depth$V12;con_depth$co=con_depth$V8+con_depth$V10

#ks.test(focal_depth$co,con_depth$co,alternative= "greater")

focal<-focal_depth;control<-con_depth
rec_RNA<-focal[focal$TE_type == "RNA",];rec_DNA<-focal[focal$TE_type == "DNA",];
rec_LTR<-focal[focal$class == "LTR",];rec_LINE<-focal[focal$class == "LINE",]
focal$group<-"foc";control$group<-"con"; 
focal$snp=focal$V7+focal$V9;control$snp=control$V5+control$V7;
focal$depth=(focal$V6*focal$V7+focal$V8*focal$V9)/(focal$V7+focal$V9);control$depth=(control$V4*control$V5+control$V6*control$V7)/(control$V5+control$V7);
#focal$depth=(focal$V6+focal$V8)/2;control$depth=(control$V4+control$V6)/2;

q_con<-quantile(control$depth,probs = c(0.05,0.95));q_foc<-quantile(focal$depth,probs = c(0.05,1));
control<-control[control$depth>q_con[1] & control$depth<q_con[2],];focal<-focal[focal$depth>q_foc[1] & focal$depth<q_foc[2],];

foc<-data.frame(focal$depth,focal$co,focal$group,focal$snp)
con<-data.frame(control$depth,control$co,control$group,control$snp)
total<-rbind(foc, setNames(con, names(foc)))
wilcox.test(focal$depth,control$depth,  exact=FALSE, correct=FALSE, paired = FALSE);median(focal$depth);median(control$depth)
wilcox.test(focal$snp,control$snp,  exact=FALSE, correct=FALSE, paired = FALSE);median(focal$snp);median(control$snp)
wilcox.test(focal$co,control$co,alternative= "less",  exact=FALSE, correct=FALSE, paired = FALSE);mean(focal$co);mean(control$co)

### Fig. S6
#CC79A7
#56B4E9
total$focal.group <- factor(total$focal.group, levels = c('foc','con'),ordered = TRUE)##CC79A7,#56B4E9
dp <- ggplot(total, aes(x=focal.group, y=focal.snp, fill=focal.group,color=focal.group)) + scale_fill_manual(values=c('#CC79A7','grey'))+
  #geom_boxplot(color = "black")+
  geom_violin(trim=FALSE,color = "black")+
  labs(x="TE or not", y = "SNP")
dp<-dp+mytheme;dp
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/SNP_to_CO_5kb_control_A6.pdf", sep = ""), plot = dp, width = 7, height = 5)
###

library(MASS);library("dplyr") 
m2 <- glm.nb(focal.co ~ focal.group+focal.depth+focal.snp,  data=total);summary(m2)
m2 <- glm.nb(focal.co ~ focal.group+focal.snp,  data=total);summary(m2)

m_class<-glm.nb(co ~ TE_type, data=focal);summary(m_class); #anova(m_class); 
m_class<-glm.nb(co ~ class, data=focal);summary(m_class); #anova(m_class); 
##check the dispersion of the model
E2 <- resid(m2, type = "pearson")
N  <- nrow(total)
p  <- length(coef(m2))+1   # '+1' is for variance parameter in NB
sum(E2^2) / (N - p)

100*sum(focal$co == 0)/nrow(focal);100*sum(control$co == 0)/nrow(control)

#### Fig. 3B
total$focal.group <- factor(total$focal.group, levels = c('focal','con'),ordered = TRUE)##CC79A7,#56B4E9
library(ggplot2)
library(gridExtra)

mytheme = theme(
  #	plot.background = element_rect(fill = 'black', color = "black"), 
  axis.ticks.y = element_line(color = "black", linewidth = 0.45), 
  axis.ticks.x = element_line(color = "black", linewidth = 0.45), 
  axis.ticks.length = unit(0.2, "cm"), 
  axis.line.x = element_line(color = "black", linewidth = 0.8), 
  axis.line.y = element_line(color = "black", linewidth = 0.8), 
  axis.text.x = element_text(color = 'black', size  = 12, margin = margin(rep(4, 4))), 
  axis.text.y = element_text(color = 'black', size =  12, margin = margin(rep(2, 4))), 
  axis.title.x = element_text(color = "black", size = 7, margin = margin(rep(15, 4))), 
  axis.title.y = element_text(color = "black", size = 7,  margin = margin(rep(15, 4))), 
  plot.title = element_text(color = "black", size =   7, margin = margin(rep(15, 4))), 
  legend.position = c(1.15, 0.9),
  aspect.ratio = 0.8,
  panel.background = element_rect(fill = "white"), 
  #	panel.border = element_blank(),
  panel.grid.major.y = element_blank(), 
  panel.grid.minor.y = element_blank(), 
  panel.grid.major.x = element_blank(), 
  panel.grid.minor.x = element_blank(),
  legend.background = element_blank())

dp <- ggplot(total, aes(x=focal.group, y=focal.co,color=focal.group)) +
  #geom_count()+
  geom_count(aes(size = after_stat(prop)))+ scale_size_area()+ 
  scale_color_manual(values=c('#56B4E9','grey'))+
  ylim(-0.5,4)+
  mytheme
dp
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/CO_within_strain_A7.pdf", sep = ""), plot = p2, width = 7, height = 5)


####
rec<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7_crossover_num_TE_left_right_window_magnitude_5000_5000_include_ambiguous.txt",header=F)#CO & magnitude
foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V13),'/',fixed=TRUE)));rec$class <- foo$X1;rec<-rec[rec$class != "Unknown", ];rec$TE_type<-ifelse(rec$class == "DNA","DNA","RNA")
rec$K9_mag_left<-(rec$V17+rec$V18+rec$V19)/3
#rec$K9_mag_right<-(rec$V20+rec$V21+rec$V22)/3
rec$K9_mag_right<-(rec$V21+rec$V22+rec$V23)/3
rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2
rec$co<-rec$V5+rec$V7
rec$snp<-rec$V2+rec$V4;rec$cov = (rec$V2*rec$V1+rec$V4*rec$V3)/(rec$snp)
rec<-rec[rec$V2 >= 5 & rec$V4 >= 5, ]

rec$leng = rec$V12-rec$V11+1
rec_as<-rec[rec$leng >= 200, ]
rec_RNA<-rec_as[rec_as$TE_type == "RNA",];rec_DNA<-rec_as[rec_as$TE_type == "DNA",];
rec_LTR<-rec_as[rec_as$class == "LTR",];rec_LINE<-rec_as[rec_as$class == "LINE",]
mean(rec_DNA$K9);mean(rec_RNA$K9);mean(rec_LINE$K9);mean(rec_LTR$K9)
wilcox.test(rec_LINE$K9,rec_DNA$K9, exact=FALSE, correct=FALSE, paired = F)
wilcox.test(rec_DNA$K9,rec_LTR$K9, exact=FALSE, correct=FALSE, paired = F)
wilcox.test(rec_LINE$K9,rec_LTR$K9, exact=FALSE, correct=FALSE, paired = F)

rec$group<-"foc";rec$depth<-rec$cov
rec_as<-rec[rec$leng >= 200, ]
rec_as<-rec_as[rec_as$K9<1,]

foc<-data.frame(rec_as$depth,rec_as$co,rec_as$group,rec_as$snp)
con<-data.frame(control$depth,control$co,control$group,control$snp)
total<-rbind(foc, setNames(con, names(foc)))
m2 <- glm.nb(rec_as.co ~ rec_as.group+rec_as.depth+rec_as.snp,  data=total);summary(m2)

m2 <- glm.nb(co ~ K9+class+K9:class,  data=rec_as);summary(m2)
m2 <- glm.nb(co ~ K9,  data=rec_as);summary(m2)
m2 <- glm.nb(co ~ K9,  data=rec_DNA);summary(m2)
m2 <- glm.nb(co ~ K9,  data=rec_RNA);summary(m2)
m2 <- glm.nb(co ~ K9,  data=rec_LINE);summary(m2)
m2 <- glm.nb(co ~ K9,  data=rec_LTR);summary(m2)
m2 <- lm(K9 ~ class,  data=rec_as);summary(m2)
m2 <- lm(K9 ~ TE_type,  data=rec_as);summary(m2)

K9_cutoff=1;rec_as$K9_type<-ifelse(rec_as$K9>K9_cutoff,"sign_en","no_en")
rec_enrich<-rec_as[rec_as$K9_type == "sign_en",];rec_un<-rec_as[rec_as$K9_type == "no_en",];wilcox.test(rec_enrich$co,rec_un$co, alternative ="less",exact=FALSE, correct=FALSE, paired = F);mean(rec_enrich$co);mean(rec_un$co)
m2 <- glm.nb(co ~ K9_type,  data=rec_as);summary(m2)
m2 <- glm.nb(co ~ K9+snp+cov+class,  data=rec_as);summary(m2)
m2 <- glm.nb(co ~ leng,  data=rec_un);summary(m2)

cor.test(rec_as$K9,rec_as$co,method="spearman")
cor.test(rec_RNA$K9,rec_RNA$co,method="spearman")
cor.test(rec_DNA$K9,rec_DNA$co,method="spearman")
cor.test(rec_LTR$K9,rec_LTR$co,method="spearman")
cor.test(rec_LINE$K9,rec_LINE$co,method="spearman")
cor.test(rec_un$leng,rec_un$co,method="spearman")


DNA_TE<-rec_as[rec_as$class == "DNA", ];LINE_TE<-rec_as[rec_as$class == "LINE", ];LTR_TE<-rec_as[rec_as$class == "LTR", ];
DNA<-rec_as[rec_as$TE_type == "DNA", ];RNA<-rec_as[rec_as$TE_type == "RNA", ];


#### comparing TEs with and without K9 enrichment
foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V13),'/',fixed=TRUE)));rec$class <- foo$X1
#rec_as<-rec_as[rec_as$class != "Unknown" & rec_as$class != "DNA", ];
rec_as<-rec_as[rec_as$class != "Unknown", ];rec$TE_type<-ifelse(rec$class == "DNA","DNA","RNA")

rec<-rec[rec$V2 >= 5 & rec$V4 >= 5, ]
rec$K9_mag_left<-(rec$V17+rec$V18+rec$V19)/3
#rec$K9_mag_right<-(rec$V20+rec$V21+rec$V22)/3
rec$K9_mag_right<-(rec$V21+rec$V22+rec$V23)/3

rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2
rec$cov<-(rec$V3+rec$V1)/2
rec$co<-rec$V5+rec$V7
rec$snp<-rec$V2+rec$V4
rec$leng = rec$V12-rec$V11+1
rec_as<-rec[rec$leng >= 200, ]
rec_as<-rec_as[rec_as$V9 != "roo", ]
rec_as$K9_type<-ifelse(rec_as$K9>1,"sign_en","no_en")
rec_enrich<-rec_as[rec_as$K9_type == "sign_en",];rec_un<-rec_as[rec_as$K9_type == "no_en",];wilcox.test(rec_enrich$co,rec_un$co, exact=FALSE, correct=FALSE, paired = F)
cor.test(rec_as$K9,rec_as$snp,method="spearman");cor.test(rec_as$K9,rec_as$cov,method="spearman");#K9_m<-lm(K9 ~ cov+snp, data=rec_as);summary(K9_m); 

K9_m<-lm(K9 ~ TE_type+ depth + snp_den, data=rec_as);summary(K9_m); 
K9_m<-lm(log2(K9) ~ class+ depth + snp_den, data=rec_as);summary(K9_m); 
K9_m<-glm(K9 ~ class+ depth + snp_den, data=rec_as,family = Gamma("inverse"));summary(K9_m); 

a1 <- aov(K9 ~ class , data=rec_as);summary(a1)
posthoc <- TukeyHSD(x=a1, 'class', conf.level=0.95);posthoc

#m1 <- glm(co ~ class, family = poisson(link = "log"), data=rec_as);summary(m1)
#m1 <- glm(co ~ K9, family = poisson(link = "log"), data=rec_as);summary(m1)
m2 <- glm.nb(co ~ K9+cov+snp,  data=rec_as);summary(m2)
m2 <- glm.nb(co ~ class,  data=rec_as);summary(m2)

rec_nonK9_TE<-rec_as[rec_as$K9<1,]
m2 <- glm.nb(co ~ leng,  data=rec_as);summary(m2)

##### Fig. 3D
#CC79A7 #56B4E9 +xlim(0,20)
p<-ggplot(aes(x = K9, y = co), data = rec)+
  xlim(0,16) + 
  geom_point(size=3,alpha = 0.5,color="#CC79A7")+ 
  #theme(axis.line.x = element_line(color="black", size = 0.7), axis.line.y = element_line(color="black", size = 0.7),axis.title.x = element_text(face="bold",  size=20),axis.title.y = element_text(angle=90,size=20,face="bold"), axis.text.x  = element_text(vjust=1, size=20,colour="Black"),axis.text.y = element_text(hjust=1, size=20,colour="Black"), legend.text = element_text(hjust=20,size=20,face="italic"),legend.title = element_text(hjust=0,size=20,face="bold"),panel.background = element_rect(fill = "white"))+
  ylab("CO number")+xlab("TE epigenetic magnitude")
p<-p+mytheme; p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/A6_within_genome_K9_TE_200bp_remove_nearby_ambi.pdf", sep = ""), plot = p, width = 6.2, height = 5)

##############
### Distance to the nearest CO
rec<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7_nearest_second_distance_TE_included_window_exclude_TE_overlapped_average_5000.txt",header=F)
foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V5),'/',fixed=TRUE)));rec$class <- foo$X1;rec<-rec[rec$class != "Unknown", ];rec$TE_type<-ifelse(rec$class == "DNA","DNA","RNA")

rec$K9_mag_left<-abs(rec$V10+rec$V11+rec$V12)/3
rec$K9_mag_right<-(rec$V13+rec$V14+rec$V15)/3
rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2

rec$log2_dis=log2(rec$V6)
rec$snp<-rec$V9;rec$snp_den=rec$snp/rec$V6;rec$cov = rec$V8
rec$leng = rec$V4-rec$V3+1
rec_as<-rec[rec$leng >= 200, ]
rec_as<-rec_as[rec_as$K9 < 3,]
cor.test(rec$cov,rec$log2_dis,method="spearman")


con<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/control_A7_nearest_second_distanc_exclude_TE_overlapped_average_10000.txt",header=F)
con$log2_dis=log2(con$V4);
wilcox.test(rec_as$log2_dis,con$log2_dis, alternative ="greater",exact=FALSE, correct=FALSE, paired = F); median(rec_as$log2_dis);median(con$log2_dis)
rec_as$win<-"focal";con$win<-"without_TE";
foc<-data.frame(rec_as$log2_dis,rec_as$snp_den,rec_as$cov,rec_as$win)
con<-data.frame(con$log2_dis,con$V7,con$V6,con$win)
total<-rbind(foc, setNames(con, names(foc))); 
m2 <- lm(rec_as.log2_dis ~ rec_as.win+rec_as.cov+rec_as.cov+rec_as.snp_den,  data=total);summary(m2)
mean(rec_as$log2_dis);mean(con$con.log2_dis)
cor.test(con$V6,con$log2_dis,method= "spearman")

####Fig. 3C
#56B4E9 #CC79A7
dp <- ggplot(total, aes(x=rec_as.win, y=rec_as.log2_dis, fill=rec_as.win,color=rec_as.win)) + scale_fill_manual(values=c('#CC79A7','grey'))+
  #geom_boxplot(color = "black")+
  geom_violin(trim = F, color = "black")+
  
  labs(x="TE or not", y = "log2 distance to CO")
dp<-dp+mytheme;dp
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/log2_distance_to_CO_A6_violin.pdf", sep = ""), plot = dp, width = 7, height = 5)
#####
#compare different TE classes
rec_RNA<-rec_as[rec_as$TE_type == "RNA",];rec_DNA<-rec_as[rec_as$TE_type == "DNA",];
rec_LTR<-rec_as[rec_as$class == "LTR",];rec_LINE<-rec_as[rec_as$class == "LINE",]
mean(rec_as$log2_dis);
mean(rec_DNA$K9);mean(rec_RNA$K9);mean(rec_LINE$K9);mean(rec_LTR$K9);
wilcox.test(rec_LINE$K9,rec_DNA$K9, exact=FALSE, correct=FALSE, paired = F)
wilcox.test(rec_DNA$K9,rec_LTR$K9, exact=FALSE, correct=FALSE, paired = F)
wilcox.test(rec_LINE$K9,rec_LTR$K9, exact=FALSE, correct=FALSE, paired = F)

wilcox.test(rec_RNA$log2_dis,rec_DNA$log2_dis, alternative ="greater",exact=FALSE, correct=FALSE, paired = F); mean(rec_RNA$log2_dis);mean(rec_DNA$log2_dis)
mean(rec_DNA$log2_dis);mean(rec_RNA$log2_dis);mean(rec_LTR$log2_dis);mean(rec_LINE$log2_dis)

m2 <- lm(log2_dis ~ K9+cov+snp_den,  data=rec_DNA);summary(m2)
m2 <- lm(log2_dis ~ K9+class+K9:class,  data=rec_as);summary(m2)
m2 <- lm(log2_dis ~ K9+TE_type+K9:TE_type,  data=rec_as);summary(m2)
m2 <- lm(log2_dis ~ K9,  data=rec_as);summary(m2)

m2 <- lm(log2_dis ~ TE_type+cov+snp_den,  data=rec_as);summary(m2)
m2 <- lm(log2_dis ~ class+cov+snp_den,  data=rec_as);summary(m2)
m2 <- lmer(log2_dis ~ K9+cov+snp_den + (1|V1),  data=rec_as);summary(m2)

m2 <- lm(log2_dis ~ K9+cov+snp_den,  data=rec_as);summary(m2)
m2 <- lm(log2_dis ~ K9+cov+snp_den,  data=rec_DNA);summary(m2)
m2 <- lm(log2_dis ~ K9+cov+snp_den,  data=rec_RNA);summary(m2)
m2 <- lm(log2_dis ~ K9+cov+snp_den,  data=rec_LINE);summary(m2)
m2 <- lm(log2_dis ~ K9+cov+snp_den,  data=rec_LTR);summary(m2)

cor.test(rec_as$K9,rec_as$log2_dis,method="spearman")
cor.test(rec_DNA$log2_dis,rec_DNA$K9,method="spearman")
cor.test(rec_RNA$log2_dis,rec_RNA$K9,method="spearman")
cor.test(rec_LINE$log2_dis,rec_LINE$K9,method="spearman")
cor.test(rec_LTR$log2_dis,rec_LTR$K9,method="spearman")
m2 <- lm(K9 ~ TE_type,  data=rec_as);summary(m2);m2 <- lm(K9 ~ class,  data=rec_as);summary(m2)


DNA_TE<-rec_as[rec_as$class == "DNA", ];LINE_TE<-rec_as[rec_as$class == "LINE", ];LTR_TE<-rec_as[rec_as$class == "LTR", ];
DNA<-rec_as[rec_as$TE_type == "DNA", ];RNA<-rec_as[rec_as$TE_type == "RNA", ];
wilcox.test(RNA$dis,DNA$dis,exact=FALSE, correct=FALSE, paired = F);median(RNA$dis);median(DNA$dis)
wilcox.test(RNA$depth,DNA$depth,exact=FALSE, correct=FALSE, paired = F);median(RNA$depth);median(DNA$depth)
wilcox.test(LINE_TE$dis,LTR_TE$dis, exact=FALSE, correct=FALSE, paired = F);median(LINE_TE$dis);median(LTR_TE$dis);median(DNA_TE$dis)
#m2 <- lm(log2(dis) ~ TE_type,  data=rec_as);summary(m2) #not sig
m2 <- lm(log2(dis) ~ class+depth+snp_den,  data=rec_as);summary(m2)
m2 <- lm(log2(dis) ~ TE_type+depth+snp_den,  data=rec_as);summary(m2)
m0 <- lm(log2(dis) ~ depth+snp_den,  data=rec_as);summary(m0); anova(m2, m0)
posthoc <- TukeyHSD(x=a1, 'class', conf.level=0.95);posthoc
m0 <- lm(snp_den ~ class, data = rec_as); summary(m0)
sum(rec_as$class == "DNA")/nrow(rec_as)


#compare TE with high and low K9 mass
rec_as$K9_type<-ifelse(rec_as$K9>median(rec_as$K9,na.rm = T),"sign_en","no_en")
rec_enrich<-rec_as[rec_as$K9_type == "sign_en",];rec_un<-rec_as[rec_as$K9_type == "no_en",];wilcox.test(rec_enrich$log2_dis,rec_un$log2_dis, alternative ="greater",exact=FALSE, correct=FALSE, paired = F);mean(rec_enrich$log2_dis);mean(rec_un$log2_dis)
m2 <- lm(log2_dis ~ K9+leng+cov+snp_den,  data=rec_as);summary(m2)
m2 <- lm(log2_dis ~ K9,  data=rec_as);summary(m2)
m2 <- lm(log2_dis ~ leng,  data=rec_un);summary(m2)
cor.test(rec_un$leng,rec_un$log2_dis,method="spearman")

#### Fig. 3E
cols <- c("grey","#CC79A7")
cols <- c("grey","#56B4E9")
p<-ggplot(aes(x = K9, y = log2_dis), data = rec_as)+
  #xlim(0,16) + 
  geom_point(size=3,alpha = 0.5,color="#56B4E9")+stat_smooth(method = "lm",    formula = y ~ x,  geom = "smooth",se = FALSE,color="#56B4E9") + 
  ylab("log2 distance to CO")+xlab("TE epigenetic magnitude")
p<-p+mytheme; p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/A7_within_genome_K9_TE_200bp_log2_distance.pdf", sep = ""), plot = p, width = 6.2, height = 5)

m2 <- lm(log2_dis ~ log2(K9)+cov+snp_den,  data=rec_as);summary(m2)
m2 <- lm(log2(dis) ~ K9+depth+snp_den,  data=rec_as);summary(m2)
m3 <- lm(log2(dis) ~ leng+depth+snp_den,  data=rec_un);summary(m3)
cor.test(log2(rec_un$dis),rec_un$leng,method = "spearman")
#A7_spreading<-rec_as$K9;A6_mag<-rec_as$K9;cor.test(A6_mag,rec_as$K9,method="spearman")
#######

cor.test(rec_as$K9,rec_as$snp,method="spearman");cor.test(rec_as$K9,rec_as$cov,method="spearman");#K9_m<-lm(K9 ~ cov+snp, data=rec_as);summary(K9_m); 
K9_m<-lm(K9 ~ TE_type+ depth + snp_den, data=rec_as);summary(K9_m); 
K9_m<-lm(log2(K9) ~ class+ depth + snp_den, data=rec_as);summary(K9_m); 
K9_m<-glm(K9 ~ class+ depth + snp_den, data=rec_as,family = Gamma("inverse"));summary(K9_m); 

a1 <- aov(K9 ~ class , data=rec_as);summary(a1)
posthoc <- TukeyHSD(x=a1, 'class', conf.level=0.95);posthoc

