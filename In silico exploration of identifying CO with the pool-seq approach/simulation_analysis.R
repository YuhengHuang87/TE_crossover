
#simulation

############# calulating the number of reads for downsampling
depth<-c(0.1,0.25,0.5,1,2,3,4,6,8)
depth=3
scale_factor = 10.6/6.6
#without screening, CO reads: mean depth = 10.605146572872,  median depth = 9
#after screening, CO reads: mean depth = 6.62852202450794, median depth = 4
screen_depth=15/scale_factor
#n_cutoff #0bp: 16041;1.5kb: 13261; 2kb:  12385; 2.5kb: 11491
n_cutoff<-c(16041,13261,12385,11491)
read=depth/screen_depth*n_cutoff
round(read)


###### Fig. 2A & C
CO_recover<-read.csv("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/CO_event_recovered_simulation.csv",header=T)
#scaleFUN <- function(x) sprintf("%.3f", x)
demo_log10(
  c(0.0, 0.2, 0.4, 0.6, 0.8,1.0),
  labels = label_number(drop0trailing = TRUE)
)
p<-ggplot(CO_recover, aes(x = depth, y = prop_recover)) + ylab("proportion of CO recalled")+ xlab("average depth per haplotype")+geom_point(size = 4)+ geom_smooth(method = "loess",span = 0.5,color = "grey",size =1.5,se=F)
p1<-p+mytheme+scale_y_continuous(limits = c(0,1),labels = label_number(accuracy = 0.1))+geom_vline(xintercept=3, linetype="solid", color = "#e69f00", size=1.2)+geom_vline(xintercept=0.09, linetype="solid", color = "#CC79A7", size=1.2)+geom_vline(xintercept=0.15, linetype="solid", color = "#56B4E9", linewidth=1.2)
p1
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/proportion_CO_recalled_solid_line.pdf", sep = ""), plot = p1, width = 7.5, height = 6)

p<-ggplot(CO_recover, aes(x = depth, y = total_depth_to_detect_a_event)) + ylab("relative cost per recalled CO")+ xlab("average depth per haplotype")+geom_point(size = 4)+ geom_smooth(method = "loess",span = 0.5,color = "grey",size =1.5,se=F)
p<-p+geom_vline(xintercept=3, linetype="solid", color = "#e69f00", size=1.2)+geom_vline(xintercept=0.09, linetype="solid", color = "#CC79A7", size=1.2)+geom_vline(xintercept=0.15, linetype="solid", color = "#56B4E9", size=1.2)
p1<-p+mytheme
p1
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/cost_per_CO_recalled_solid_line.pdf", sep = ""), plot = p1, width = 7.5, height = 6)
#######


### Fig. 2B
CO_reads<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/range_depth_event_identified_read_num.txt",header=F)
p<-ggplot(CO_reads, aes(x = V1, y = V7)) + ylab("number of reads per CO")+ xlab("average depth per haplotype")
#+geom_jitter(size = 1.5,width = 0.1, alpha = 0.1)
p<-p+geom_vline(xintercept=3, linetype="solid", color = "#e69f00", size=1.2)+geom_vline(xintercept=0.09, linetype="solid", color = "#CC79A7", size=1.2)+geom_vline(xintercept=0.15, linetype="solid", color = "#56B4E9", size=1.2)
p1<-p+mytheme+geom_boxplot(aes(group = V1),width = 0.1, alpha = 0.5, outlier.shape = NA)+ geom_smooth(method = "loess",span = 0.63,color = "grey",size =1.5,se=F)+ylim(1,11)
p1


library(dplyr)
data_subset <- CO_reads |>  filter(V1 %in% c(0.1, 0.25, 0.5))
data_rest <- CO_reads |>  filter(!V1 %in% c(0.1, 0.25, 0.5))
CO_reads_summary <- CO_reads |>
  group_by(V1, V7) |>
  summarise(count = n(), .groups = 'drop')
p1<-ggplot() +
  geom_boxplot(data = data_rest, aes(x = V1, y = V7, group = V1),
               width = 0.4, alpha = 0.3, outlier.shape = NA,
               lwd = 1, fatten = 1) +
  geom_point(data = CO_reads_summary, aes(x = V1, y = V7, size = count),
             color = "black", alpha = 0.5, stroke = 0) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8)) +
  labs(x = "average depth per haplotype",
       y = "number of reads per CO") +
  geom_vline(xintercept=3, linetype= "solid", color = "#E69F00", linewidth= 1.2) +
  geom_vline(xintercept=0.09, linetype= "solid", color = "#CC79A7", linewidth= 1.2) +
  geom_vline(xintercept=0.15, linetype= "solid", color = "#56B4E9", linewidth= 1.2) +
  mytheme 

  #theme(axis.title.x = element_blank(),axis.title.y = element_blank())
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/reads_per_CO_recalled_boxplot_KLversion.pdf", sep = ""), plot = p1, width = 8.3, height = 6)
##############


###########

n=538;
mat<-rmultinom(n=1000, size=736, prob=rep(1/n,n))

mat<-rmultinom(n=1000, size=578*0.69, prob=rep(1/192,192))
mat<-rmultinom(n=100, size=1410, prob=rep(1/17852,17852))
mat<-rmultinom(n=1000, size=600, prob=runif(192))

write.table(mat, file="/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/equal_pool_CO_benchmark_Illumina_resample.txt",sep="\t",eol="\n",row.names = F,col.names = F)

m=1;trial<-mat[,m];
length(trial[trial==0]);length(trial[trial==1]);length(trial[trial==2]);length(trial[trial==3]);length(trial[trial==4]);
posi<-trial[trial>0]
hist(posi)
posi=as.data.frame(posi)
p<-ggplot(posi, aes(x=posi)) + 
  geom_histogram(color="black", fill="white",binwidth = 1)
p
sum(mat[,2])

mylist=c()
for (i in 1:ncol(mat)){
  mySample=0
  for (j in 1:nrow(mat)){
    if (mat[j,i]>0){
      mySample=mySample+mat[j,i]
    }
  }
  mylist=c(mylist,mySample)
}
mylist;mean(mylist)

#####Fig. S4
bin_co_sim<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/15_co_region_false_negative_rate_nonoverlap_30.txt",header=F)#
bin_co_sim<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/3_co_region_false_negative_rate_nonoverlap_21.txt",header=F)#
#chr<-bin_co_sim[grep("CM010569.1", bin_co_sim$V1),]
chr<-bin_co_sim
#co_sum=sum(chr$V2);event_sum=sum(chr$V3)
myP=c()
for (i in 1:nrow(chr)){
  #x <- matrix(c(chr[i,2],       chr[i,3],       co_sum-chr[i,2],       event_sum-chr[i,3]), ncol = 2) 
  x <- matrix(c(chr[i,6],       chr[i,7],       chr[i,8],       chr[i,9]), ncol = 2) 
  p_value=chisq.test(x)$p.value
  myP=c(myP,p_value)
}
p_adj<-p.adjust (myP, method="fdr") 
sum(p_adj<0.05)
which(p_adj<0.05)
sig<-p_adj[p_adj<0.05];sig

share<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/3_event_identified_2000",header=F)#
unique<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/3_event_unidentified_2000",header=F)#

share$V5<-"share";
unique$V5<-"unique";
share$V6<-"1";
unique$V6<-"0";

total=rbind(share,unique)
data<-total[grep("CM010569.1", total$V2),]
nrow(data)
p<-ggplot(data, aes(x=V3, y = 0, color=V5,size=V5))+scale_color_manual(values=c("#D55E00", "#009E73")) + geom_point(shape="|")+scale_size_manual(values=c(7,6))#+scale_y_continuous(expand = c(0, 0))
p1<-p+theme_minimal()+geom_hline(aes(yintercept = 1))+ theme(axis.title.x = element_text(face="bold",  size=18),axis.title.y = element_text(angle=90,size=18,face="bold"), axis.text.x  = element_text(vjust=1, size=16,colour="Black"), legend.text = element_text(hjust=20,size=20,face="italic"),legend.title = element_text(hjust=0,size=20,face="bold"))+
  xlab("X")
p1
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/PBSIM_X_co_gc_nMaxL.pdf", sep = ""), plot = p1, width = 20, height = 3)

snp_den<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/CM010570.1_snp_density_aln_A6_refer_A4",header=F)
#snp_den<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/CM010570.1_snp_density_aln_A6_refer_A4_sliding_windows",header=F)
p<-ggplot(snp_den, aes(x = V2, y = V4/100000)) + ylab("SNP density")+ xlab("2L")+geom_point(size = 3)
p1<-p+mytheme+ylim(0,0.0155)
p1
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/PBSIM_2L_snp_den_Unsliding.pdf", sep = ""), plot = p1, width = 30, height = 3)


##TE read recall
x <- matrix(c(219,	86,	615,	400), ncol = 2)
x <- matrix(c(210,	98,	591,	495), ncol = 2)
w_TE_detect=219+210; w_TE_undetect=86+98
wo_TE_detect=615+591; wo_TE_undetect=400+495
w_TE_detect/(w_TE_detect+w_TE_undetect)
wo_TE_detect/(wo_TE_detect+wo_TE_undetect)

x <- matrix(c(wo_TE_detect,	w_TE_undetect,	wo_TE_detect,	wo_TE_undetect), ncol = 2)
chisq.test(x)
fisher.test(x)

##SV read recall
x <- matrix(c(155,	44,	679,	442), ncol = 2)
x <- matrix(c(138,	73,	663,	520), ncol = 2)
w_TE_detect=155+138; w_TE_undetect=44+73
wo_TE_detect=679+663; wo_TE_undetect=442+520

x <- matrix(c(425,	163,	409,	323), ncol = 2)
x <- matrix(c(419,	182,	382,	411), ncol = 2)
w_TE_detect=425+419; w_TE_undetect=163+182
wo_TE_detect=409+382; wo_TE_undetect=323+411
w_TE_detect/(w_TE_detect+w_TE_undetect)
wo_TE_detect/(wo_TE_detect+wo_TE_undetect)
########

#### for Fig. S5 A
repeat_samp=c();un_samp=c()
#effective_depth=775;total_indiv=7166
#effective_depth=966;total_indiv=7141
#effective_depth=578;total_indiv=192
effective_depth=1000;total_indiv=1000


for (k in 1:95){
  test_indiv=total_indiv*k*0.1
  mat<-rmultinom(n=1000, size=effective_depth, prob=rep(1/test_indiv,test_indiv))
  mylist=c();myUnlist=c()
  for (i in 1:ncol(mat)){
    mySample=0; myUnsample=0;
    for (j in 1:nrow(mat)){
      if (mat[j,i]>1){
        extra=mat[j,i]-1
        mySample=mySample+extra
      }
      if (mat[j,i]==0){
        myUnsample=myUnsample+1
      }
    }
    mylist=c(mylist,mySample)
    myUnlist=c(myUnlist,myUnsample)
  }
  repeat_samp=c(repeat_samp,mean(mylist)/effective_depth)
  un_samp=c(un_samp,mean(myUnlist)/test_indiv)
}
repeat_samp;un_samp

#my_resamp=(myA6_resamp+myA7_resamp)/2
ind_depth_ratio=seq(0.1,9.5,by = 0.1)
depth_to_ind_ratio=1/ind_depth_ratio
#+scale_color_manual(values=c("#CC79A7", "#56B4E9"))
subsample<-data.frame(ind_depth_ratio,depth_to_ind_ratio,repeat_samp,un_samp)
data<-subsample[subsample$depth_to_ind_ratio<8,]
p<-ggplot(data, aes(x = depth_to_ind_ratio, y = repeat_samp)) + ylab("proportion of repeatedly sampled individuals")+ xlab("depth per haplotype")+geom_point(size = 0.5)+ geom_smooth(method = "loess",span = 0.1,color = "black",size =1.5)
#p<-p#+geom_vline(xintercept=0.332, linetype="dashed", color = "black", size=1.5)+geom_vline(xintercept=9.25, linetype="dashed", color = "#CC79A7", size=1.5)+geom_vline(xintercept=7.39, linetype="dashed", color = "#56B4E9", size=1.5)
p<-p+geom_vline(xintercept=3, linetype="dashed", color = "#e69f00", size=1.5)+geom_vline(xintercept=0.09, linetype="dashed", color = "#CC79A7", size=1.5)+geom_vline(xintercept=0.16, linetype="dashed", color = "#56B4E9", size=1.5)
p1<-p+mytheme
p1
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/proportion_resampled_depth_to_haplotype_ratio_benchmark_A6_A7.pdf", sep = ""), plot = p1, width = 5, height = 4)

######
