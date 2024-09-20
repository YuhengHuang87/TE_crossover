#######CO number

###TE window pairs vs control pairs
rec<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A6_crossover_num_TE_windows_left_right_window_magnitude_focal_homolog.txt",header=F)
foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V13),'/',fixed=TRUE)));rec$class <- foo$X1;rec<-rec[rec$class != "Unknown", ];rec$TE_type<-ifelse(rec$class == "DNA","DNA","RNA")
#rec<-rec[rec$TE_type == "RNA", ]
#rec<-rec[rec$V12 != "roo",]
rec<-rec[rec$V5 >= 5 & rec$V7 >= 5 & rec$V24 >=5 & rec$V26 >=5, ]
#rec<-rec[rec$V12 != "nearby_ambiguous",]
rec$K9_mag_left<-(rec$V14+rec$V15+rec$V16)/3;rec$K9_mag_right<-(rec$V17+rec$V18+rec$V19)/3;rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2
rec$K9_alt_left<-(rec$V31+rec$V32+rec$V33)/3;rec$K9_alt_right<-(rec$V34+rec$V35+rec$V36)/3;rec$alt_K9<-(rec$K9_alt_left+rec$K9_alt_right)/2
rec$nor_K9=rec$K9/rec$alt_K9
rec$co<-rec$V8+rec$V10;rec$cov<-(rec$V4+rec$V6)/2;rec$alt_co<-rec$V27+rec$V29
rec$snp<-rec$V5+rec$V7;rec$TE <- paste(rec$V1, rec$V2, sep="_");rec$leng = rec$V3-rec$V2+1
#roo_en<-rec[rec$V12 =="roo" & rec$co > 0, ];rec<-rec[!(rec$TE %in% roo_en$TE),]
rec$alt_snp<-rec$V24+rec$V26;
rec$foc.cov=round((rec$V4+rec$V6)/4);rec$alt.cov=round((rec$V23+rec$V25)/4);
rec<-rec[rec$leng>=200,]
K9_cutoff=1
#common_enrich<-rec[rec$alt_K9>1 & rec$K9>1,];rec<-rec[!(rec$TE %in% common_enrich$TE),]
rec_K9<-rec[rec$K9 > K9_cutoff, ];
rec_nonK9<-rec[rec$K9 < K9_cutoff, ];

rec<-rec[rec$K9 < K9_cutoff & rec$alt_K9 < K9_cutoff, ];
enrich<-rec[rec$K9>1,];nonenrich<-rec[rec$K9<1,];wilcox.test(enrich$snp/enrich$alt_snp,nonenrich$snp/nonenrich$alt_snp, exact=FALSE, correct=FALSE, paired = F)
rec$ratio_snp=rec$snp/rec$alt_snp
wilcox.test(rec_K9$ratio_snp,rec_nonK9$ratio_snp, exact=FALSE, correct=FALSE, paired = F);
wilcox.test(rec$snp,rec$alt_snp, exact=FALSE, correct=FALSE, paired = T);rec$dif_snp=rec$snp-rec$alt_snp

control_1<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7with_rec_local_control_window_location_include_nearby_ambig_10000_30kb",header=F)
control_2<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/control_focal_A7_alter_A6_TE_height_HMD_25_5000_50_mean_include_ambiguous.txt",header=F)
control_1$K9=(control_1$V19+control_1$V20+control_1$V21+control_1$V22+control_1$V23+control_1$V24)/6 #switch to match with A6 (start with 23) 
control_2$K9=(control_2$V19+control_2$V20+control_2$V21+control_2$V23+control_2$V24+control_2$V25)/6

#control_1$K9=(control_1$V19+control_1$V20+control_1$V21+control_1$V22+control_1$V23+control_1$V24)/6 #switch to match with A7 (start with 22)
#control_2$K9=(control_2$V19+control_2$V20+control_2$V21+control_2$V23+control_2$V24+control_2$V25)/6

control_1<-control_1[control_1$K9 < 1, ];control_2<-control_2[control_2$K9 < 1, ]
control<-control_1[control_1$V2 %in% intersect(control_1$V2, control_2$V2),];control$win <- paste(control$V1, control$V2, sep="_")

control$foc.cov = round(control$V4/2); control$alt.cov = round(control$V9/2);control$dif_cov = control$foc.cov-control$alt.cov;control$ratio_cov = control$foc.cov/control$alt.cov
control$foc_snp = control$V5;control$alt_snp = control$V10;control$dif_snp=control$foc_snp-control$alt_snp

control$co=control$V11+control$V13;control$alt_co=control$V15+control$V17;
control$ratio_snp=control$foc_snp/control$alt_snp;

a=seq(0.1,1,0.1);bin<-c(0,quantile(rec$ratio_snp,probs=a))
bin_counts<-hist(control$ratio_snp,breaks = c(bin,Inf), plot = FALSE)$counts[1:10]
hist(control$ratio_snp,breaks = c(bin,Inf), plot = FALSE)$counts
hist(rec$ratio_snp,breaks = c(bin,Inf), plot = FALSE)$counts

var_unique <- data.frame()
#cont<-cont[cont$ratio_snp<5.2,]
for (j in 1:10){
  sub_set<-control[control$ratio_snp>bin[j] & control$ratio_snp<bin[j+1],]
  sub_sample<- sub_set[sample(1:nrow(sub_set), min(bin_counts), replace=F),];
  var_unique<-rbind(var_unique,sub_sample)
}
cont<-var_unique

wilcox.test(rec$snp,rec$alt_snp, exact=FALSE, correct=FALSE, paired = T);mean(rec$snp);mean(cont$alt_snp)
wilcox.test(cont$foc_snp,cont$alt_snp, exact=FALSE, correct=FALSE, paired = T);mean(cont$foc_snp);mean(cont$alt_snp)
wilcox.test(rec$ratio_snp,var_unique$ratio_snp, exact=FALSE, correct=FALSE, paired = F);

wilcox.test(rec$dif_snp,cont$dif_snp, exact=FALSE, correct=FALSE, paired = F);mean(rec$dif_snp);mean(cont$dif_snp)
wilcox.test(rec$snp,cont$foc_snp, exact=FALSE, correct=FALSE, paired = F)
wilcox.test(rec$alt_snp,cont$alt_snp, exact=FALSE, correct=FALSE, paired = F)

rec_as<-data.frame(rec$V1,rec$V2, rec$co,rec$alt_co, rec$foc.cov, rec$alt.cov)

con_as<-data.frame(cont$V1,cont$V2,cont$co,cont$alt_co,cont$foc.cov,cont$alt.cov)

### You might save the subsampled control windows
write.table(con_as, file="/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A6_control_CO_number_betw_strain_K9_SNP_ratio_filtered.txt",sep="\t",eol="\n",row.names = F,col.names = F)
con_as<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A6_control_CO_number_betw_strain_K9_SNP_ratio_filtered.txt",header=F)
colnames(con_as) <- c("cont.V1", "cont.V2","cont.co","cont.alt_co","cont.foc.cov","cont.alt.cov")
###

co_sam=c();alt_sam=c();co_sam_con=c();alt_sam_con=c()
dif_co_sam=c();dif_co_sam_K9=c();dif_co_sam_nonK9=c();
dif_cont_sam=c()
p_MWU_TE_K9=c();
p_MWU_TE_K9_con=c();p_MWU_TE_nonK9_con=c()
p_MWU_TE_cont=c();
rho=c();p_cor=c();

for (j in 1:1000){
  for (i in 1:nrow(rec_as)){
    if (rec_as$rec.foc.cov[i] > rec_as$rec.alt.cov[i]){
      if (rec_as$rec.co[i] > 0){
        s_list<-sample.int(rec_as$rec.foc.cov[i], size = rec_as$rec.alt.cov[i], replace = FALSE)
        rec_as$sam_rec_co[i] = length(s_list[s_list<=rec_as$rec.co[i]])
        rec_as$sam_cov[i] = rec_as$rec.alt.cov[i]
        rec_as$sam_alt_co[i] = rec_as$rec.alt_co[i]
        rec_as$sam_cov[i] = rec_as$rec.alt.cov[i]
      }else{
        rec_as$sam_rec_co[i] = rec_as$rec.co[i]
        rec_as$sam_cov[i] = rec_as$rec.alt.cov[i]
        rec_as$sam_alt_co[i] = rec_as$rec.alt_co[i]
        rec_as$sam_cov[i] = rec_as$rec.alt.cov[i]
      }
    }else{
      if (rec_as$rec.alt_co[i] > 0){
        s_list<-sample.int(rec_as$rec.alt.cov[i], size =rec_as$rec.foc.cov[i] , replace = FALSE)
        rec_as$sam_rec_co[i] = rec_as$rec.co[i]
        rec_as$sam_cov[i] = rec_as$rec.foc.cov[i]
        rec_as$sam_alt_co[i] = length(s_list[s_list<=rec_as$rec.alt_co[i]])
        rec_as$sam_cov[i] = rec_as$rec.foc.cov[i]
      }else{
        rec_as$sam_rec_co[i] = rec_as$rec.co[i]
        rec_as$sam_cov[i] = rec_as$rec.foc.cov[i]
        rec_as$sam_alt_co[i] = rec_as$rec.alt_co[i]
        rec_as$sam_cov[i] = rec_as$rec.foc.cov[i]
      }
    }}
  rec_as$sam_co_dif = rec_as$sam_rec_co-rec_as$sam_alt_co;
  dif_co_sam=c(dif_co_sam, mean(rec_as$sam_co_dif));
  co_sam=c(co_sam,mean(rec_as$sam_rec_co));alt_sam=c(alt_sam,mean(rec_as$sam_alt_co))
  
  
  for (i in 1:nrow(con_as)){
    if (con_as$cont.foc.cov[i] > con_as$cont.alt.cov[i]){
      if (con_as$cont.co[i] > 0){
        s_list<-sample.int(con_as$cont.foc.cov[i], size = con_as$cont.alt.cov[i], replace = FALSE)
        con_as$focal_co[i] = length(s_list[s_list<=con_as$cont.co[i]])
        con_as$focal_cov[i] = con_as$cont.alt.cov[i]
        con_as$alter_co[i] = con_as$cont.alt_co[i]
        con_as$alter_cov[i] = con_as$cont.alt.cov[i]
      }else{
        con_as$focal_co[i] = con_as$cont.co[i]
        con_as$focal_cov[i] = con_as$cont.alt.cov[i]
        con_as$alter_co[i] = con_as$cont.alt_co[i]
        con_as$alter_cov[i] = con_as$cont.alt.cov[i]
      }
    }else{
      if (con_as$cont.alt_co[i] > 0){
        s_list<-sample.int(con_as$cont.alt.cov[i], size =con_as$cont.foc.cov[i] , replace = FALSE)
        con_as$focal_co[i] = con_as$cont.co[i]
        con_as$focal_cov[i] = con_as$cont.foc.cov[i]
        con_as$alter_co[i] = length(s_list[s_list<=con_as$cont.alt_co[i]])
        con_as$alter_cov[i] = con_as$cont.foc.cov[i]
      }else{
        con_as$focal_co[i] = con_as$cont.co[i]
        con_as$focal_cov[i] = con_as$cont.foc.cov[i]
        con_as$alter_co[i] = con_as$cont.alt_co[i]
        con_as$alter_cov[i] = con_as$cont.foc.cov[i]
      }
    }}
  con_as$sam_dif = con_as$focal_co-con_as$alter_co;
  co_sam_con=c(co_sam_con,mean(con_as$focal_co));alt_sam_con=c(alt_sam_con,mean(con_as$alter_co))
  
  dif_cont_sam=c(dif_cont_sam, mean(con_as$sam_dif));
  p.value_mwu<-wilcox.test(rec_as$sam_co_dif,con_as$sam_dif, alternative = "less", exact=FALSE, correct=FALSE, paired = F)$p.value;#median(rec$sam_dif);median(cont$sam_dif)
  #p.K9_mwu<-wilcox.test(rec_K9_TE$sam_co_dif,con_as$sam_dif, alternative = "less", exact=FALSE, correct=FALSE, paired = F)$p.value;#median(rec$sam_dif);median(cont$sam_dif)
  #p.nonK9_mwu<-wilcox.test(rec_nonK9_TE$sam_co_dif,con_as$sam_dif, alternative = "less", exact=FALSE, correct=FALSE, paired = F)$p.value;#median(rec$sam_dif);median(cont$sam_dif)
  p_MWU_TE_cont=c(p_MWU_TE_cont,p.value_mwu);
  #p_MWU_TE_K9_con=c(p_MWU_TE_K9_con,p.K9_mwu);p_MWU_TE_nonK9_con=c(p_MWU_TE_nonK9_con,p.nonK9_mwu)
}
mean(co_sam);mean(alt_sam)
mean(co_sam_con);mean(alt_sam_con)

median(dif_co_sam);median(dif_cont_sam)
median(p_MWU_TE_cont);

p_MWU=as.data.frame(p_MWU_TE_cont)
p<-ggplot(p_MWU, aes(x=p_MWU_TE_cont)) + scale_x_continuous()+
  geom_histogram(color="black", fill="#CC79A7")+xlab("Mann-Whitney U test p-value")+ylab("counts")
p<-p+mytheme
p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/histogram_TE_control_MWU_pvalue_A6.pdf", sep = ""), plot = p, width = 7, height = 6)

cols <- c("grey","#af3d81","#d993bc")
cols <- c("grey","#346db2","#9cbce2")
TE_K9<-data.frame(dif_co_sam_K9); TE_nonK9<-data.frame(dif_co_sam_nonK9); control_win<-data.frame(dif_cont_sam)
TE_K9$group<-"TE_K9";TE_nonK9$group<-"TE_nonK9";control_win$group<-"control"
data<-rbind(TE_K9, setNames(TE_nonK9, names(TE_K9)))
dat<-rbind(data,setNames(control_win,names(data)))

cols <- c("grey","#CC79A7")
cols <- c("grey","#56B4E9")
TE<-data.frame(dif_co_sam); control_win<-data.frame(dif_cont_sam)
TE$group<-"TE";control_win$group<-"control"
dat<-rbind(TE,setNames(control_win,names(TE)))

colnames(dat)[1] <- "Dif_CO"
p<-ggplot(dat, aes(x=Dif_CO, colour=group)) + geom_density(linewidth=2,bw=0.008)+ scale_color_manual(values = cols)+xlim(-0.025,0.0605) #+xlim(-0.16,0.08)
#p<-ggplot(dat, aes(x=Dif_CO, fill=group)) + geom_histogram(alpha=.2,binwidth=0.01)+ scale_fill_manual(values = cols)#+xlim(-0.025,0.055) 
p<-p+mytheme_black; p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/A6_CO_number_density_betw_strains_TEVscontrol_K9_SNPratio_filtered_black.pdf", sep = ""), plot = p, width = 7, height = 5)
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/A7_CO_number_density_betw_strains_TE_K9VsNonK9.pdf", sep = ""), plot = p, width = 7, height = 5)


####TE window pairs with and without K9
rec<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7_crossover_num_TE_windows_left_right_window_magnitude_focal_homolog.txt",header=F)
foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V13),'/',fixed=TRUE)));rec$class <- foo$X1;rec<-rec[rec$class != "Unknown", ];rec$TE_type<-ifelse(rec$class == "DNA","DNA","RNA")
#rec<-rec[rec$TE_type == "RNA", ]
#rec<-rec[rec$class == "LTR", ]
#rec<-rec[rec$V12 != "roo",]
rec<-rec[rec$V5 >= 5 & rec$V7 >= 5 & rec$V24 >=5 & rec$V26 >=5, ]

rec$K9_mag_left<-(rec$V14+rec$V15+rec$V16)/3;rec$K9_mag_right<-(rec$V17+rec$V18+rec$V19)/3;rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2
rec$K9_alt_left<-(rec$V31+rec$V32+rec$V33)/3;rec$K9_alt_right<-(rec$V34+rec$V35+rec$V36)/3;rec$alt_K9<-(rec$K9_alt_left+rec$K9_alt_right)/2
rec$nor_K9=rec$K9/rec$alt_K9
rec$co<-rec$V8+rec$V10;rec$cov<-(rec$V4+rec$V6)/2;rec$alt_co<-rec$V27+rec$V29
rec$snp<-rec$V5+rec$V7;rec$TE <- paste(rec$V1, rec$V2, sep="_");rec$leng = rec$V3-rec$V2+1
#roo_en<-rec[rec$V12 =="roo" & rec$co > 0, ];rec<-rec[!(rec$TE %in% roo_en$TE),]
rec$alt_snp<-rec$V24+rec$V26;
rec$foc.cov=round((rec$V4+rec$V6)/4);rec$alt.cov=round((rec$V23+rec$V25)/4);
rec<-rec[rec$leng>=200,]
#rec<-rec[rec$V12 != "nearby_ambiguous",]
rec$TE_fam<-ifelse(rec$V12 == "nearby_ambiguous","TE_complex","simple")

#cluster<-rec[rec$V12 == "nearby_ambiguous" & rec$leng > 10000,];rec<-rec[!(rec$TE %in% cluster$TE),]

K9_cutoff=1
rec<-rec[rec$alt_K9<K9_cutoff,]#for studying cor K9 and CO
#common_enrich<-rec[rec$alt_K9>1 & rec$K9>1,];rec<-rec[!(rec$TE %in% common_enrich$TE),]
#rec<-rec[rec$K9<K9_cutoff & rec$alt_K9<K9_cutoff,]#for studying cor length and CO

rec_K9<-rec[rec$K9 > K9_cutoff, ];
rec_nonK9<-rec[rec$K9 < K9_cutoff, ];
#alter_K9<-rec[rec$alt_K9>1,]
#enrich<-rec[rec$K9>1,];nonenrich<-rec[rec$K9<1,];wilcox.test(enrich$snp/enrich$alt_snp,nonenrich$snp/nonenrich$alt_snp, exact=FALSE, correct=FALSE, paired = F)
#wilcox.test(rec$snp,rec$alt_snp, exact=FALSE, correct=FALSE, paired = T);rec$dif_snp=rec$snp-rec$alt_snp
LTR<-rec[rec$class == "LTR",];mean(LTR$K9);var(LTR$K9)
LINE<-rec[rec$class == "LINE",];mean(LINE$K9);var(LINE$K9)
DNA<-rec[rec$class == "DNA",];mean(DNA$K9);var(DNA$K9)

rec_as<-data.frame(rec$V1,rec$V2, rec$co,rec$alt_co, rec$foc.cov, rec$alt.cov, rec$K9, rec$alt_K9,rec$leng,rec$class,rec$TE_type,rec$TE_fam)
#rec_as<-data.frame(rec_K9_TE$V1,rec_K9_TE$V2, rec_K9_TE$co,rec_K9_TE$alt_co, rec_K9_TE$foc.cov, rec_K9_TE$alt.cov)
#rec_as<-data.frame(rec_nonK9_TE$V1,rec_nonK9_TE$V2, rec_nonK9_TE$co,rec_nonK9_TE$alt_co, rec_nonK9_TE$foc.cov, rec_nonK9_TE$alt.cov)
co_sam=c();alt_sam=c();co_sam_con=c();alt_sam_con=c()
dif_co_sam=c();dif_co_sam_K9=c();dif_co_sam_nonK9=c();
dif_cont_sam=c()
p_MWU_TE_K9=c();
p_MWU_TE_K9_con=c();p_MWU_TE_nonK9_con=c()
p_MWU_TE_cont=c();
rho=c();p_cor=c();

for (j in 1:1000){
  for (i in 1:nrow(rec_as)){
    if (rec_as$rec.foc.cov[i] > rec_as$rec.alt.cov[i]){
      if (rec_as$rec.co[i] > 0){
        s_list<-sample.int(rec_as$rec.foc.cov[i], size = rec_as$rec.alt.cov[i], replace = FALSE)
        rec_as$sam_rec_co[i] = length(s_list[s_list<=rec_as$rec.co[i]])
        rec_as$sam_cov[i] = rec_as$rec.alt.cov[i]
        rec_as$sam_alt_co[i] = rec_as$rec.alt_co[i]
        rec_as$sam_cov[i] = rec_as$rec.alt.cov[i]
      }else{
        rec_as$sam_rec_co[i] = rec_as$rec.co[i]
        rec_as$sam_cov[i] = rec_as$rec.alt.cov[i]
        rec_as$sam_alt_co[i] = rec_as$rec.alt_co[i]
        rec_as$sam_cov[i] = rec_as$rec.alt.cov[i]
      }
    }else{
      if (rec_as$rec.alt_co[i] > 0){
        s_list<-sample.int(rec_as$rec.alt.cov[i], size =rec_as$rec.foc.cov[i] , replace = FALSE)
        rec_as$sam_rec_co[i] = rec_as$rec.co[i]
        rec_as$sam_cov[i] = rec_as$rec.foc.cov[i]
        rec_as$sam_alt_co[i] = length(s_list[s_list<=rec_as$rec.alt_co[i]])
        rec_as$sam_cov[i] = rec_as$rec.foc.cov[i]
      }else{
        rec_as$sam_rec_co[i] = rec_as$rec.co[i]
        rec_as$sam_cov[i] = rec_as$rec.foc.cov[i]
        rec_as$sam_alt_co[i] = rec_as$rec.alt_co[i]
        rec_as$sam_cov[i] = rec_as$rec.foc.cov[i]
      }
    }}
  rec_as$sam_co_dif = rec_as$sam_rec_co-rec_as$sam_alt_co;
  dif_co_sam=c(dif_co_sam, mean(rec_as$sam_co_dif));
  co_sam=c(co_sam,mean(rec_as$sam_rec_co));alt_sam=c(alt_sam,mean(rec_as$sam_alt_co))
  
  rec_K9_TE<-rec_as[rec_as$rec.K9 > K9_cutoff, ];rec_nonK9_TE<-rec_as[rec_as$rec.K9 < K9_cutoff, ];
  #rec_K9_TE<-rec_as[rec_as$rec.TE_type == "RNA", ];rec_nonK9_TE<-rec_as[rec_as$rec.TE_type == "DNA", ];
  #rec_K9_TE<-rec_as[rec_as$rec.class == "LTR", ];rec_nonK9_TE<-rec_as[rec_as$rec.class == "DNA", ];
  
  dif_co_sam_K9=c(dif_co_sam_K9, mean(rec_K9_TE$sam_co_dif));
  dif_co_sam_nonK9=c(dif_co_sam_nonK9, mean(rec_nonK9_TE$sam_co_dif));
  p_mwu<-wilcox.test(rec_K9_TE$sam_co_dif,rec_nonK9_TE$sam_co_dif,alternative = "less", exact=FALSE, correct=FALSE, paired = F)$p.value;#mean(rec_K9_TE$sam_dif);mean(rec_nonK9_TE$sam_dif)
  p_MWU_TE_K9=c(p_MWU_TE_K9,p_mwu)
  rho=c(rho,cor.test(rec_as$rec.K9,rec_as$sam_co_dif,method = "spearman")$estimate);p_cor=c(p_cor,cor.test(rec_as$rec.K9,rec_as$sam_co_dif,alternative = "less",method = "spearman")$p.value);
  #rho=c(rho,cor.test(rec_as$rec.nor_K9,rec_as$sam_co_dif,method = "spearman")$estimate);p_cor=c(p_cor,cor.test(rec_as$rec.nor_K9,rec_as$sam_co_dif,alternative = "less",method = "spearman")$p.value);
  #rho=c(rho,cor.test(rec_as$rec.leng,rec_as$sam_co_dif,method = "spearman")$estimate);p_cor=c(p_cor,cor.test(rec_as$rec.leng,rec_as$sam_co_dif,alternative = "less",method = "spearman")$p.value);
}
median(dif_co_sam_K9);median(dif_co_sam_nonK9)
median(p_MWU_TE_K9)

median(rho);median(p_cor)
nrow(rec_K9_TE);nrow(rec_nonK9_TE)

cols <- c("#af3d81","#d993bc")
cols <- c("#346db2","#9cbce2")
TE_K9<-data.frame(dif_co_sam_K9); TE_nonK9<-data.frame(dif_co_sam_nonK9); control_win<-data.frame(dif_cont_sam)
TE_K9$group<-"TE_K9";TE_nonK9$group<-"TE_nonK9";
dat<-rbind(TE_K9, setNames(TE_nonK9, names(TE_K9)))

colnames(dat)[1] <- "Dif_CO"
p<-ggplot(dat, aes(x=Dif_CO, colour=group)) + geom_density(linewidth=2,bw=0.008)+ scale_color_manual(values = cols)+xlim(-0.21,-0.06)#+xlim(-0.16,0.07)
#p<-ggplot(dat, aes(x=Dif_CO, fill=group)) + geom_histogram(alpha=.2,binwidth=0.01)+ scale_fill_manual(values = cols)#+xlim(-0.025,0.055) 
p<-p+mytheme; p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/A7_CO_number_density_betw_strains_TE_K9VsNonK9.pdf", sep = ""), plot = p, width = 7, height = 5)

cols <- c("grey","#CC79A7")
cols <- c("grey","#56B4E9")
p_MWU=as.data.frame(p_MWU_TE_K9)
p<-ggplot(p_MWU, aes(x=p_MWU_TE_K9)) + scale_x_continuous()+
  geom_histogram(color="black", fill="#56B4E9")+xlab("Mann-Whitney U test p-value")+ylab("counts")+xlim(0.15,0.75)
p<-p+mytheme
p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/histogram_TE_K9_MWU_pvalue_A7.pdf", sep = ""), plot = p, width = 7, height = 6)

p_rho=as.data.frame(rho)
p<-ggplot(p_rho, aes(x=rho)) + scale_x_continuous()+
  geom_histogram(color="black", fill="#CC79A7",binwidth=0.02)+xlab("spearman correlation coefficient")+ylab("counts")
p<-p+mytheme
p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/histogram_TE_length_rho_A6.pdf", sep = ""), plot = p, width = 7, height = 6)

p1<-ggplot(aes(x = rec.leng, y = sam_co_dif,shape = rec.TE_fam,colour= rec.TE_fam), data = rec_as) + geom_point(size=4,alpha = 0.5)+ xlab("TE length")+ylab("Dif_CO") #theme(axis.line.x = element_line(size = 1), axis.line.y = element_line(size = 1),axis.title.x = element_text(size=18),axis.title.y = element_text(angle=90,size=18), axis.text.x  = element_text(vjust=1, size=18,color="white"),axis.text.y = element_text(hjust=1, size=18,color="white"), legend.text = element_text(hjust=20,size=20,face="italic"),legend.title = element_text(hjust=0,size=20,face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())
p1<-p1+mytheme;p1
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/Dif_CO_distance_one_downsampling_scatter_TE_length_simple_complex_A7.pdf", sep = ""), plot = p1, width = 8, height = 6)
#####

##### Distance to the nearest COs
### compare between window pairs with and without TEs
rec<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7_crossover_nearest_distance_coverage_SNP_TE_left_right_window_mass_5000_exclude_distance_5000",header=F)
foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V5),'/',fixed=TRUE)));rec$class <- foo$X1;rec<-rec[rec$class != "Unknown", ];rec$TE_type<-ifelse(rec$class == "DNA","DNA","RNA")
#rec<-rec[rec$V1 != "roo", ]

rec$leng = rec$V4-rec$V3+1;
rec$snp_den=rec$V9/rec$V6;rec$con_snp_den=rec$V22/rec$V19;
rec$ratio_snp=rec$snp_den/rec$con_snp_den

rec$K9_mag_left<-abs(rec$V10+rec$V11+rec$V12)/3
rec$K9_mag_right<-abs(rec$V13+rec$V14+rec$V15)/3
#rec$K9_mag_left<-abs(rec$V30+rec$V31+rec$V32)/3 #for filtering out overlapping region in alternative strains
#rec$K9_mag_right<-abs(rec$V33+rec$V34+rec$V35)/3
rec$TE <- paste(rec$V2, rec$V3, sep="_");

rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2
rec<-rec[rec$leng>=200,]

rec$K9_alter_left<-abs(rec$V23+rec$V24+rec$V25)/3
rec$K9_alter_right<-abs(rec$V26+rec$V27+rec$V28)/3
rec$K9_alter<-(rec$K9_alter_left+rec$K9_alter_right)/2
rec$foc.cov<-round(rec$V8/2);rec$alt.cov<-round(rec$V21/2)

#median(rec$K9);median(rec$K9_alter)
rec<-rec[rec$K9_alter< 3,]#need to change if comparing different TEs or with control
#common_enrich<-rec[rec$K9_alter>3 & rec$K9>3,];rec<-rec[!(rec$TE %in% common_enrich$TE),]

#alter_K9<-rec[rec$K9_alter>3,]
K9_cutoff=median(rec$K9)
rec<-rec[rec$K9< K9_cutoff,]
#K9_cutoff=3
LTR<-rec[rec$class == "LTR",]
LINE<-rec[rec$class == "LINE",]
DNA<-rec[rec$class == "DNA",]

control_1<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/control_focal_K9_magnitude_A7_alter_A6_nearest_CO_distance_mean_Mass_TE_HMD_25_1000_10_measures_include_ambiguous.txt",header=F)
control_2<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/control_focal_A7_alter_K9_magnitude_A6_nearest_CO_distance_TE_HMD_25_1000_10_measures_include_ambiguous.txt",header=F)
control_1$K9=(control_1$V15+control_1$V16+control_1$V17+control_1$V18+control_1$V19+control_1$V20)/6
control_2$K9=(control_2$V15+control_2$V16+control_2$V17+control_2$V18+control_2$V19+control_2$V20)/6
control_1<-control_1[control_1$K9 < 1, ];control_2<-control_2[control_2$K9 < 1, ]
control<-control_1[control_1$V2 %in% intersect(control_1$V2, control_2$V2),];control$win <- paste(control$V1, control$V2, sep="_")

#con<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/control_A6_alter_A7_nearest_second_distanc_exclude_TE_overlapped_average_10000.txt",header=F)
#con$win <- paste(con$V1, con$V2, sep="_");con$alt_win <-paste(con$V8, con$V9, sep="_")
#cont<-con
#cont<-con[con$win %in% intersect(con$win, control$win),]

control$foc.cov = control$V6/2; control$alt.cov = control$V13/2;
control$foc_cov = control$V6/2; control$alter_cov = control$V13/2;
control$snp_den = control$V7/control$V4;control$alt_snp_den = control$V14/control$V11;control$ratio_snp=control$snp_den/control$alt_snp_den;
wilcox.test(rec$ratio_snp,control$ratio_snp, exact=FALSE, correct=FALSE, paired = F);median(rec$ratio_snp);median(control$ratio_snp)

a=seq(0.1,1,0.1);bin<-c(0,quantile(rec$ratio_snp,probs=a))
#high_ratio_snp<-rec[rec$ratio_snp>5,]
hist(control$ratio_snp,breaks = c(bin,Inf), plot = FALSE);#hist(rec$ratio_snp,breaks = c(bin,Inf), plot = FALSE)
bin_counts<-hist(control$ratio_snp,breaks = c(bin,Inf), plot = FALSE)$counts[1:10]

var_unique <- data.frame()
for (j in 1:10){
  sub_set<-control[control$ratio_snp>bin[j] & control$ratio_snp<bin[j+1],]
  sub_sample<- sub_set[sample(1:nrow(sub_set), min(bin_counts), replace=F),];#140 A6; 131 A7
  var_unique<-rbind(var_unique,sub_sample)
}
wilcox.test(rec$ratio_snp,var_unique$ratio_snp, exact=FALSE, correct=FALSE, paired = F);

cont<-var_unique
con_as<-data.frame(cont$V1,cont$V4,cont$V5,cont$V11,cont$V12,cont$foc.cov,cont$alt.cov)

### You may save the subsampled control windows
write.table(con_as, file="/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A6_control_nearest_distance_betw_strain_K9_SNP_ratio_filtered.txt",sep="\t",eol="\n",row.names = F,col.names = T)
cont<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7_control_nearest_distance_betw_strain_K9_SNP_ratio_filtered.txt",header=T)
colnames(cont) <- c("V1", "V4","V5","V11","V12","foc_cov","alter_cov")
###

dif_rec_sam=c();dif_rec_sam_K9=c();dif_rec_sam_nonK9=c();
dif_cont_sam=c()
p_KS_TE_K9=c();p_MWU_TE_K9=c();p_t_TE_K9=c()
rho=c();p_cor=c();
p_KS_TE_cont=c();p_MWU_TE_cont=c();
p_MWU_TE_K9_con=c();p_MWU_TE_nonK9_con=c()

for (j in 1:1000){
  for (i in 1:nrow(rec)){
    if (rec$foc.cov[i] > rec$alt.cov[i]){
      s_list<-sample.int(rec$foc.cov[i], size = rec$alt.cov[i], replace = FALSE)
      if (1%in% s_list){
        rec$foc_dis[i] = rec$V6[i]
        rec$alt_dis[i] = rec$V19[i]
      }else{
        rec$foc_dis[i] = rec$V7[i]
        rec$alt_dis[i] = rec$V19[i]
      }
    }else{
      s_list<-sample.int(rec$alt.cov[i], size =rec$foc.cov[i] , replace = FALSE)
      if (1%in% s_list){
        rec$foc_dis[i] = rec$V6[i]
        rec$alt_dis[i] = rec$V19[i]
      }else{
        rec$foc_dis[i] = rec$V6[i]
        rec$alt_dis[i] = rec$V20[i]
      }
    }}
  rec$sam_dif = log2(rec$foc_dis/rec$alt_dis);
  #rec$sam_dif = rec$foc_dis/rec$alt_dis;
  
  #rec_K9_TE<-rec[rec$K9 > K9_cutoff, ];rec_nonK9_TE<-rec[rec$K9 < K9_cutoff, ]
  #rec_K9_TE<-rec[rec$dif_K9 > K9_cutoff, ];rec_nonK9_TE<-rec[rec$dif_K9 < K9_cutoff, ]
  #rec_K9_TE<-rec[rec$TE_type == "RNA", ];rec_nonK9_TE<-rec[rec$TE_type == "DNA", ];
  #rec_K9_TE<-rec[rec$class == "LTR", ];rec_nonK9_TE<-rec[rec$class == "DNA", ];
  
  dif_rec_sam=c(dif_rec_sam, median(rec$sam_dif));
  #rho=c(rho,cor.test(rec$K9,rec$sam_dif,method = "spearman")$estimate);p_cor=c(p_cor,cor.test(rec$K9,rec$sam_dif,alternative = "greater",method = "spearman")$p.value);
  #rho=c(rho,cor.test(rec$dif_K9,rec$sam_dif,method = "spearman")$estimate);p_cor=c(p_cor,cor.test(rec$dif_K9,rec$sam_dif,alternative = "greater",method = "spearman")$p.value);
  #rho=c(rho,cor.test(rec$leng,rec$sam_dif,method = "spearman")$estimate);p_cor=c(p_cor,cor.test(rec$leng,rec$sam_dif,alternative = "greater",method = "spearman")$p.value);
  
  
  for (i in 1:nrow(cont)){
    if (cont$foc_cov[i] > cont$alter_cov[i]){
      s_list<-sample.int(cont$foc_cov[i], size = cont$alter_cov[i], replace = FALSE)
      if (1%in% s_list){
        cont$foc_dis[i] = cont$V4[i]
        cont$alt_dis[i] = cont$V11[i]
      }else{
        cont$foc_dis[i] = cont$V5[i]
        cont$alt_dis[i] = cont$V11[i]
      }
    }else{
      s_list<-sample.int(cont$alter_cov[i], size =cont$foc_cov[i] , replace = FALSE)
      if (1%in% s_list){
        cont$foc_dis[i] = cont$V4[i]
        cont$alt_dis[i] = cont$V11[i]
      }else{
        cont$foc_dis[i] = cont$V4[i]
        cont$alt_dis[i] = cont$V12[i]
      }
    }}
  #cont$sam_dif = cont$foc_dis-cont$alt_dis;
  cont$sam_dif = log2(cont$foc_dis/cont$alt_dis);
  #cont$sam_dif = cont$foc_dis/cont$alt_dis;
  
  
  #cont$focal_dis = (cont$V8+cont$V10)/2; cont$alter_dis = (cont$V12 + cont$V14)/2;cont$dif_dis = (cont$focal_dis-cont$alter_dis);
  p.value_mwu<-wilcox.test(rec$sam_dif,cont$sam_dif, alternative = "greater", exact=FALSE, correct=FALSE, paired = F)$p.value;#median(rec$sam_dif);median(cont$sam_dif)
  p_MWU_TE_cont=c(p_MWU_TE_cont,p.value_mwu);
  dif_cont_sam=c(dif_cont_sam, median(cont$sam_dif));
}
median(dif_cont_sam);median(dif_rec_sam)
median(p_MWU_TE_cont);


p_MWU=as.data.frame(p_MWU_TE_cont)
p<-ggplot(p_MWU, aes(x=p_MWU_TE_cont)) + scale_x_continuous()+
  geom_histogram(color="black", fill="#CC79A7")+xlab("Mann-Whitney U test p-value")+ylab("counts")
p<-p+mytheme
p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/histogram_dif_log2_distance_CO_TE_control_MWU_pvalue_A6.pdf", sep = ""), plot = p, width = 7, height = 6)

cor.test(rec_focal$sam_co_dif,rec_focal$rec_s.K9,method="spearman")
cor.test(rec_nonK9_TE$rec_s.leng,rec_nonK9_TE$sam_co_dif,method="spearman")


cols <- c("grey","#af3d81","#d993bc")
cols <- c("grey","#346db2","#9cbce2")
#TE_K9<-data.frame(dif_co_K9); TE_nonK9<-data.frame(dif_co_nonK9); control_win<-data.frame(dif_co_control)
TE_K9<-data.frame(dif_rec_sam_K9); TE_nonK9<-data.frame(dif_rec_sam_nonK9);control_win<-data.frame(dif_cont_sam)
TE_K9$group<-"TE_K9";TE_nonK9$group<-"TE_non";control_win$group<-"control"
data<-rbind(TE_K9, setNames(TE_nonK9, names(TE_K9)))
dat<-rbind(data,setNames(control_win,names(data)))

#CC79A7 #56B4E9
cols <- c("grey","#CC79A7")
cols <- c("grey","#56B4E9")
TE<-data.frame(dif_rec_sam); control_win<-data.frame(dif_cont_sam)
#TE<-data.frame(rec$sam_dif); control_win<-data.frame(cont$sam_dif)
TE$group<-"TE";control_win$group<-"control"
dat<-rbind(TE,setNames(control_win,names(TE)))

cols <- c("grey","#d993bc")
cols <- c("grey","#9cbce2")
TE_nonK9<-data.frame(dif_rec_sam_nonK9); control_win<-data.frame(dif_cont_sam)
TE_nonK9$group<-"TE_nonK9";control_win$group<-"control"
dat<-rbind(TE_nonK9,setNames(control_win,names(TE_nonK9)))

cols <- c("grey")
control_win<-data.frame(dif_cont_sam);control_win$group<-"control"
dat<-control_win;


colnames(dat)[1] <- "log2_ratio_of_distance"
p<-ggplot(dat, aes(x=log2_ratio_of_distance, colour=group)) + geom_density(linewidth=2,bw=0.05)+ scale_color_manual(values = cols)#+xlim(-0.08,0.8)#+xlim(-0.8,0.4)##+xlim(-12000,8000)#
#p<-ggplot(dat, aes(x=log2_ratio_of_distance,colour = group,fill = group)) + geom_histogram(alpha=0.5, position="identity")+ scale_color_manual(values = cols)+ scale_fill_manual(values = cols)#+xlim(-0.025,0.055) alpha=.2,binwidth=0.01
p<-p+mytheme_black; p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/one_downsampling_A7_dif_nearest_distance_historm_TE_vs_control.pdf", sep = ""), plot = p, width = 8, height = 5)
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/one_downsampling_A7_dif_nearest_distance_historm_TE_K9orNot.pdf", sep = ""), plot = p, width = 8, height = 5)
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/density_A6_dif_nearest_distance_betw_strains_mean_K9orNot.pdf", sep = ""), plot = p, width = 8, height = 5)
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/distance_to_CO/density_A6_dif_nearest_distance_betw_strains_TE_control_black.pdf", sep = ""), plot = p, width = 7, height = 5)

#ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/density_A6_log2_ratio_of_nearest_distance_betw_strains_TE_control_alter_K9_nearby_excluded_5000.pdf", sep = ""), plot = p, width = 8, height = 5)
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/density_A6_log2_ratio_of_nearest_distance_betw_strains_TE_control_K9_both_allele_SNP_density.pdf", sep = ""), plot = p, width = 8, height = 5)
#ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/density_A6_log2_ratio_of_nearest_distance_betw_strains_TE_K9orNot_alter_K9_excluded_200bp.pdf", sep = ""), plot = p, width = 8, height = 5)


cont$group<-"con"; 
foc<-data.frame(rec_s$dif_dis,rec_s$group,rec_s$dif_depth,rec_s$ratio_depth,rec_s$dif_snp_den)
foc<-data.frame(rec_K9_TE$dif_dis,rec_K9_TE$group,rec_K9_TE$dif_depth,rec_K9_TE$ratio_depth,rec_K9_TE$dif_snp_den)
foc<-data.frame(rec_nonK9_TE$dif_dis,rec_nonK9_TE$group,rec_nonK9_TE$dif_depth,rec_nonK9_TE$ratio_depth,rec_nonK9_TE$dif_snp_den)
con<-data.frame(cont$dif_dis,cont$group,cont$dif_depth,cont$ratio_depth,cont$dif_snp_den)
total<-rbind(con, setNames(foc, names(con)))
m2 <- lm(log2(cont.dif_dis) ~ cont.group + cont.dif_depth + cont.dif_snp_den ,  data=total);summary(m2)

#### compare window pairs with high and low K9 mass
rec<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7_crossover_nearest_distance_coverage_SNP_TE_left_right_window_mass_5000_exclude_distance_5000",header=F)
#rec<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/focal_A6_alter_nearest_CO_distance_TE_mass_height_HMD_25_5000_50_mean_each_measures_include_ambiguous.txt",header=F)

foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V5),'/',fixed=TRUE)));rec$class <- foo$X1;rec<-rec[rec$class != "Unknown", ];rec$TE_type<-ifelse(rec$class == "DNA","DNA","RNA")
#rec<-rec[rec$V1 != "roo", ]
#rec<-rec[rec$TE_type == "RNA", ]
#rec<-rec[rec$class == "LTR", ]
rec$leng = rec$V4-rec$V3+1;
rec$snp_den=rec$V9/rec$V6;rec$con_snp_den=rec$V22/rec$V19;
rec$ratio_snp=rec$snp_den/rec$con_snp_den

rec$K9_mag_left<-abs(rec$V10+rec$V11+rec$V12)/3
rec$K9_mag_right<-abs(rec$V13+rec$V14+rec$V15)/3
rec$TE <- paste(rec$V2, rec$V3, sep="_");

rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2
rec<-rec[rec$leng>=200,]
#rec<-rec[rec$V1 != "nearby_ambiguous",]
#cluster<-rec[rec$V1 == "nearby_ambiguous" & rec$leng > 10000,];rec<-rec[!(rec$TE %in% cluster$TE),]

rec$K9_alter_left<-abs(rec$V23+rec$V24+rec$V25)/3
rec$K9_alter_right<-abs(rec$V26+rec$V27+rec$V28)/3
rec$K9_alter<-(rec$K9_alter_left+rec$K9_alter_right)/2
rec$foc.cov<-round(rec$V8/2);rec$alt.cov<-round(rec$V21/2)

rec<-rec[rec$K9_alter<3,]#for distance and K9
#rec<-rec[rec$K9<3 & rec$K9_alter<3,]#for distance and length
#common_enrich<-rec[rec$K9_alter>3 & rec$K9>3,];rec<-rec[!(rec$TE %in% common_enrich$TE),]
K9_cutoff=median(rec$K9)

#alter_K9<-rec[rec$K9_alter>3,]
#rec<-rec[rec$K9_alter<2 & rec$K9<2,]
rec$dif_K9=log2(rec$K9/rec$K9_alter)
#high_mass<-rec[rec$K9>50,]
#K9_cutoff=3

LTR<-rec[rec$class == "LTR",]
LINE<-rec[rec$class == "LINE",]
DNA<-rec[rec$class == "DNA",]

dif_rec_sam=c();dif_rec_sam_K9=c();dif_rec_sam_nonK9=c();
dif_cont_sam=c()
p_KS_TE_K9=c();p_MWU_TE_K9=c();p_t_TE_K9=c()
rho=c();p_cor=c();
p_KS_TE_cont=c();p_MWU_TE_cont=c();
p_MWU_TE_K9_con=c();p_MWU_TE_nonK9_con=c()

for (j in 1:1000){
  for (i in 1:nrow(rec)){
    if (rec$foc.cov[i] > rec$alt.cov[i]){
      s_list<-sample.int(rec$foc.cov[i], size = rec$alt.cov[i], replace = FALSE)
      if (1%in% s_list){
        rec$foc_dis[i] = rec$V6[i]
        rec$alt_dis[i] = rec$V19[i]
      }else{
        rec$foc_dis[i] = rec$V7[i]
        rec$alt_dis[i] = rec$V19[i]
      }
    }else{
      s_list<-sample.int(rec$alt.cov[i], size =rec$foc.cov[i] , replace = FALSE)
      if (1%in% s_list){
        rec$foc_dis[i] = rec$V6[i]
        rec$alt_dis[i] = rec$V19[i]
      }else{
        rec$foc_dis[i] = rec$V6[i]
        rec$alt_dis[i] = rec$V20[i]
      }
    }}
  rec$sam_dif = log2(rec$foc_dis/rec$alt_dis);#rec$sam_dif = rec$foc_dis/rec$alt_dis;
  
  rec_K9_TE<-rec[rec$K9 > K9_cutoff, ];rec_nonK9_TE<-rec[rec$K9 < K9_cutoff, ]
  #rec_K9_TE<-rec[rec$dif_K9 > K9_cutoff, ];rec_nonK9_TE<-rec[rec$dif_K9 < K9_cutoff, ]
  #rec_K9_TE<-rec[rec$TE_type == "RNA", ];rec_nonK9_TE<-rec[rec$TE_type == "DNA", ];
  #rec_K9_TE<-rec[rec$class == "LINE", ];rec_nonK9_TE<-rec[rec$class == "DNA", ];
  
  p_mwu<-wilcox.test(rec_K9_TE$sam_dif,rec_nonK9_TE$sam_dif,alternative = "greater", exact=FALSE, correct=FALSE, paired = F)$p.value;#mean(rec_K9_TE$sam_dif);mean(rec_nonK9_TE$sam_dif)
  #p_mwu<-wilcox.test(rec_K9_TE$sam_dif,rec_nonK9_TE$sam_dif, exact=FALSE, correct=FALSE, paired = F)$p.value;#mean(rec_K9_TE$sam_dif);mean(rec_nonK9_TE$sam_dif)
  p_MWU_TE_K9=c(p_MWU_TE_K9,p_mwu)
  dif_rec_sam=c(dif_rec_sam, median(rec$sam_dif));
  dif_rec_sam_K9=c(dif_rec_sam_K9, median(rec_K9_TE$sam_dif));dif_rec_sam_nonK9=c(dif_rec_sam_nonK9, median(rec_nonK9_TE$sam_dif));
  rho=c(rho,cor.test(rec$K9,rec$sam_dif,method = "spearman")$estimate);p_cor=c(p_cor,cor.test(rec$K9,rec$sam_dif,alternative = "greater",method = "spearman")$p.value);
  #rho=c(rho,cor.test(rec$dif_K9,rec$sam_dif,method = "spearman")$estimate);p_cor=c(p_cor,cor.test(rec$dif_K9,rec$sam_dif,alternative = "greater",method = "spearman")$p.value);
  #rho=c(rho,cor.test(rec$leng,rec$sam_dif,method = "spearman")$estimate);p_cor=c(p_cor,cor.test(rec$leng,rec$sam_dif,alternative = "greater",method = "spearman")$p.value);
}
median(dif_rec_sam_K9);median(dif_rec_sam_nonK9)
median(p_MWU_TE_K9);
median(rho);median(p_cor);
nrow(rec_K9_TE);nrow(rec_nonK9_TE)

cols <- c("#af3d81","#d993bc")
cols <- c("#346db2","#9cbce2")
TE_K9<-data.frame(dif_rec_sam_K9); TE_nonK9<-data.frame(dif_rec_sam_nonK9)
#TE_K9<-data.frame(rec_K9_TE$sam_dif); TE_nonK9<-data.frame(rec_nonK9_TE$sam_dif)
TE_K9$group<-"TE_K9";TE_nonK9$group<-"TE_non"
dat<-rbind(TE_K9, setNames(TE_nonK9, names(TE_K9)))

colnames(dat)[1] <- "log2_ratio_of_distance"
p<-ggplot(dat, aes(x=log2_ratio_of_distance, colour=group)) + geom_density(linewidth=2,bw=0.05)+ scale_color_manual(values = cols)+xlim(-0.15,0.8)##+xlim(-0.8,0.4)
#p<-ggplot(dat, aes(x=log2_ratio_of_distance,colour = group,fill = group)) + geom_histogram(alpha=0.5, position="identity")+ scale_color_manual(values = cols)+ scale_fill_manual(values = cols)#+xlim(-0.025,0.055) alpha=.2,binwidth=0.01
p<-p+mytheme_black; p
#ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/density_A6_log2_ratio_of_nearest_distance_betw_strains_control_TE_K9orNot_nearby_excluded.pdf", sep = ""), plot = p, width = 8, height = 5)
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/density_A7_log2_ratio_of_nearest_distance_betw_strains_control_TE_K9vsNonK9_alter_K9_excluded_black.pdf", sep = ""), plot = p, width = 8, height = 5)


cols <- c("grey","#CC79A7")
cols <- c("grey","#56B4E9")

p_MWU=as.data.frame(p_MWU_TE_K9)
p<-ggplot(p_MWU, aes(x=p_MWU_TE_K9)) + scale_x_continuous()+
  geom_histogram(color="black", fill="#56B4E9")+xlab("Mann-Whitney U test p-value")+ylab("counts")
p<-p+mytheme
p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/histogram_dif_log2_distance_CO_TE_K9_MWU_pvalue_A7.pdf", sep = ""), plot = p, width = 7, height = 6)

p_rho=as.data.frame(rho)
p<-ggplot(p_rho, aes(x=rho)) + scale_x_continuous()+
  geom_histogram(color="black", fill="#CC79A7",binwidth=0.01)+xlab("spearman correlation coefficient")+ylab("counts")
p<-p+mytheme
p
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/histogram_dif_log2_distance_CO_TE_length_rho_A6.pdf", sep = ""), plot = p, width = 7, height = 6)


