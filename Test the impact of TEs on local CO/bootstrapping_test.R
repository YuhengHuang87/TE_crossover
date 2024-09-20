#######CO number
#### TE window pairs vs control window pairs
rec<-read.table("/dfs7/grylee/yuhenh3/recombination_rate/A6_crossover_num_TE_windows_left_right_window_magnitude_focal_homolog.txt",header=F)
foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V13),'/',fixed=TRUE)));rec$class <- foo$X1
#rec<-rec[rec$V1 != "roo", ]
rec<-rec[rec$V5 >= 5 & rec$V7 >= 5 & rec$V24 >=5 & rec$V26 >=5, ]
nearby_ambi<-rec[rec$V12 == "nearby_ambiguous",]
rec$K9_mag_left<-(rec$V14+rec$V15+rec$V16)/3;rec$K9_mag_right<-(rec$V17+rec$V18+rec$V19)/3;rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2
rec$K9_alt_left<-(rec$V31+rec$V32+rec$V33)/3;rec$K9_alt_right<-(rec$V34+rec$V35+rec$V36)/3;rec$alt_K9<-(rec$K9_alt_left+rec$K9_alt_right)/2
rec$chr<-rec$V1;rec$start<-rec$V2
rec$co<-rec$V8+rec$V10;rec$cov<-(rec$V4+rec$V6)/2;rec$alt_co<-rec$V27+rec$V29
rec$snp<-rec$V5+rec$V7;rec$TE <- paste(rec$V1, rec$V2, sep="_");rec$leng = rec$V3-rec$V2+1
#roo_en<-rec[rec$V12 =="roo" & rec$co > 0, ];rec<-rec[!(rec$TE %in% roo_en$TE),]
rec$alt_snp<-rec$V24+rec$V26;
rec$foc.cov=round((rec$V4+rec$V6)/4);rec$alt.cov=round((rec$V23+rec$V25)/4);
rec<-rec[rec$leng>=200,]
K9_cutoff=1
#rec<-rec[rec$K9 < K9_cutoff & rec$alt_K9 < K9_cutoff, ];

rec_as<-data.frame(rec$chr,rec$start, rec$co, rec$alt_co, rec$foc.cov, rec$alt.cov)

#rec_as<-data.frame(rec_K9_TE$V1,rec_K9_TE$V2, rec_K9_TE$co,rec_K9_TE$alt_co, rec_K9_TE$foc.cov, rec_K9_TE$alt.cov)
#rec_as<-data.frame(rec_nonK9_TE$V1,rec_nonK9_TE$V2, rec_nonK9_TE$co,rec_nonK9_TE$alt_co, rec_nonK9_TE$foc.cov, rec_nonK9_TE$alt.cov)

#control_1<-read.table("/dfs7/grylee/yuhenh3/recombination_rate/A7with_rec_local_control_window_location_include_nearby_ambig_10000_30kb",header=F)
#control_2<-read.table("/dfs7/grylee/yuhenh3/CUT_tag_experiment/control_focal_A7_alter_A6_TE_height_HMD_25_5000_50_mean_include_ambiguous.txt",header=F)
#control_1$K9=(control_1$V19+control_1$V20+control_1$V21+control_1$V22+control_1$V23+control_1$V24)/6 #switch to match with A6 or A7 indexs
#control_2$K9=(control_2$V19+control_2$V20+control_2$V21+control_2$V23+control_2$V24+control_2$V25)/6
#control_1<-control_1[control_1$K9 < 1, ];control_2<-control_2[control_2$K9 < 1, ]
#control<-control_1[control_1$V2 %in% intersect(control_1$V2, control_2$V2),];control$win <- paste(control$V1, control$V2, sep="_")
#cont<-control
#con_as<-data.frame(cont$V1,cont$V2,cont$co,cont$alt_co,cont$foc.cov,cont$alt.cov)
####or use input control file previously saved
con_as<-read.table("/dfs7/grylee/yuhenh3/recombination_rate/permutation/A6_control_CO_number_betw_strain_K9_SNP_ratio_filtered.txt",header=F)
colnames(con_as) <- c("cont.V1", "cont.V2","cont.co","cont.alt_co","cont.foc.cov","cont.alt.cov")

dif_CO_median=c()
rep=1
for (k in rep:(rep+1000-1)){print (k)
  chrX<-con_as[con_as$cont.V1=="CM010569.1",];focal_X<- chrX[sample(1:nrow(chrX), sum(rec_as$rec.chr=="CM010569.1"), replace=T),];
  chr2L<-con_as[con_as$cont.V1=="CM010570.1",];focal_2L<- chr2L[sample(1:nrow(chr2L), sum(rec_as$rec.chr=="CM010570.1"), replace=T),];
  chr2R<-con_as[con_as$cont.V1=="CM010571.1",];focal_2R<- chr2R[sample(1:nrow(chr2R), sum(rec_as$rec.chr=="CM010571.1"), replace=T),];
  chr3L<-con_as[con_as$cont.V1=="CM010572.1",];focal_3L<- chr3L[sample(1:nrow(chr3L), sum(rec_as$rec.chr=="CM010572.1"), replace=T),];
  chr3R<-con_as[con_as$cont.V1=="CM010573.1",];focal_3R<- chr3R[sample(1:nrow(chr3R), sum(rec_as$rec.chr=="CM010573.1"), replace=T),];
  
  #chrX<-con_as[con_as$cont.V1=="CM010576.1",];focal_X<- chrX[sample(1:nrow(chrX), sum(rec_as$rec.chr=="CM010576.1"), replace=T),];
  #chr2L<-con_as[con_as$cont.V1=="CM010577.1",];focal_2L<- chr2L[sample(1:nrow(chr2L), sum(rec_as$rec.chr=="CM010577.1"), replace=T),];
  #chr2R<-con_as[con_as$cont.V1=="CM010578.1",];focal_2R<- chr2R[sample(1:nrow(chr2R), sum(rec_as$rec.chr=="CM010578.1"), replace=T),];
  #chr3L<-con_as[con_as$cont.V1=="CM010579.1",];focal_3L<- chr3L[sample(1:nrow(chr3L), sum(rec_as$rec.chr=="CM010579.1"), replace=T),];
  #chr3R<-con_as[con_as$cont.V1=="CM010580.1",];focal_3R<- chr3R[sample(1:nrow(chr3R), sum(rec_as$rec.chr=="CM010580.1"), replace=T),];
  
  per_focal<-rbind(focal_X,focal_2L,focal_2R,focal_3L,focal_3R)
  
  dif_rec_sam=c();dif_rec_sam_K9=c();dif_rec_sam_nonK9=c();
  dif_cont_sam=c()
  
  for (j in 1:100){
    for (i in 1:nrow(per_focal)){
      if (per_focal$cont.foc.cov[i] > per_focal$cont.alt.cov[i]){
        if (per_focal$cont.co[i] > 0){
          s_list<-sample.int(per_focal$cont.foc.cov[i], size = per_focal$cont.alt.cov[i], replace = FALSE)
          per_focal$focal_co[i] = length(s_list[s_list<=per_focal$cont.co[i]])
          per_focal$focal_cov[i] = per_focal$cont.alt.cov[i]
          per_focal$alter_co[i] = per_focal$cont.alt_co[i]
          per_focal$alter_cov[i] = per_focal$cont.alt.cov[i]
        }else{
          per_focal$focal_co[i] = per_focal$cont.co[i]
          per_focal$focal_cov[i] = per_focal$cont.alt.cov[i]
          per_focal$alter_co[i] = per_focal$cont.alt_co[i]
          per_focal$alter_cov[i] = per_focal$cont.alt.cov[i]
        }
      }else{
        if (per_focal$cont.alt_co[i] > 0){
          s_list<-sample.int(per_focal$cont.alt.cov[i], size =per_focal$cont.foc.cov[i] , replace = FALSE)
          per_focal$focal_co[i] = per_focal$cont.co[i]
          per_focal$focal_cov[i] = per_focal$cont.foc.cov[i]
          per_focal$alter_co[i] = length(s_list[s_list<=per_focal$cont.alt_co[i]])
          per_focal$alter_cov[i] = per_focal$cont.foc.cov[i]
        }else{
          per_focal$focal_co[i] = per_focal$cont.co[i]
          per_focal$focal_cov[i] = per_focal$cont.foc.cov[i]
          per_focal$alter_co[i] = per_focal$cont.alt_co[i]
          per_focal$alter_cov[i] = per_focal$cont.foc.cov[i]
        }
      }}
    per_focal$sam_dif = per_focal$focal_co-per_focal$alter_co;
    dif_cont_sam=c(dif_cont_sam, mean(per_focal$sam_dif));
  }
  
  dif_CO_median=c(dif_CO_median,median(dif_cont_sam))
}
outputname<- paste0("/dfs7/grylee/yuhenh3/recombination_rate/permutation/bootstrap_control_nonK9_TE_homolog_filter_SNP_ratio_dif_CO_num_A6_",rep, ".txt")
write.table(dif_CO_median, file=outputname,sep="\t",row.names = F,col.names = F)

##### comparing TE window pairs with and without K9
rec<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7_crossover_num_TE_windows_left_right_window_magnitude_focal_homolog.txt",header=F)
foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V13),'/',fixed=TRUE)));rec$class <- foo$X1
#rec<-rec[rec$V1 != "roo", ]
rec<-rec[rec$V5 >= 5 & rec$V7 >= 5 & rec$V24 >=5 & rec$V26 >=5, ]

rec$K9_mag_left<-(rec$V14+rec$V15+rec$V16)/3;rec$K9_mag_right<-(rec$V17+rec$V18+rec$V19)/3;rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2
rec$K9_alt_left<-(rec$V31+rec$V32+rec$V33)/3;rec$K9_alt_right<-(rec$V34+rec$V35+rec$V36)/3;rec$alt_K9<-(rec$K9_alt_left+rec$K9_alt_right)/2
rec$co<-rec$V8+rec$V10;rec$cov<-(rec$V4+rec$V6)/2;rec$alt_co<-rec$V27+rec$V29
rec$snp<-rec$V5+rec$V7;rec$TE <- paste(rec$V1, rec$V2, sep="_");rec$leng = rec$V3-rec$V2+1
rec$alt_snp<-rec$V24+rec$V26;
rec$foc.cov=round((rec$V4+rec$V6)/2);rec$alt.cov=round((rec$V23+rec$V25)/2);
#rec_unenriched<-rec[rec$K9<1,]
rec<-rec[rec$leng>=200,]
rec<-rec[rec$alt_K9<1,]
K9_cutoff=1
rec_K9_TE<-rec[rec$K9 > K9_cutoff, ];rec_nonK9_TE<-rec[rec$K9 < K9_cutoff, ];


rec_as<-data.frame(rec_K9_TE$V1,rec_K9_TE$V2, rec_K9_TE$co,rec_K9_TE$alt_co, rec_K9_TE$foc.cov, rec_K9_TE$alt.cov)

con_as<-data.frame(rec_nonK9_TE$V1,rec_nonK9_TE$V2,rec_nonK9_TE$co,rec_nonK9_TE$alt_co,rec_nonK9_TE$foc.cov,rec_nonK9_TE$alt.cov)
colnames(rec_as) <- c("rec.V1", "rec.V2","rec.co","rec.alt_co","rec.foc.cov","rec.alt.cov")
colnames(con_as) <- c("rec.V1", "rec.V2","rec.co","rec.alt_co","rec.foc.cov","rec.alt.cov")

dif_CO_median=c()
rep=1
for (k in rep:(rep+1000-1)){print (k)
  #chrX<-con_as[con_as$rec.V1=="CM010569.1",];focal_X<- chrX[sample(1:nrow(chrX), sum(rec_as$rec.V1=="CM010569.1"), replace=T),];
  #chr2L<-con_as[con_as$rec.V1=="CM010570.1",];focal_2L<- chr2L[sample(1:nrow(chr2L), sum(rec_as$rec.V1=="CM010570.1"), replace=T),];
  #chr2R<-con_as[con_as$rec.V1=="CM010571.1",];focal_2R<- chr2R[sample(1:nrow(chr2R), sum(rec_as$rec.V1=="CM010571.1"), replace=T),];
  #chr3L<-con_as[con_as$rec.V1=="CM010572.1",];focal_3L<- chr3L[sample(1:nrow(chr3L), sum(rec_as$rec.V1=="CM010572.1"), replace=T),];
  #chr3R<-con_as[con_as$rec.V1=="CM010573.1",];focal_3R<- chr3R[sample(1:nrow(chr3R), sum(rec_as$rec.V1=="CM010573.1"), replace=T),];
  
  chrX<-con_as[con_as$rec.V1=="CM010576.1",];focal_X<- chrX[sample(1:nrow(chrX), sum(rec_as$rec.V1=="CM010576.1"), replace=T),];
  chr2L<-con_as[con_as$rec.V1=="CM010577.1",];focal_2L<- chr2L[sample(1:nrow(chr2L), sum(rec_as$rec.V1=="CM010577.1"), replace=T),];
  chr2R<-con_as[con_as$rec.V1=="CM010578.1",];focal_2R<- chr2R[sample(1:nrow(chr2R), sum(rec_as$rec.V1=="CM010578.1"), replace=T),];
  chr3L<-con_as[con_as$rec.V1=="CM010579.1",];focal_3L<- chr3L[sample(1:nrow(chr3L), sum(rec_as$rec.V1=="CM010579.1"), replace=T),];
  chr3R<-con_as[con_as$rec.V1=="CM010580.1",];focal_3R<- chr3R[sample(1:nrow(chr3R), sum(rec_as$rec.V1=="CM010580.1"), replace=T),];
  
  per_focal<-rbind(focal_X,focal_2L,focal_2R,focal_3L,focal_3R)
  
  dif_rec_sam=c();dif_rec_sam_K9=c();dif_rec_sam_nonK9=c();
  dif_cont_sam=c()
  
  for (j in 1:100){
    for (i in 1:nrow(per_focal)){
      if (per_focal$rec.foc.cov[i] > per_focal$rec.alt.cov[i]){
        if (per_focal$rec.co[i] > 0){
          s_list<-sample.int(per_focal$rec.foc.cov[i], size = per_focal$rec.alt.cov[i], replace = FALSE)
          per_focal$focal_co[i] = length(s_list[s_list<=per_focal$rec.co[i]])
          per_focal$focal_cov[i] = per_focal$rec.alt.cov[i]
          per_focal$alter_co[i] = per_focal$rec.alt_co[i]
          per_focal$alter_cov[i] = per_focal$rec.alt.cov[i]
        }else{
          per_focal$focal_co[i] = per_focal$rec.co[i]
          per_focal$focal_cov[i] = per_focal$rec.alt.cov[i]
          per_focal$alter_co[i] = per_focal$rec.alt_co[i]
          per_focal$alter_cov[i] = per_focal$rec.alt.cov[i]
        }
      }else{
        if (per_focal$rec.alt_co[i] > 0){
          s_list<-sample.int(per_focal$rec.alt.cov[i], size =per_focal$rec.foc.cov[i] , replace = FALSE)
          per_focal$focal_co[i] = per_focal$rec.co[i]
          per_focal$focal_cov[i] = per_focal$rec.foc.cov[i]
          per_focal$alter_co[i] = length(s_list[s_list<=per_focal$rec.alt_co[i]])
          per_focal$alter_cov[i] = per_focal$rec.foc.cov[i]
        }else{
          per_focal$focal_co[i] = per_focal$rec.co[i]
          per_focal$focal_cov[i] = per_focal$rec.foc.cov[i]
          per_focal$alter_co[i] = per_focal$rec.alt_co[i]
          per_focal$alter_cov[i] = per_focal$rec.foc.cov[i]
        }
      }}
    per_focal$sam_dif = per_focal$focal_co-per_focal$alter_co;
    dif_rec_sam=c(dif_rec_sam, mean(per_focal$sam_dif));
    
  }
  
  dif_CO_median=c(dif_CO_median,median(dif_rec_sam))
}
outputname<- paste0("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/bootstrap_dif_CO_num_TE_K9_control_A7_",rep, ".txt")
write.table(dif_CO_median, file=outputname,sep="\t",row.names = F,col.names = F)
#########

######### Distance to the nearest CO
#### compare between TE window pairs and control pairs
rec<-read.table("/dfs7/grylee/yuhenh3/recombination_rate/A6_crossover_nearest_distance_coverage_SNP_TE_left_right_window_mass_5000_exclude_distance_5000",header=F)
foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V5),'/',fixed=TRUE)));rec$class <- foo$X1
rec<-rec[rec$class == "DNA", ]
rec<-rec[rec$class != "Unknown", ];rec$TE_type<-ifelse(rec$class == "DNA","DNA","RNA")

rec$leng = rec$V4-rec$V3+1;
rec$snp_den=rec$V9/rec$V6;rec$con_snp_den=rec$V22/rec$V21;
rec$K9_mag_left<-abs(rec$V10+rec$V11+rec$V12)/3
rec$K9_mag_right<-abs(rec$V13+rec$V14+rec$V15)/3
rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2
rec<-rec[rec$leng>=200,]

rec$K9_alter_left<-abs(rec$V23+rec$V24+rec$V25)/3
rec$K9_alter_right<-abs(rec$V26+rec$V27+rec$V28)/3
rec$K9_alter<-(rec$K9_alter_left+rec$K9_alter_right)/2
rec<-rec[rec$K9_alter<3,]

rec$dif_K9=log2(rec$K9/rec$K9_alter)
rec$foc.cov<-round(rec$V8/2);rec$alt.cov<-round(rec$V21/2)
K9_cutoff=median(rec$K9)

#rec<-rec[rec$K9 > K9_cutoff, ];
rec<-rec[rec$K9 < K9_cutoff, ]

control_1<-read.table("/dfs7/grylee/yuhenh3/CUT_tag_experiment/control_focal_K9_magnitude_A6_alter_A7_nearest_CO_distance_mean_Mass_TE_HMD_25_1000_10_measures_include_ambiguous.txt",header=F)
control_2<-read.table("/dfs7/grylee/yuhenh3/CUT_tag_experiment/control_focal_A6_alter_K9_magnitude_A7_nearest_CO_distance_TE_HMD_25_1000_10_measures_include_ambiguous.txt",header=F)

control_1$K9=(control_1$V15+control_1$V16+control_1$V17+control_1$V18+control_1$V19+control_1$V20)/6
control_2$K9=(control_2$V15+control_2$V16+control_2$V17+control_2$V18+control_2$V19+control_2$V20)/6

control_1<-control_1[control_1$K9 < 1, ];control_2<-control_2[control_2$K9 < 1, ]
control<-control_1[control_1$V2 %in% intersect(control_1$V2, control_2$V2),];control$win <- paste(control$V1, control$V2, sep="_")
cont<-control


cont$foc.cov = round(cont$V6); cont$alt.cov = round(cont$V13);cont$dif_cov = cont$foc.cov-cont$alt.cov;cont$ratio_cov = cont$foc.cov/cont$alt.cov
cont$snp_den = cont$V7/cont$V6;cont$alt_snp_den = cont$V14/cont$V13;cont$dif_snp_den=cont$snp_den-cont$alt_snp_den;

rec_as<-data.frame(rec$V2, rec$V6, rec$V7, rec$V19, rec$V20, rec$foc.cov, rec$alt.cov)
#con_as<-data.frame(cont$V1,cont$V4,cont$V5,cont$V11,cont$V12,cont$foc.cov,cont$alt.cov)
####or use input control file previously saved
con_as<-read.table("/dfs7/grylee/yuhenh3/recombination_rate/permutation/A6_control_nearest_distance_betw_strain_K9_SNP_ratio_filtered.txt",header=T)
####

dif_dis_median=c()
rep=1
for (k in rep:(rep+1000-1)){print (k)
  chrX<-con_as[con_as$cont.V1=="CM010569.1",];focal_X<- chrX[sample(1:nrow(chrX), sum(rec_as$rec.V2=="CM010569.1"), replace=T),];
  chr2L<-con_as[con_as$cont.V1=="CM010570.1",];focal_2L<- chr2L[sample(1:nrow(chr2L), sum(rec_as$rec.V2=="CM010570.1"), replace=T),];
  chr2R<-con_as[con_as$cont.V1=="CM010571.1",];focal_2R<- chr2R[sample(1:nrow(chr2R), sum(rec_as$rec.V2=="CM010571.1"), replace=T),];
  chr3L<-con_as[con_as$cont.V1=="CM010572.1",];focal_3L<- chr3L[sample(1:nrow(chr3L), sum(rec_as$rec.V2=="CM010572.1"), replace=T),];
  chr3R<-con_as[con_as$cont.V1=="CM010573.1",];focal_3R<- chr3R[sample(1:nrow(chr3R), sum(rec_as$rec.V2=="CM010573.1"), replace=T),];
  
  #chrX<-con_as[con_as$cont.V1=="CM010576.1",];focal_X<- chrX[sample(1:nrow(chrX), sum(rec_as$rec.V2=="CM010576.1"), replace=T),];
  #chr2L<-con_as[con_as$cont.V1=="CM010577.1",];focal_2L<- chr2L[sample(1:nrow(chr2L), sum(rec_as$rec.V2=="CM010577.1"), replace=T),];
  #chr2R<-con_as[con_as$cont.V1=="CM010578.1",];focal_2R<- chr2R[sample(1:nrow(chr2R), sum(rec_as$rec.V2=="CM010578.1"), replace=T),];
  #chr3L<-con_as[con_as$cont.V1=="CM010579.1",];focal_3L<- chr3L[sample(1:nrow(chr3L), sum(rec_as$rec.V2=="CM010579.1"), replace=T),];
  #chr3R<-con_as[con_as$cont.V1=="CM010580.1",];focal_3R<- chr3R[sample(1:nrow(chr3R), sum(rec_as$rec.V2=="CM010580.1"), replace=T),];
  
  per_focal<-rbind(focal_X,focal_2L,focal_2R,focal_3L,focal_3R)
  
  dif_rec_sam=c();dif_rec_sam_K9=c();dif_rec_sam_nonK9=c();
  dif_cont_sam=c()
  
  for (j in 1:100){
    for (i in 1:nrow(per_focal)){
      if (per_focal$cont.foc.cov[i] > per_focal$cont.alt.cov[i]){
        s_list<-sample.int(per_focal$cont.foc.cov[i], size = per_focal$cont.alt.cov[i], replace = FALSE)
        if (1%in% s_list){
          per_focal$foc_dis[i] = per_focal$cont.V4[i]
          per_focal$alt_dis[i] = per_focal$cont.V11[i]
        }else{
          per_focal$foc_dis[i] = per_focal$cont.V5[i]
          per_focal$alt_dis[i] = per_focal$cont.V11[i]
        }
      }else{
        s_list<-sample.int(per_focal$cont.alt.cov[i], size =per_focal$cont.foc.cov[i] , replace = FALSE)
        if (1%in% s_list){
          per_focal$foc_dis[i] = per_focal$cont.V4[i]
          per_focal$alt_dis[i] = per_focal$cont.V11[i]
        }else{
          per_focal$foc_dis[i] = per_focal$cont.V4[i]
          per_focal$alt_dis[i] = per_focal$cont.V12[i]
        }
      }}
    #rec$sam_dif = rec$foc_dis-rec$alt_dis;
    per_focal$sam_dif = log2(per_focal$foc_dis/per_focal$alt_dis);
    dif_rec_sam=c(dif_rec_sam, median(per_focal$sam_dif));
  }
  
  dif_dis_median=c(dif_dis_median,median(dif_rec_sam))
}
outputname<- paste0("/dfs7/grylee/yuhenh3/recombination_rate/permutation/bootstrap_log2_control_TE_nonK9_homolog_extra_SNPratio_A6_",rep, ".txt")
write.table(dif_dis_median, file=outputname,sep="\t",row.names = F,col.names = F)

#### compare window pairs with high and low K9 mass
rec<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/A7_crossover_nearest_distance_coverage_SNP_TE_left_right_window_mass_5000_exclude_distance_5000",header=F)

foo <- data.frame(do.call('rbind', strsplit(as.character(rec$V5),'/',fixed=TRUE)));rec$class <- foo$X1
#rec<-rec[rec$V1 != "roo", ]
rec$leng = rec$V4-rec$V3+1;
rec$snp_den=rec$V9/rec$V6;rec$con_snp_den=rec$V22/rec$V21;
rec$K9_mag_left<-abs(rec$V10+rec$V11+rec$V12)/3
rec$K9_mag_right<-abs(rec$V13+rec$V14+rec$V15)/3
rec$K9<-(rec$K9_mag_left+rec$K9_mag_right)/2
rec<-rec[rec$leng>=200,]
#rec_as<-rec_as[rec_as$K9 > 0, ]

rec$K9_alter_left<-abs(rec$V23+rec$V24+rec$V25)/3
rec$K9_alter_right<-abs(rec$V26+rec$V27+rec$V28)/3
rec$K9_alter<-(rec$K9_alter_left+rec$K9_alter_right)/2
#median(rec$K9)
rec<-rec[rec$K9_alter<3,]
rec$dif_K9=log2(rec$K9/rec$K9_alter)
high_mass<-rec[rec$K9>50,]
rec$foc.cov<-round(rec$V8/2);rec$alt.cov<-round(rec$V21/2)
K9_cutoff=median(rec$K9)
rec_K9_TE<-rec[rec$K9 > K9_cutoff, ];
rec_nonK9_TE<-rec[rec$K9 < K9_cutoff, ]

rec_as<-data.frame(rec_K9_TE$V2, rec_K9_TE$V6, rec_K9_TE$V7, rec_K9_TE$V19, rec_K9_TE$V20, rec_K9_TE$foc.cov, rec_K9_TE$alt.cov)
con_as<-data.frame(rec_nonK9_TE$V2, rec_nonK9_TE$V6, rec_nonK9_TE$V7, rec_nonK9_TE$V19, rec_nonK9_TE$V20, rec_nonK9_TE$foc.cov, rec_nonK9_TE$alt.cov)

dif_dis_median=c()
for (k in 1:1000){print (k)
  
  #chrX<-con_as[con_as$rec_nonK9_TE.V2=="CM010569.1",];focal_X<- chrX[sample(1:nrow(chrX), sum(rec_as$rec_K9_TE.V2=="CM010569.1"), replace=T),];
  #chr2L<-con_as[con_as$rec_nonK9_TE.V2=="CM010570.1",];focal_2L<- chr2L[sample(1:nrow(chr2L), sum(rec_as$rec_K9_TE.V2=="CM010570.1"), replace=T),];
  #chr2R<-con_as[con_as$rec_nonK9_TE.V2=="CM010571.1",];focal_2R<- chr2R[sample(1:nrow(chr2R), sum(rec_as$rec_K9_TE.V2=="CM010571.1"), replace=T),];
  #chr3L<-con_as[con_as$rec_nonK9_TE.V2=="CM010572.1",];focal_3L<- chr3L[sample(1:nrow(chr3L), sum(rec_as$rec_K9_TE.V2=="CM010572.1"), replace=T),];
  #chr3R<-con_as[con_as$rec_nonK9_TE.V2=="CM010573.1",];focal_3R<- chr3R[sample(1:nrow(chr3R), sum(rec_as$rec_K9_TE.V2=="CM010573.1"), replace=T),];
  
  chrX<-con_as[con_as$rec_nonK9_TE.V2=="CM010576.1",];focal_X<- chrX[sample(1:nrow(chrX), sum(rec_as$rec_K9_TE.V2=="CM010576.1"), replace=T),];
  chr2L<-con_as[con_as$rec_nonK9_TE.V2=="CM010577.1",];focal_2L<- chr2L[sample(1:nrow(chr2L), sum(rec_as$rec_K9_TE.V2=="CM010577.1"), replace=T),];
  chr2R<-con_as[con_as$rec_nonK9_TE.V2=="CM010578.1",];focal_2R<- chr2R[sample(1:nrow(chr2R), sum(rec_as$rec_K9_TE.V2=="CM010578.1"), replace=T),];
  chr3L<-con_as[con_as$rec_nonK9_TE.V2=="CM010579.1",];focal_3L<- chr3L[sample(1:nrow(chr3L), sum(rec_as$rec_K9_TE.V2=="CM010579.1"), replace=T),];
  chr3R<-con_as[con_as$rec_nonK9_TE.V2=="CM010580.1",];focal_3R<- chr3R[sample(1:nrow(chr3R), sum(rec_as$rec_K9_TE.V2=="CM010580.1"), replace=T),];
  
  per_focal<-rbind(focal_X,focal_2L,focal_2R,focal_3L,focal_3R)
  
  dif_rec_sam=c();dif_rec_sam_K9=c();dif_rec_sam_nonK9=c();
  dif_cont_sam=c()
  
  for (j in 1:100){
    for (i in 1:nrow(per_focal)){
      if (per_focal$rec_nonK9_TE.foc.cov[i] > per_focal$rec_nonK9_TE.alt.cov[i]){
        s_list<-sample.int(per_focal$rec_nonK9_TE.foc.cov[i], size = per_focal$rec_nonK9_TE.alt.cov[i], replace = FALSE)
        if (1%in% s_list){
          per_focal$foc_dis[i] = per_focal$rec_nonK9_TE.V6[i]
          per_focal$alt_dis[i] = per_focal$rec_nonK9_TE.V19[i]
        }else{
          per_focal$foc_dis[i] = per_focal$rec_nonK9_TE.V7[i]
          per_focal$alt_dis[i] = per_focal$rec_nonK9_TE.V19[i]
        }
      }else{
        s_list<-sample.int(per_focal$rec_nonK9_TE.alt.cov[i], size =per_focal$rec_nonK9_TE.foc.cov[i] , replace = FALSE)
        if (1%in% s_list){
          per_focal$foc_dis[i] = per_focal$rec_nonK9_TE.V6[i]
          per_focal$alt_dis[i] = per_focal$rec_nonK9_TE.V19[i]
        }else{
          per_focal$foc_dis[i] = per_focal$rec_nonK9_TE.V6[i]
          per_focal$alt_dis[i] = per_focal$rec_nonK9_TE.V20[i]
        }
      }}
    #rec$sam_dif = rec$foc_dis-rec$alt_dis;
    per_focal$sam_dif = log2(per_focal$foc_dis/per_focal$alt_dis);
    dif_rec_sam=c(dif_rec_sam, median(per_focal$sam_dif));
  }
  dif_dis_median=c(dif_dis_median,median(dif_rec_sam))
}
#write.table(dif_dis_median_dif, file="/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/dif_dis_TE_types_A7.txt",sep="\t",row.names = T,col.names = T)
write.table(dif_dis_median, file="/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/bootstrap_log2_TE_K9vsNonK9_nearby_excluded_5000_A7.txt",sep="\t",row.names = F,col.names = F)
