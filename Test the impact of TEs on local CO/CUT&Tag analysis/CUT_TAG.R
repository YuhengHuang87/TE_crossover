library(dplyr);library(ggplot2);library(corrplot)

sampleList = c("A6_1", "A6_2", "A6_4")
sampleList = c("A7_2","A7_3","A7_4")
histList = c("A6","A7")

## Collect the alignment results from the bowtie2 alignment summary files; check correlation in H3K9me3 across replicates
alignResult = c()
for(hist in sampleList){
  alignRes = read.table(paste0("/Users/yuhenghuang/Documents/Unix_and_Perl/Code/Recombination/CUT_tag/output_", hist, "_bowtie2.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo = strsplit(hist, "_")[[1]]
  alignResult = data.frame(Histone = histInfo[1], Replicate = histInfo[2], 
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                           MappedFragNum_hg38 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric, 
                           AlignmentRate_hg38 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}
alignResult$Histone = factor(alignResult$Histone, levels = histList)
alignResult %>% mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"))

reprod = c()
fragCount = NULL
for(hist in sampleList){
  if(is.null(fragCount)){
    fragCount = read.table(paste0("/Users/yuhenghuang/Documents/Unix_and_Perl/Code/Recombination/CUT_tag/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
  }else{
    fragCountTmp = read.table(paste0("/Users/yuhenghuang/Documents/Unix_and_Perl/Code/Recombination/CUT_tag/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
  }
}

#cor.test(log2(fragCount$A7_2),log2(fragCount$A7_3),method = "pearson")
#cor.test(log2(fragCount$A6_2),log2(fragCount$A6_1),method = "spearman")

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 

corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))


p<-read.table("/Users/yuhenghuang/Documents/Unix_and_Perl/Code/Recombination/CUT_tag/A7_average_TE_height_HMD_25_1kb_10measures.txt",header=F)#
cor.test(p$V7,p$V8,method="spearman")
plot(p$V6,p$V7,xlim=c(0,10),ylim=c(0,10))
nrow(p);
sum(p$V6>1);sum(p$V7>1);

######## plot average H3K9me collate between strains Fig. S7
enrich_1<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/K9_Collate_betw_strain_A6_1_TE_200bp_50kb_win_1000_10000_local_K9_standardized_mean_include_ambiguous.txt",header=F)#A7 #2, 3, 5
enrich_2<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/K9_Collate_betw_strain_A6_2_TE_200bp_50kb_win_1000_10000_local_K9_standardized_mean_include_ambiguous.txt",header=F)#
enrich_3<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/K9_Collate_betw_strain_A6_4_TE_200bp_50kb_win_1000_10000_local_K9_standardized_mean_include_ambiguous.txt",header=F)#
enrich_p<-(enrich_1+enrich_2+enrich_3)/3

enrich_1<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/K9_Collate_betw_strain_A6_alternative_A7_2_TE_200bp_50kb_win_1000_10000_local_K9_standardized_mean_include_ambiguous.txt",header=F)#
enrich_2<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/K9_Collate_betw_strain_A6_alternative_A7_3_TE_200bp_50kb_win_1000_10000_local_K9_standardized_mean_include_ambiguous.txt",header=F)#
enrich_3<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/K9_Collate_betw_strain_A6_alternative_A7_4_TE_200bp_50kb_win_1000_10000_local_K9_standardized_mean_include_ambiguous.txt",header=F)#
enrich_a<-(enrich_1+enrich_2+enrich_3)/3

enrich_p$group <- "A6"
enrich_a$group <- "A7"
enrich <- rbind(enrich_p, enrich_a)
enrich<-enrich[(abs(enrich$V1)<21),]
#p <- ggplot(enrich, aes(x=V1, y=V2, group=group, col=factor(group), fill=factor(group)))+scale_color_manual(values=c("orange", "blue")) +scale_fill_manual(values=c("orange", "blue")) + geom_smooth(size=1,method = "loess",span = 0.15) + theme(axis.line.x = element_line(color="black", size = 0.7), axis.line.y = element_line(color="black", size = 0.7),axis.title.x = element_text(face="bold",  size=22),axis.title.y = element_text(angle=90,size=16,face="bold"), axis.text.x  = element_text(vjust=1, size=18,colour="Black"),axis.text.y = element_text(hjust=1, size=18,colour="Black"), legend.text = element_text(hjust=20,size=10,face="italic"),legend.title = element_text(hjust=0,size=0,face="bold"),panel.background = element_rect(fill = "white"))+ ylab("HMD")+xlab("distance to TE (kb)")
p <- ggplot(enrich, aes(x=V1, y=V2, group=group, col=factor(group), fill=factor(group)))+scale_color_manual(values=c("#CC79A7", "#56B4E9")) +scale_fill_manual(values=c("#CC79A7", "#56B4E9")) + geom_smooth(size=3,method = "loess",span = 0.2, se = F) + theme(axis.line.x = element_line(color="white", size = 0.7), axis.line.y = element_line(color="white", size = 0.7),axis.title.x = element_text(face="bold",  size=22),axis.title.y = element_text(angle=90,size=16,face="bold"), axis.text.x  = element_text(vjust=1, size=18,color="white"),axis.text.y = element_text(hjust=1, size=18,color="white"), legend.text = element_text(hjust=20,size=10,face="italic"),legend.title = element_text(hjust=0,size=0,face="bold"))+ ylab("K9 level")+xlab("distance to TE (kb)")
#p <- ggplot(enrich, aes(x=V1, y=V2, group=group, col=factor(group)))+scale_color_manual(values=c("blue", "black")) + geom_smooth(size=1,method = "loess",span = 0.15, se=F) + theme(axis.line.x = element_line(color="black", size = 0.7), axis.line.y = element_line(color="black", size = 0.7),axis.title.x = element_text(face="bold",  size=22),axis.title.y = element_text(angle=90,size=16,face="bold"), axis.text.x  = element_text(vjust=1, size=18,colour="Black"),axis.text.y = element_text(hjust=1, size=18,colour="Black"), legend.text = element_text(hjust=20,size=10,face="italic"),legend.title = element_text(hjust=0,size=0,face="bold"),panel.background = element_rect(fill = "white"))+ ylab("HMD")+xlab("distance to TE (kb)")
#geom_point(size=0) +
p1<-p+mytheme;p1
ggsave(file = paste("/Users/yuhenghuang/Documents/Postdoc_UCI/recombination/result/plot/K9_collate_A6_focal_A7_alternative_20_40kb_local_normalized.pdf", sep = ""), plot = p1, width = 6, height = 4)


