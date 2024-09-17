I. Pipeline of identifying CO in PacBio pool-seq sequence data
1. map PacBio long reads to both parental genome assembles: map_pacbio_pool-seq.sub
2. Filter out reads with low mapping quality: perl select_reads_pass_mapping_quality.pl
3. Classify the parental origins of SNPs along each read. This step takes the longest, so consider splitting the sam file into small parts first: split -C 300m --numeric-suffixes aln_A6_benchmark_pool_pbmm2_sorted_que20_4kb.sam aln_A6_benchmark_pool_que20_4kb_prefix &
   perl match_read_SNP_location.pl
   cat *_aln_A6_benchmark_pool_sorted_que20_read_posi > aln_A6_benchmark_pool_sorted_que20_read_posi
5. Identify reads show a switch of parental origins



II. Compare the distribution of CO from our approach with an existing map from the standard genotyping approach


