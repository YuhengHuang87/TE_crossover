I. Pipeline of identifying CO in PacBio pool-seq sequence data
1. map PacBio long reads to both parental genome assembles: map_pacbio_pool-seq.sub
2. Filter out reads with low mapping quality: perl select_reads_pass_mapping_quality.pl
3. Classify the parental origins of SNPs along each read. This step takes the longest, so consider splitting the sam file into small parts first: split -C 300m --numeric-suffixes aln_A6_benchmark_pool_pbmm2_sorted_que20_4kb.sam aln_A6_benchmark_pool_que20_4kb_prefix
   perl match_read_SNP_location.pl
   cat *_aln_A6_benchmark_pool_sorted_que20_read_posi > aln_A6_benchmark_pool_sorted_que20_read_posi
4. Identify reads that show switches of parental origins (process separately for reference A6/A7 and A4): perl sort_sliding_window_reads_chr_position.pl
5. Identify candidate CO reads and reads that show more than one parental origin switch: perl reads_share_refers_breakpoints.pl
6. Zoom in the location of CO and filter out reads when the CO breakpoint is located near the edge of the alignment: perl CO_location_exclu_double_CO_close_edge.pl



II. Compare the distribution of CO from our approach with an existing map from the standard genotyping approach
There are three files containing CO information: The map for the standard genotyping approach was obtained from Comeron et al. 2012: Comeron_2012_PLoSGent_inR6.txt; COs identified in A6 experimental pool: A6_cross_output.txt; and COs identified in A7 experimental pool: A7_cross_output.txt
1. Normalize CO numbers in each experimental pool by depth and convert coordinates to iso-1 release 6: perl the convert_v6_map_A6_A7_CO_depth_normalized.pl
2. combine outputs from different experimental pools: join -j 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.4,2.5,2.6,2.7,2.8 <(sort -k2 CO_broad_A6_focal_window_depth.txt) <(sort -k2 CO_broad_A7_focal_window_depth.txt) > comeron_pool_A6_A7_100kb_rate_depth_normalized
3. perform sliding windows analysis to show distributions of CO along chromosome arms: perl sliding_win_CO_cor_Comeron_A6_A7.pl


