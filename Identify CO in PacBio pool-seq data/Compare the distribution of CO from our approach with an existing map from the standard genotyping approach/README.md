
II. Compare the distribution of CO from our approach with an existing map from the standard genotyping approach
There are three files containing CO information: The map for the standard genotyping approach was obtained from Comeron et al. 2012: Comeron_2012_PLoSGent_inR6.txt; COs identified in A6 experimental pool: A6_cross_output.txt; and COs identified in A7 experimental pool: A7_cross_output.txt
1. Estimate depth in benchmark, A6, and A7 experimental pools after read filtering: perl depth_per_SNP_pass_screening_combine.pl
2. Normalize CO numbers in each experimental pool by depth and convert coordinates to iso-1 release 6: perl the convert_v6_map_A6_A7_CO_depth_normalized.pl
3. combine outputs from different experimental pools: join -j 2 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.4,2.5,2.6,2.7,2.8 <(sort -k2 CO_broad_A6_focal_window_depth.txt) <(sort -k2 CO_broad_A7_focal_window_depth.txt) > comeron_pool_A6_A7_100kb_rate_depth_normalized
4. perform sliding windows analysis to show distributions of CO along chromosome arms: perl sliding_win_CO_cor_Comeron_A6_A7.pl
5. analysis and plot the data: rec_map_comeron_A6_A7.R


