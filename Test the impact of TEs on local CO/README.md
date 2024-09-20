I. Annotate TEs to family level with scripts in folder /assign TE family

II. Analyze the CUT&Tag sequence data with scripts in folder /CUT&Tag analysis

III. Within strains analysis

a. CO number
1. Obtain the SNP number and depth for the flanking windows around TEs and control windows: perl assign_SNP_num_depth_window_flanking_TE_control_within_strains.pl
2. Obtain the CO number for the TE flanking windows and control window: perl assign_CO_TE_flanking_control_within_strains.pl


b. Distance to the nearest CO
1. Calculate the distance to the 1st and 2nd nearest COs for TEs and control location: perl distance_nearest_CO_TE_flanking_control.pl
2. Estimate the SNP density and depth for the region for TEs and control location to the 1st and 2nd nearest COs: perl assign_SNP_depth_nearest_CO_distance_flanking_TE_control_within_strains.pl

C. Analyze and plot the data: Within_strain_analysis.R


IV. Between strains analysis
a. Identify focal and homologous alleles
1. Identify homologous alleles in alternative strains: perl TE_control_window_alignment_convert.pl
2. Exclude homologous alleles with TEs located nearby: perl focal_homolog_allele_exclude_nearby_TE.pl

b. CO number
1. Obtain the SNP number and depth for the homologous allele without TEs and control alleles in both strains: perl assign_SNP_num_depth_window_flanking_TE_control_between_strains.pl
2. calculate the CO numbers for flanking regions of alleles without TEs and control alleles in both strains: perl assign_CO_TE_flanking_control_between_strains.pl

c. Distance to the nearest CO
1. Use the same script to calculate the distance to the 1st and 2nd nearest COs in homologous alleles: perl distance_nearest_CO_TE_flanking_control.pl
2. Estimate the SNP density and depth for regions for these alleles to the 1st and 2nd nearest COs: perl assign_SNP_depth_nearest_CO_distance_flanking_TE_control_between_strains.pl

d. Analyze and plot the data
1. Make figures and perform Mann–Whitney U tests: Between_strains_analysis.R
2. Perform bootstrapping tests: bootstrapping_test.R



