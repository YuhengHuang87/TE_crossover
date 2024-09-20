I. run scripts in folder /assign TE family

II. Analyze the CUT&Tag sequence data with scripts in folder /CUT&Tag analysis

III. Within strains analysis

a. CO number
1. Obtain the SNP number and depth for the flanking windows around TEs and control windows: perl assign_SNP_num_depth_window_flanking_TE_control_within_strains.pl
2. Obtain the CO number for the TE flanking windows and control window: perl assign_CO_TE_flanking_control_within_strains.pl


b. Distance to the nearest CO
1. Calculate the distance to the 1st and 2nd nearest COs for TEs and control location: perl distance_nearest_CO_TE_flanking_control.pl
2. Estimate the SNP density and depth for the region for TEs and control location to the 1st and 2nd nearest COs: perl assign_SNP_depth_nearest_CO_distance_flanking_TE_control_within_strains.pl

C. Analyze and plot the data: With_strain_analysis.R


IV. Between strains analysis
