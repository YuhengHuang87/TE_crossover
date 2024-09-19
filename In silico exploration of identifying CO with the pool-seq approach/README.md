I. explore the false positive rate with in silico simulation
1. sbatch PBSIM_minimap_FP.sub
2. run the Identify CO in PacBio pool-seq data pipeline to see if any reads were called falsely (starting step 2).

II. Test the impact of sequencing, depth, and structural variants on the false negative CO recall rate
1. generate recombinant haplotypes (one CO per chromosome arm): perl generate_one_recombiant_individual_chro_arm.pl
2. run PBSIM to simulate PacBio long reads, extract those read with obtain_PBSIM_reads_recombinant.pl, and map those CO reads: sbatch PBSIM_minimap_False_negative.sub
