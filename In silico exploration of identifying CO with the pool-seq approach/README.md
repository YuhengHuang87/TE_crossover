I. explore the false positive rate with in silico simulation
1. sbatch PBSIM_minimap_FP.sub
2. run the Identify CO in PacBio pool-seq data pipeline to see if any reads were called falsely (starting step 2).

II. Test the impact of sequencing, depth, and structural variants on the false negative CO recall rate
1. generate recombinant haplotypes (one CO per chromosome arm): perl generate_one_recombiant_individual_chro_arm.pl
2. run PBSIM to simulate PacBio long reads, extract those read with obtain_PBSIM_reads_recombinant.pl, and map those CO reads: sbatch PBSIM_minimap_False_negative.sub
3. run the pipeline as above to recall CO events using our approach
4. perform downsampling on the identified CO reads to explore the impacts of read depth on recall rate, e.g., shuf -n <<number of reads>> all_15_co_gc_shared_each_read_fragment_cutoff_2000_10_10_10_2 > ${depth}_co_gc_shared_each_read_fragment_cutoff_2000_10_10_10_2
5. calculate recall rate per event and binning them: perl percentage_recombinant_events_recall.pl
6. Identify reads overlap with TEs or SVs (master_table_d10_100bp.txt contains SV profiles, obtained from Chakraborty et al. 2019): perl classify_PBSIM_CO_reads_contrain_TE_SV.pl
7. check recall rates for reads w/ and w/o TEs or SVs: perl percentage_pbsim_CO_reads_recalled_TE_SV_within.pl

