1. Map the CUT&Tag Illumina reads and estimate the genome-wide H3K9me3 fragment: sbatch cut_tag_mapping_bowtie2_bed.sub
2. Estimate local H3K9me3 level 20 - 40 kb from TEs as controls: perl fragment_CUT_Tag_bedgraph_local_control.pl
3. Estimate 5 kb adjacent K9 enrichment around TEs: perl adjacent_H3K9me3_enrichment_estimate.pl
4. Estimate the total K9 mass around TEs: perl H3K8me3_mass_estimate.pl
5. Average K9 across alleles with TEs and homologous alleles without TEs between strains: perl collate_plot_avg_TEs_alternative_strain_mean_median.pl
6. plot the results: CUT_TAG.R
