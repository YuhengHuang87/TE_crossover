
1. generate variant files (vcf files) between the parental genome assemblies:
   sbatch align_Pacbio_genome_assemblies_SNP_calling.sub

2. map Illumina sequences from parental strains to check residual heterozygosity and map Illumina sequences from F2 individual fly:
   sbatch map_Illumina_reads_SNP_calling_pipline.sub

3. polish the SNP between parental genomes using their Illumina sequences:
   perl polish_vcf_snp_file.pl

4. identify CO in each F2 individual:
   perl call_CO_Illumina_genotype_parental_prop.pl
