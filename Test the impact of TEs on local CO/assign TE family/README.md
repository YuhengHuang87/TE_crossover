1. Run repeatmodeler to annotate TEs and Blast TE sequences to TE family information (dmel-all-transposon-r6.32_canonicalBLAST, publicly available TE family information for D. mel): sbatch repeatmodeler.sub
2. Assign TE family based on Blast output: perl assign_TE_family.pl
3. Merge or exclude nearby TEs based on family information: perl TE_nearby_distance_merge_by_family.pl
