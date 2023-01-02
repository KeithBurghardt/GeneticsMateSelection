## Shared SNPs  

Here we show all traits associated with SNPs shared between married individuals who share slices.

Slices are denoted by CHR_ID (chomosome ID) plus slice_bp1 (basepair position of the start of the slice) and slice_bp2 (basepair position of end of the slice). The common allele is under the "base" column.

All features are extracted from the GWAS catalog, e.g., https://www.ebi.ac.uk/gwas/regions/chr10:122598339-123215850.

- AllSNPs.csv are shared SNPs with known traits among the full cohort
- AllSNPs_one_ethnic.csv are shared SNPs with known traits amont Europeans
- CommonSNPs_abs_cog27_all.zip show all shared SNPs in our cohort, where files are named "common_married_SNPs_[chromosome]_[beginning of slice basepair position]-[end of slice basepair position][OPTIONAL: 'one_ethnic' implies European sub-cohort]_cog27.csv"
    Important columns are 
     - chrom: SNP chromosome
     - snp_name: SNP name
     - cM: centiMorgan position
     - bp: basepair position
     - base: shared allele between married respondents

