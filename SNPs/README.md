## Shared SNPs  

Here we show all traits associated with SNPs shared between married individuals who share slices.

Slices are denoted by CHR_ID (chomosome ID) plus slice_bp1 (basepair position of the start of the slice) and slice_bp2 (basepair position of end of the slice). The common allele is under the "base" column.


- AllSNPs.csv are shared SNPs with known traits among the full cohort
- AllSNPs_one_ethnic.csv are shared SNPs with known traits amont Europeans
- SNP_GWAS_abs_cog27_all.zip contains all known traits associated with SNPs within a given slice. These features are extracted from the GWAS catalog, e.g., https://www.ebi.ac.uk/gwas/regions/chr10:122598339-123215850.
- CommonSNPs_abs_cog27_all.zip show all shared SNPs in our cohort, where files are named "common_married_SNPs_[chromosome]_[beginning of slice basepair position]-[end of slice basepair position][OPTIONAL: 'one_ethnic' implies European sub-cohort]_cog27.csv"
    Important columns are 
     - chrom: SNP chromosome
     - snp_name: SNP name
     - cM: centiMorgan position
     - bp: basepair position
     - base: shared allele between married respondents

- Code used in this analysis are in GeneticsMateSelection/Analyze Data/ExtractSNPs:
  - find_slices.py uses study_slices.py to find the common SNPs
  - Final processing of these data is with TraitsOfSNPs.py
  - SNP names come from: 'https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1319403889_9EQurcZePlSE3YnxbtdlokRG59T0&clade=mammal&org=Human&db=hg19&hgta_group=varRep&hgta_track=knownGene&hgta_table=0&hgta_regionType=range&position=chr'+chrom+'%3A'+http_bp1+'-'+http_bp2+'&hgta_outputType=primaryTable&hgta_outFileName='+new_slice+'.csv'
    - Replace chrom with the slice chromosome, and new_slice with the slice of the form [Chromosome]_[basepair position start]-[basepair position end], e.g., 10_122598339-123215850
    - http_bp1 = '%2C'.join([bp1[::-1][i:i+n][::-1] for i in range(0,len(bp1), n)][::-1])
    - http_bp2 = '%2C'.join([bp2[::-1][i:i+n][::-1] for i in range(0,len(bp2), n)][::-1])
  - SNP trait data comes from: 'https://www.ebi.ac.uk/gwas/regions/chr'+new_slice.replace('_',':')  
