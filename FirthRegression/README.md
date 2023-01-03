# Firth logistic regression
- p-value files are of the form: firth_[coefficient added]_prob_full_fold=0_abs_no_cog-bmi_cog27.csv
- Coefficient files are of the form: firth_[coefficient added]_coefs_full_fold=0_abs_no_cog-bmi_cog27.csv

Each file will have the following columns:

- chrom: chromosome of slice
- bp1: beginning basepair of slice 
- bp2: end basepair of slice	
- (Intercept): p-value or coefficient of intercept of Firth logistic regression	
- slice: p-value or coefficient of slice with Firth logistic regression
