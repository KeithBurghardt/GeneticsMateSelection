# Code to replicate analysis from "Shared Genetics Drive Mate Selection"

- ILASH: this is code to run ILASH on each chromosome of the HRS dataset. The code to run ILASH can be found here: https://github.com/roohy/iLASH
- Analyze Data: this is code to extract features and find the Firth logistic regression for each genetic slice. The code is decribed in more detail below

## Code to analyze data (Analyze Data folder)
- create_pairwise_features.py: This code extracts features from pairs of participants.
- covariates.py: This code finds the correlations between features for each married participant, and creates some of the plots in the paper, such as PCA embeddings of participants and mean PC distance between married and random participants.
- prepare_firth_regression.py: This code collects and cleans pairwise features extracted in create_pairwise_features.py. This cleaned data is fed into firth_fit.r, which is run through this code
- firth_fit.r: This code extracts logistic regression coefficients and p-values for each slice, both with and without controlling for covariates. 
- prs_score.py: the Slice Polygeneic Relationship Score (SliPRS) for each pair of participants


## Software requirements

### The following software were used for data analysis:
- Python version 3.7.6 
- R version 4.0.0

### The following packages were used for data analysis: 
Python packages:
- os,random,glob,pickle,multiprocessing (packages come with Python)
- numpy 1.21.4
- pandas 1.0.5
- scipy 1.2.0

R packages:
- (logistf)[https://cran.r-project.org/web/packages/logistf/index.html] 1.24.1
- data.table (package comes with R)


## License

The project is licensed under the [BSD 3-Clause license](https://github.com/pysal/spaghetti/blob/main/LICENSE.txt).


## BibTeX Citation
```
@article{Burghardt2022,
    year      = {2022},
    publisher = {},
    author    = {Keith Burghardta and Thalida Em Arpawongb and José Luis Ambite},
    title     = {Shared Genetics Drive Mate Selection},
    journal   = {}
}
```


## Funding
The HRS (U.S. Health and Retirement Study) is sponsored by the National Institute on Aging through a cooperative agreement (grant number NIA U01AG009740) The datasets are produced by the University of Michigan, Ann Arbor. The HRS phenotypic data files are public use datasets (https://hrs.isr.umich.edu). The HRS genotype data is available to authorized researchers (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000428.v2.p2). This work funded J.L.A. and K.B. under NIH grant R01HG010297 (PIs: T.C. Matise, C.R. Gignoux) and T.E.A under NIH grant P30AG017265 (PI: E. Crimmins).
