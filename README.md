Project Overview:

Genome-wide association study (GWAS) summary statistics from the MILK (n = 540), CHILD (n = 980), and Lifelines NEXT (n = 524) cohorts were harmonized and meta-analyzed using METAL to identify genetic determinants of HMO composition.


The analysis includes:

- Harmonization of GWAS summary statistics across cohorts
- Genome build conversion (GRCh37 to GRCh38)
- Meta-analysis using METAL
- Generation of summary tables
- Manhattan, QQ, forest, and locus zoom plots
- Downstream visualization and result summarization in R


Results:


Implications:

In total, eight loci were associated with HMO composition, including three putatively novel associations. 

These findings demonstrate that meta-analysis increases statistical power to identify genetic variants influencing HMO composition and provide candidate loci for future functional investigation.


Limitations:

Identified associations are statistical and do not establish causal relationships or identify the underlying causal variants or genes.

Biological mechanisms linking identified genetic loci to HMO biosynthesis remain unresolved.

Differences in cohort characteristics, sample sizes, and phenotype measurement protocols may have introduced residual heterogeneity, and due to the predominantly European ancestry of the participating cohorts, the generalizability of these findings to other ancestral populations remains to be determined.



Reproducibility:

The GWAS summary statistics used in this project originate from multiple independent cohorts and are subject to data-sharing restrictions. Consequently, the raw input data cannot be made publicly available through this repository.

Because the original summary statistics are unavailable, the complete analysis cannot be reproduced directly from this repository. However, all analysis scripts, plotting code, and workflow used to process the data and generate the reported results are provided for transparency and as a reference for similar meta-analysis projects.

Researchers with access to comparable GWAS summary statistics should be able to adapt the workflow to their own data.
