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


References:

1.	Kellman, Benjamin P et al. “Elucidating Human Milk Oligosaccharide biosynthetic genes through network-based multi-omics integration.” Nature communications vol. 13,1 2455. 4 May. 2022, doi:10.1038/s41467-022-29867-4

2.	Johnson, Kelsey E et al. “Human milk variation is shaped by maternal genetics and impacts the infant gut microbiome.” Cell genomics vol. 4,10 (2024): 100638. doi:10.1016/j.xgen.2024.100638

3.	Ambalavanan, A., Chang, L., Choi, J. et al. “Human milk oligosaccharides are associated with maternal genetics and respiratory health of human milk-fed children.” Nat Commun 15, 7735 (2024). https://doi.org/10.1038/s41467-024-51743-6

4.	Spreckels J, Kurilshikov A, Fernández-Pato A et al. “Host and environmental determinants of human milk oligosaccharides and microbiota in the Lifelines NEXT cohort” Cell Reports, 2025; 44

5.	Hinrichs AS, Karolchik D, Baertsch R, Barber GP, Bejerano G, Clawson H, Diekhans M, Furey TS, Harte RA, Hsu F et al. The UCSC Genome Browser Database: update 2006. Nucleic Acids Res. 2006 Jan 1;34(Database issue):D590-8.

6.	Willer, Cristen J et al. “METAL: fast and efficient meta-analysis of genomewide association scans.” Bioinformatics (Oxford, England) vol. 26,17 (2010): 2190-1. doi:10.1093/bioinformatics/btq340

7.	Qiu, M., Gao, Z., Tian, M. et al. LIMCH1 inhibits antitumor immunity by upregulating PDL1 in triple-negative breast cancer. Cell Commun Signal 23, 419 (2025). https://doi.org/10.1186/s12964-025-02387-6

8.	Fang, Chao et al. “Unveiling Genetic Markers for Milk Yield in Xinjiang Donkeys: A Genome-Wide Association Study and Kompetitive Allele-Specific PCR-Based Approach.” International journal of molecular sciences vol. 26,7 2961. 25 Mar. 2025, doi:10.3390/ijms26072961

9.	Arun, Sondur J et al. “Targeted Analysis Reveals an Important Role of JAK-STAT-SOCS Genes for Milk Production Traits in Australian Dairy Cattle.” Frontiers in genetics vol. 6 342. 15 Dec. 2015, doi:10.3389/fgene.2015.00342

