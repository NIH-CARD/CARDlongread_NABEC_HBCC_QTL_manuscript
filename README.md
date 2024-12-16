# **Long-read sequencing of hundreds of diverse brains provides insight into the impact of structural variation on gene expression and DNA methylation**

![Sample Image](Workflow.png)

CARD ❤️ Open Science 😍

**Written By:** *Kimberley Billingsley, *Melissa Meredith, *Kensuke Daida

**Last Updated**: December 2024 

**Quick Description**: Evaluating the impact of SNVs and SVs on gene expression and methylation in the human brain using long-read sequencing data. 

**Link to Manuscript:** 

### **Summary:**
This is the online repository for the article titled "Long-read sequencing of hundreds of diverse brains provides insight into the impact of structural variation on gene expression and DNA methylation". Structural variants (SVs) drive gene expression in the human brain and are causative of many neurological conditions. However, most existing genetic studies have been based on short-read sequencing methods, which capture fewer than half of the SVs present in any one individual. Long-read sequencing (LRS) enhances our ability to detect disease-associated and functionally relevant structural variants (SVs); however, its application in large-scale genomic studies has been limited by challenges in sample preparation and high costs. Here, we leverage a new scalable wet-lab protocol and computational pipeline for whole-genome Oxford Nanopore Technologies sequencing and apply it to neurologically normal control samples from the North American Brain Expression Consortium (NABEC) (European ancestry) and Human Brain Collection Core (HBCC) (African or African admixed ancestry) cohorts. Through this work, we present a publicly available long-read resource from 351 human brain samples (median N50: 27 Kbp and at an average depth of ~40x genome coverage). We discover approximately 250,000 SVs and produce locally phased assemblies that cover 95% of all protein-coding genes in GRCh38. Utilizing matched expression datasets for these samples, we apply quantitative trait locus (QTL) analyses and identify SVs that impact gene expression in post-mortem frontal cortex brain tissue. Further, we determine haplotype-specific methylation signatures at millions of CpGs and, with this data, identify cis-acting SVs. In summary, these results highlight that large-scale LRS can identify complex regulatory mechanisms in the brain that were inaccessible using previous approaches. We believe this new resource provides a critical step toward understanding the biological effects of genetic variation in the human brain. 

### **Repository Orientation:**


**1.Variant calling and QC** 
- 1.1 Sv_merging.py: script to localize and merge SVs on Terra using Sniffles and Truvari. 
- 1.2 truvari_resolve_symbolic_svs.py: script to replace symbolic SVs with reference sequence for merging. 
- 1.3 sv_figure_ploting.py: script to plot figures describing the merged SV set.


**2. eQTL analyses** 
- 2.1 prep_quants_RNA.ipynb:  notebook to format expression data for tensorQTL
- 2.2 cis_eqtl_tensorqtl.ipynb:  notebook to run tensorQTL


**3. Methyaltion analyses**

- 3.1 unphased_bedmethyl_merge.py: localizes and merges aggregated methylation for samples in each cohort. 
- 3.2 phased_methylation_merge_funcitons.py: functions used in the merging process
- 3.3 GenomeWindowingMethylation_and_AgeRegressions_COV_gender_PMI.py: Merging aggregated windows of methylation across the genome, running linear regressions with age and covarites, plotting methylation figure images. 

**4.  mQTL analyses**

- 4.1  cis_mqtl_tensorqtl.ipynb:  notebook to format methylation data for tensorQTL
- 4.2  cis_mqtl_tensorqtl.ipynb:  notebook to run tensorQTL
