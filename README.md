# DEG_LDSC

Authors: Alice Yu and Hyejung Won 

Created: 7/11/2021 

Updated: 7/20/2021

Each folder contains Stratified LD Score Regression (S-LDSC) annotation file for each transcriptomic feature. 

Each folder contains two subdirectories: 
  - Hi-C: annotation file generated based on Hi-C evidence 
  - QTL: annotation file generated based on QTL evidence 

## Cell-type specific DEG
  - Reference: Ruzicka, W.B. et al.(2020) medRxiv doi: https://doi.org/10.1101/2020.11.06.20225342

## DEG/DET
  - Reference: Gandal, M.J. et al.(2018) PMID: 30545856
  - Differentially expressed genes (DEGs) or differentially expressed transcriptions (DETs) with FDR < 0.05 were stratified into up- and downregulated genes based on log2 fold change (logFC)

## Gene Co-expression Modules
  - PsychENCODE Consortium (PEC): Gandal, M.J. et al.(2018) PMID: 30545856
  - CommonMind consortium (CMC): Fromer, M. et al.(2016) PMID: 27668389

## Isoform-level Co-expression Modules
  - Reference: Gandal, M.J. et al.(2018) PMID: 30545856

## kME/kTotal
  - kME (PEC): Gandal, M.J. et al.(2018) PMID: 30545856
  - kTotal (CMC): Fromer, M. et al.(2016) PMID: 27668389
  - Genes classified into 10 groups based on centrality/connectivity values 

## Transcription Factor-Target Gene
  - Reference: Wang, D. et al.(2018) PMID: 30545857
### fishers_TF-Target_genes.R
  - Script runs Fisher's Exact Test with protein-coding DEGs and target genes of schizophrenia-associated transcription factors (TFs) 
  - Schizophrenia-associated TFs were defined using TFs of high-confidence schizophrenia risk genes from PEC 

## .results
  - Contains heritability enrichment and significance of each transcriptomic feature
  - Cols
    - Prop._SNPs: Proportion of SNPs
    - Prop._h2: Proportion of heritability
    - Prop._h2_std_error: Standard Error of heritability
    - Enrichment: Heritability Enrichment
    - Enrichment_std_error: Standard Error of Heritability Enrichment
    - Enrichment_p: P-value of enrichment 

## Functional Genomic Data
  - Hi-C, eQTL, and isoQTL: Wang, D. et al.(2018) PMID: 30545857
     - Used data generated from dorsolateral prefrontal cortex (DLPFC) 


## Reference
Please cite this paper: Yu, A.W.; Peery, J.D.; Won, H. Limited Association between Schizophrenia Genetic Risk Factors and Transcriptomic Features. Genes 2021, 12, 1062. https://doi.org/10.3390/genes12071062

SCZ GWAS: Pardinas, A.F. et al.(2018) PMID: 29483656

S-LDSC (v1.0.1): Finucane, H.K. et al.(2015) PMID: 26414678 
  - Baseline model (v1.1.0) was used to assess heritability enrichment

H-MAGMA input file: Sey, N.Y.A. et al.(2020) PMID: 32152537
  - Used to identify SNPs mapped to each gene 

