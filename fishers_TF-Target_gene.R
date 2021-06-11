
load("~/PEC_TF_E_P_Target_cutoff=0.1_elastic_net.RData")
sczrisk = read.csv("~/ppi_interactors.csv")

## define schizophrenia-associated TFs by identifying TFs out of high-confidence schizophrenia risk genes
transcript_factor = intersect(edgelist$RF, sczrisk$x)
target_genes = edgelist[edgelist$RF %in% transcript_factor, "TG"]
target_genes = unique(target_genes)

#psychencode DEG
gene_list = read.csv("~/Combined_DGEsumstats_freeze1+2_4SVs_29seqPCs_LME.csv")

##create subsets of genes to create transcription factor-target gene linkages 
pc_genes = gene_list[gene_list$gene_type == "protein_coding", "gene_name"]
deg = gene_list[gene_list$gene_type == "protein_coding" & gene_list$SCZ.fdr < 0.05, "gene_name"]
non_deg = gene_list[gene_list$gene_type == "protein_coding" & gene_list$SCZ.fdr > 0.05, "gene_name"]
deg_up = gene_list[gene_list$gene_type == "protein_coding" & gene_list$SCZ.fdr < 0.05 & gene_list$SCZ.log2FC>0, "gene_name"]
non_deg_up = setdiff(pc_genes, deg_up)
deg_down = gene_list[gene_list$gene_type == "protein_coding" & gene_list$SCZ.log2FC<0 & gene_list$SCZ.fdr < 0.05 , "gene_name"]
non_deg_down = setdiff(pc_genes, deg_down)

non_targetgenes = setdiff(pc_genes, target_genes)

#creating matrix and running Fisher's Exact Test 
run_fishers_protein = function(deg_dat, non_deg_dat) {
  deg_tf = intersect(deg_dat, target_genes)
  non_deg_tf = intersect(non_deg_dat, target_genes)
  deg_non_tf = intersect(deg_dat, non_targetgenes)
  non_degtf = setdiff(pc_genes, union(deg_dat, target_genes))
  
  A = length(deg_tf)
  B = length(non_deg_tf)
  C = length(deg_non_tf)
  D = length(non_degtf)
  
  return(fisher.test(matrix(c(A,B,C,D),2,2)))
  
}


corr_deg_tf = run_fishers_protein(deg, non_deg)
corr_deg_up_tf = run_fishers_protein(deg_up, non_deg_up)
corr_deg_down_tf = run_fishers_protein(deg_down, non_deg_down)





                     