options(stringsAsFactors = F)
library(dplyr)

genes = read.table("/proj/hyejunglab/eQTL/capstone4/eQTL/adultonly/capstone4.nom_eQTL.adult.significant.txt")

genes$V1 = gsub(pattern="\\..*",replacement="",genes$V1)
eQTL = select(genes, V1, V8)
colnames(eQTL) = c("gene", "snp")

outputdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/eQTL"

##subset genes into down- and up-regulated 
deg = read.csv("/proj/hyejunglab/expression/capstone1_disorder/Combined_DGEsumstats_freeze1+2_4SVs_29seqPCs_LME.csv")
sczdowngenes = deg[deg$SCZ.log2FC<0 & deg$SCZ.fdr<0.05,"gene_id"]
sczupgenes = deg[deg$SCZ.log2FC>0 & deg$SCZ.fdr<0.05,"gene_id"]  
sczdowngenes = unique(sczdowngenes[!is.na(sczdowngenes)])
sczupgenes = unique(sczupgenes[!is.na(sczupgenes)])

##create SNP lists by identifying SNPs mapped to each gene set 
create_snp = function(gene) {
  genesnplist = c()
  for(i in 1:length(gene)){
    genegrep = gene[i]
    genesnp = eQTL[eQTL$gene==genegrep,"snp"] 
    genesnplist = c(genesnplist, genesnp)
    genesnplist = unique(genesnplist[genesnplist!=""])
  }
  return (genesnplist)
}

sczdownsnplist = create_snp(sczdowngenes)
sczupsnplist = create_snp(sczupgenes) 


sczlist = list(sczupsnplist, sczdownsnplist)
sczlistname = c("SCZ_DEG_up", "SCZ_DEG_down")

for (i in 1:22) {
  ##Read in locations of SNPs from the annotation files
  snps.table = read.table(gzfile(paste0("/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD.",i,".annot.gz")),header=T);
  snplist = paste(snps.table$CHR, snps.table$BP, sep=":")
  
  indicator = matrix(0,nrow(snps.table),length(sczlist))
  
  for(k in 1:length(sczlist)){
    indicator[,k] <- ifelse(is.na(match(snplist,sczlist[[k]])), 0, 1)
  }
  
  outframe = cbind(snps.table[,1:4], indicator)
  colnames(outframe)[5:ncol(outframe)] = sczlistname
  gz1 = gzfile(paste0(outputdir,"/",i,".annot.gz"), "w")
  write.table(outframe,gz1,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t");
  close(gz1);
  print(i)
}