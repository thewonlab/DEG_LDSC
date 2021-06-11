options(stringsAsFactors = FALSE)
outputdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/"

det = read.csv("/proj/hyejunglab/expression/capstone1_disorder/Combined_DTEsumstats_freeze1+2_4SVs_29seqPCs_LME_v2.csv")

##subset genes into down- and up-regulated 
sczdowngenes = deg[deg$SCZ.log2FC<0 & deg$SCZ.fdr<0.05,"X"]
sczupgenes = deg[deg$SCZ.log2FC>0 & deg$SCZ.fdr<0.05,"X"]  
sczdowngenes = unique(sczdowngenes[!is.na(sczdowngenes)])
sczupgenes = unique(sczupgenes[!is.na(sczupgenes)])

n_col=max(count.fields("/proj/hyejunglab/crossdisorder/annotation/AB.transcript.annot",sep="\t")) -2
snp_gene = read.table("/proj/hyejunglab/crossdisorder/annotation/AB.transcript.annot",sep="\t",fill=T,col.names=c("ENSID","RANGE",1:n_col))

##create SNP lists by identifying SNPs mapped to each gene set 
create_snp = function(gene) {
  genegrep = paste(gene, collapse="|")
  genesnp = snp_gene[grep(genegrep, snp_gene$ENSID),]
  genesnp = genesnp[,3:ncol(genesnp)]
  genesnplist = as.vector(t(genesnp))
  genesnplist = unique(genesnplist[genesnplist!=""])
  return (genesnplist)
}

sczdownsnplist = create_snp(sczdowngenes)
sczupsnplist = create_snp(sczupgenes)

sczlist = list(sczupsnplist, sczdownsnplist)
sczlistname = c("SCZ_DEG_up", "SCZ_DEG_down")

for (i in 1:22) {
  ##Read in locations of SNPs from the annotation files
  snps.table = read.table(gzfile(paste0("/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD.",i,".annot.gz")),header=T);
  
  indicator = matrix(0,nrow(snps.table),length(sczlist))
  
  for(k in 1:length(sczlist)){
    indicator[,k] <- ifelse(is.na(match(snps.table$SNP,sczlist[[k]])), 0, 1)
  }
  
  outframe = cbind(snps.table[,1:4], indicator)
  colnames(outframe)[5:ncol(outframe)] = sczlistname
  gz1 = gzfile(paste0(outputdir,"/",i,".annot.gz"), "w")
  write.table(outframe,gz1,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t");
  close(gz1);
  print(i)
}

