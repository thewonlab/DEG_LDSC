options(stringsAsFactors = F)
library(dplyr)
outputdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/kTotal_eQTL/"

data = read.csv("/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/nn.4399-S7.csv")
data = data[order(data$ktotal),]
numentries = nrow(data)
numgenes = ceiling(numentries/10)

kdatalist = vector("list", 10)
lowbound = 1

##subset genes based on kTotal
for(i in 1:10) {
  upbound = numgenes*i
  subdata = data[lowbound:upbound, ]
  kdata = subdata[,"Ensembl.ID"]
  kdata = unique(kdata[!is.na(kdata)])
  kdatalist[[i]] = kdata
  lowbound = upbound+1
}

# read in eQTL dataset 
genes = read.table("/proj/hyejunglab/eQTL/capstone4/eQTL/adultonly/capstone4.nom_eQTL.adult.significant.txt")

genes$V1 = gsub(pattern="\\..*",replacement="",genes$V1)
eQTL = select(genes, V1, V8)
colnames(eQTL) = c("gene", "snp")

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

ktotalGroup = paste0("ktotal", 1:10)

sczlist = list() 
for(i in 1:length(ktotalGroup)) {
  dir.create(paste0(outputdir, ktotalGroup[i]))
  sczlist[[i]] = create_snp(kdatalist[[i]])
}

for (i in 1:22) {
  ##Read in locations of SNPs from the annotation files
  snps.table = read.table(gzfile(paste0("/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD.",i,".annot.gz")),header=T);
  snplist = paste(snps.table$CHR, snps.table$BP, sep=":")
  
  indicator = matrix(0,nrow(snps.table),1)
  rownames(indicator) = snplist
  
  for(j in 1:length(ktotalGroup)) {
    colnames(indicator) = ktotalGroup[j]
    ktotaldir = paste0(outputdir, ktotalGroup[j])
    indicator[match(sczlist[[j]],rownames(indicator)),] = 1
    
    outframe = cbind(snps.table[,1:4], indicator)
    gz1 = gzfile(paste0(ktotaldir,"/",i,".annot.gz"), "w")
    write.table(outframe,gz1,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t");
    close(gz1);
  }
  print(i)
}

