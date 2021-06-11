library(data.table)
options(stringsAsFactors=F)
outputdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/kME_hiC/"

data = read.csv("/proj/hyejunglab/expression/capstone1_disorder/geneModule.csv")

modules = names(table(data$Module))

##get gene IDs and subset based on module number 
kmedat = list()
for(i in 1:length(modules)){
  module = modules[i]
  modnum = unlist(strsplit(module, "geneM"))[2]
  subdata = data[data$Module==module, ]
  kme = paste0("kME", modnum)
  subdata = subdata[,c("ensembl_gene_id", kme)]
  colnames(subdata) = c("gene","kME")
  kmedat = rbind(kmedat, subdata)
  print(i)
}

kmelist = vector("list", 10)

##Categorize genes based on kME 
for(i in 1:10) {
  cri1 = 0.1*i
  cri2 = cri1-0.1 
  
  kmesubdat = kmedat[kmedat$kME<cri1 & kmedat$kME>=cri2, ]
  kmelist[[i]] = unique(kmesubdat$gene)
  print(paste(i, length(kmelist[[i]])))
  
}

##read in Hi-C data
n_col=max(count.fields("/proj/hyejunglab/crossdisorder/annotation/AB_wointron.genes.annot",sep="\t")) -2
snp_gene = read.table("/proj/hyejunglab/crossdisorder/annotation/AB_wointron.genes.annot",sep="\t",fill=T,col.names=c("ENSID","RANGE",1:n_col))

##create SNP lists by identifying SNPs mapped to each gene set 
create_snp = function(gene) {
  genegrep = paste(gene, collapse="|")
  genesnp = snp_gene[grep(genegrep, snp_gene$ENSID),] #check whether ensembl_gene_ids are the same as the gene list, read annotation file 
  genesnp = genesnp[,3:ncol(genesnp)]
  genesnplist = as.vector(t(genesnp))
  genesnplist = unique(genesnplist[genesnplist!=""])
  return (genesnplist)
}

kmeGroup = paste0("kME",1:10)

for(i in 1:length(kmeGroup)) {
  dir.create(paste0(outputdir, kmeGroup[i]))
}

kme1 = kmelist[[1]]
kme1.1 = kme1[1:2000]
kme1.2 = kme1[2001:4000]
kme1.3 = kme1[4001:length(kme1)]
kme1.1snp = create_snp(kme1.1)
kme1.2snp = create_snp(kme1.2)
kme1.3snp = create_snp(kme1.3)
kme1snp = union(kme1.1snp, kme1.2snp)
kme1snp = unique(union(kme1snp, kme1.3snp))

kme2 = kmelist[[2]]
kme2.1 = kme2[1:2000]
kme2.2 = kme2[2001:length(kme2)]
kme2.1snp = create_snp(kme2.1)
kme2.2snp = create_snp(kme2.2)
kme2snp = unique(union(kme2.1snp, kme2.2snp))

kme3 = kmelist[[3]]
kme3.1 = kme3[1:1500]
kme3.2 = kme3[1501:length(kme3)]
kme3.1snp = create_snp(kme3.1)
kme3.2snp = create_snp(kme3.2)
kme3snp = unique(union(kme3.1snp, kme3.2snp))

kme4 = kmelist[[4]]
kme4.1 = kme4[1:1400]
kme4.2 = kme4[1401:length(kme4)]
kme4.1snp = create_snp(kme4.1)
kme4.2snp = create_snp(kme4.2)
kme4snp = unique(union(kme4.1snp, kme4.2snp))

kme5 = kmelist[[5]]
kme5.1 = kme5[1:1300]
kme5.2 = kme5[1301:length(kme5)]
kme5.1snp = create_snp(kme5.1)
kme5.2snp = create_snp(kme5.2)
kme5snp = unique(union(kme5.1snp, kme5.2snp))

kme6snp = create_snp(kmelist[[6]])
kme7snp = create_snp(kmelist[[7]])
kme8snp = create_snp(kmelist[[8]])
kme9snp = create_snp(kmelist[[9]])
kme10snp = create_snp(kmelist[[10]])

sczlist = list(kme1snp, kme2snp, kme3snp, kme4snp, kme5snp, kme6snp, kme7snp, kme8snp, kme9snp, kme10snp)

for (i in 1:22) {
  ##Read in locations of SNPs from the annotation files
  snps.table = read.table(gzfile(paste0("/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD.",i,".annot.gz")),header=T);
  
  indicator = matrix(0,nrow(snps.table),1)
  rownames(indicator) = snps.table$SNP
  
  for(j in 1:length(kmeGroup)) {
    colnames(indicator) = kmeGroup[j]
    kmedir = paste0(outputdir, kmeGroup[j])
    indicator[match(sczlist[[j]],rownames(indicator)),] = 1
    
    outframe = cbind(snps.table[,1:4], indicator)
    gz1 = gzfile(paste0(kmedir,"/",i,".annot.gz"), "w")
    write.table(outframe,gz1,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t");
    close(gz1);
  }
  print(i)
}




