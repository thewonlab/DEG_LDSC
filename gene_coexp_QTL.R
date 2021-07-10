####### Generation of annot files #######
options(stringsAsFactors = FALSE)
library(GenomicRanges)
outputdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/DEG_coexp_QTL/"
setwd(outputdir)
load("/proj/hyejunglab/chr/geneAnno_allgenes.rda")
load("/proj/hyejunglab/eQTL/capstone4/eQTL/adultonly/eQTL_GRanges_colocinput.rda")
qtlranges$transcript = unlist(lapply(strsplit(qtlranges$gene, split="[.]"),'[[',1))

deg = read.csv("/proj/hyejunglab/expression/capstone1_disorder/geneModule.csv")
isoforms = paste0("geneM",1:34)
genes = vector(length=length(isoforms),mode="list")
for(i in 1:length(isoforms)){
  isoform = isoforms[i]
  genes[[i]] = unique(deg[deg$Module==isoform,"ensembl_gene_id"])
  dir.create(paste0(outputdir, isoform))
}

for (i in 1:22) {
  ##Read in locations of SNPs from the annotation files
  snps.table = read.table(gzfile(paste0("/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD.",i,".annot.gz")),header=T);
  
  snpranges = GRanges(snps.table$CHR, IRanges(snps.table$BP,snps.table$BP))
  
  for(j in 1:length(isoforms)){
    isoform = isoforms[j]
    isodir = paste0(outputdir, isoform)
    genelist = genes[[j]]
    
    degranges = qtlranges[qtlranges$transcript %in% genelist]
    olap = findOverlaps(snpranges, degranges)
    outframe = snps.table[,1:4]
    outframe$deg = 0
    outframe[unique(queryHits(olap)),"deg"] <- 1
    gz = gzfile(paste0(isodir,"/",i,".annot.gz"), "w")
    write.table(outframe,gz,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t");
    close(gz);
  }
  print(i)
}

####### Run LDSC #######
options(stringsAsFactors=F)
setwd("/proj/hyejunglab/crossdisorder/LDSC/sumstat/")
bashout = "/proj/hyejunglab/crossdisorder/LDSC/partherit/bash_temp/"
outdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/DEG_coexp_QTL/"
condition = "scz3"
isoforms = paste0("geneM",1:34)

for(i in 1:length(isoforms)){
  isoform = isoforms[i]
  annotdir = paste0(outdir, isoform, "/")
  
  system(paste0("echo '#!/bin/tcsh' > ",bashout, condition,".tcsh"))
  system(paste0("echo '/proj/hyejunglab/program/ldsc/ldsc.py ",
                "--h2 /proj/hyejunglab/crossdisorder/LDSC/sumstat/", condition,".sumstats.gz ",
                "--out ",annotdir, condition, "_results ",
                "--frqfile-chr /proj/hyejunglab/program/ldsc/LDSC/1000G_Phase3_frq/1000G.EUR.QC. ",
                "--overlap-annot --ref-ld-chr ",annotdir,",/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD. ",
                "--w-ld-chr /proj/hyejunglab/program/ldsc/LDSC/weights_hm3_no_hla/weights.",
                "' >> ",bashout, condition,".tcsh"))  
  system(paste0("sbatch -n 1 --mem=30g -o ", bashout, condition, ".out ", bashout, condition,".tcsh"))
}

####### Read results #######
options(stringsAsFactors=F)
setwd("/proj/hyejunglab/crossdisorder/LDSC/sumstat/")
bashout = "/proj/hyejunglab/crossdisorder/LDSC/partherit/bash_temp/"
outdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/DEG_coexp_QTL/"
isoforms = paste0("geneM",1:34)

heritability = c()
for(i in 1:length(isoforms)){
  isoform = isoforms[i]
  annotdir = paste0(outdir, isoform, "/")
  setwd(annotdir)
  
  scz = read.table("scz3_results.results", header=T)
  heritability = rbind(heritability, scz[1,])
  
  print(i)
}

heritability$Category = isoforms
heritability$FDR = p.adjust(heritability$Enrichment_p, "BH")
heritability[heritability$FDR<0.05,]
heritability$logp = -log10(heritability$Enrichment_p)

####### Plot results #######
library(ggplot2)
library(data.table)
library(gtools)
library(grid)
library(gridExtra)

setwd(outdir)
isoAsso = read.csv("/proj/hyejunglab/expression/capstone1_disorder/module-disease_association.csv")
isoAsso = isoAsso[grep("geneM", isoAsso$ModulName),]
isoAsso = isoAsso[isoAsso$Group=="SCZ",]
isoAsso$fdr[isoAsso$fdr>0.05] = "NA"
isoAsso$ModulName = factor(isoAsso$ModulName,levels=rev(isoforms))

heritability$Category = factor(heritability$Category,levels=rev(isoforms))

plot1 = ggplot(isoAsso,aes(Group,ModulName)) + geom_tile(aes(fill=beta), colour="grey") + coord_fixed(ratio=1) + scale_fill_gradient2(low="steelblue", mid="white", high="indianred", midpoint=0) + ylab("Isoform Module")
plot2 = ggplot(heritability,aes(x=Category,y=logp, fill=Enrichment)) + geom_bar(stat="identity",position=position_dodge(),color="grey") + geom_hline(yintercept=2,lty=2) + coord_flip() + scale_fill_gradient2(low="steelblue", mid="white", high="indianred",midpoint=0)

pdf("gene_coexp_PEC_QTL.pdf", height=10, width=10)
grid.arrange(plot1, plot2, nrow = 1, widths=c(1,2))
dev.off()