####### Generation of annot files #######
options(stringsAsFactors = FALSE)
library(GenomicRanges)
outputdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/CMC_coexp_QTL/"
setwd(outputdir)
load("/proj/hyejunglab/eQTL/capstone4/eQTL/adultonly/eQTL_GRanges_colocinput.rda")
qtlranges$transcript = unlist(lapply(strsplit(qtlranges$gene, split="[.]"),'[[',1))

deg = read.csv("/proj/hyejunglab/expression/SCZ_RNAseq/SCZ_module_CMC.csv")
isoforms = paste0("M",1:35,"c")
genes = vector(length=length(isoforms),mode="list")
for(i in 1:length(isoforms)){
  isoform = isoforms[i]
  genes[[i]] = unique(deg[deg$Module==isoform,"Ensembl.ID"])
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
  print(paste0("chr",i))
}

####### Run LDSC #######
options(stringsAsFactors=F)
setwd("/proj/hyejunglab/crossdisorder/LDSC/sumstat/")
bashout = "/proj/hyejunglab/crossdisorder/LDSC/partherit/bash_temp/"
outdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/CMC_coexp_QTL/"
condition = "scz3"
isoforms = paste0("M",1:35,"c")

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
outdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/CMC_coexp_QTL/"
isoforms = paste0("M",1:35,"c")

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

