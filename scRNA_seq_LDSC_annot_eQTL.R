options(stringsAsFactors = FALSE)
library(dplyr)
outputdir = "/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/scRNA_eQTL/"

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#read in deg up and down excel workbooks
deg_up = read_excel_allsheets("/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/scRNA_deg_data/scRNA_deg_up.xlsx")
deg_down = read_excel_allsheets("/proj/hyejunglab/crossdisorder/LDSC/partherit/SCZ_expression/scRNA_deg_data/scRNA_deg_down.xlsx")

load("/proj/hyejunglab/chr/geneAnno_allgenes.rda")

celltype = names(deg_up)

deglist = vector(mode="list", length=length(celltype))
#combine both up and down genes & extract genes
for(i in 1:length(celltype)){
  deg = rbind(deg_up[[i]], deg_down[[i]])
  degene = unique(deg$Gene)
  densg = geneAnno1[geneAnno1$hgnc_symbol %in% degene, "ensembl_gene_id"]
  densg = unique(densg[!is.na(densg) & densg!=""])
  deglist[[i]] = densg
}

##In-PV = Basket and Chandelier
in_pv = c(deglist[[13]], deglist[[14]])
deglist[[13]] = NULL
deglist[[13]] = NULL
deglist[[19]] = in_pv 

remove = c("In-PV (Basket)", "In-PV (Chandelier)")
celltype = celltype[!celltype %in% remove]
celltype[19] = "In-PV"

# read in eQTL dataset 
genes = read.table("/proj/hyejunglab/eQTL/capstone4/eQTL/adultonly/capstone4.nom_eQTL.adult.significant.txt")

genes$V1 = gsub(pattern="\\..*",replacement="",genes$V1)
eQTL = select(genes, V1, V8)
colnames(eQTL) = c("gene", "snp")

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

sczlist = list()

for(i in 1:length(celltype)) {
  dir.create(paste0(outputdir, celltype[i]))
  sczlist[[i]] = create_snp(deglist[[i]])
  
}

cellname = celltype 


for (i in 1:22) {
  ##Read in locations of SNPs from the annotation files
  snps.table = read.table(gzfile(paste0("/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD.",i,".annot.gz")),header=T);
  snplist = paste(snps.table$CHR, snps.table$BP, sep=":")
  
  indicator = matrix(0,nrow(snps.table),1)
  rownames(indicator) = snplist
  
  for(j in 1:length(cellname)) {
    colnames(indicator) = cellname[j]
    celldir = paste0(outputdir, cellname[j])
    indicator[match(sczlist[[j]],rownames(indicator)),] = 1
    
    outframe = cbind(snps.table[,1:4], indicator)
    gz1 = gzfile(paste0(celldir,"/",i,".annot.gz"), "w")
    write.table(outframe,gz1,quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t");
    close(gz1);
  }
  
  print(i)
}


