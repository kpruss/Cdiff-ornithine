######################################

# Kali Pruss
# Pruss et al, Nature Metabolism 2021

# Code for in vivo transcriptional profiling differential gene expression for experiments with R20291 WT and R20291 Tox-
# Figs. 1c, 1e, Supplementary Table 3

######################################

## Input: takes raw counts table mapped to C. difficile R20291 genome
# rows are genes, with unique row IDs and
# columns are samples, with unique column names that match sample names in metadata provided

# Input: takes metadata dataframe (library_trmts) where first column (sample IDs) matches column headers in counts table and
# additional metadata is provided in subsequent columns, in this case Treatment is the comparison of interest

# Input: annotated genome dataframe where rows = genes that match to RNA-seq mapping

library(tidyverse)
library(DESeq2)

# Organize dataframes -----------------------------------------------------


# dataframe Cd8: reads mapped to R20291 genome
Cd8

# dataframe lib_trmts_08: metadata with sample names in rows, includes column with experimental group (R20291_Wt or R20291_Tox)
lib_trmts_08 

# R20291 annotation dataframe
anno.R2 

lib_trmts <- lib_trmts_08 # rename for ease

Cd8 <- Cd8[, -1] # remove Name column so retain only counts matrix
x <- Cd8 # rename for ease


# Run differential expression ---------------------------------------------


exp.all <- x
geneIDs = row.names(exp.all)
colData <- DataFrame(lib_trmts)
ds <- DESeqDataSetFromMatrix(countData = as.data.frame(exp.all),
                             colData = colData, design = ~Treatment)
ds 
dds <- estimateSizeFactors(ds)

# remove lowly-expressed genes
dds <- dds[ rowSums(counts(dds)) > 5, ] 

dds <- DESeq(dds)
res <- results(dds)

WtvTox <- results(dds, alpha = 0.01, contrast = c("Treatment", 
                                                  "R20291_Wt", "R20291_Tox")) 


# Annotate dataframe ------------------------------------------------------


list.res <- list(WtvTox)

# subset
Cd8_subset_p01 <- lapply(list.res, function(res){
  sub <- subset(res, padj < 0.01)
  df <- as.data.frame(sub)
})

# full dataframes
Cd8_dataframes <- lapply(list.res, function(res){
  df <- as.data.frame(res)
})

# subset
Cd8_subset_p05 <- lapply(list.res, function(res){
  sub <- subset(res, padj < 0.05)
  df <- as.data.frame(sub)
})


transform <- lapply(Cd8_dataframes, function(df){
  df$padj[is.na(df$padj)] <- 1
  df$neg.log.padj <- -log10(df$padj)
  df$log2FoldChange[is.na(df$log2FoldChange)] <- 0
  df$flip = 0
  idxs = df$log2FoldChange > 0
  df$flip[idxs] = df$neg.log.padj[idxs]
  df$flip[!idxs] = -1*df$neg.log.padj[!idxs]
  df$locus_tag <- row.names(df); df
})

# rename flip columns for condition
names(transform[[1]])[names(transform[[1]])=="flip"] <- "R2_WtvTox" # neg log padj
names(transform[[1]])[names(transform[[1]])=="log2FoldChange"] <- "R2_WtvTox_lfc" # log2 fold change
names(transform[[1]])[names(transform[[1]])=="lfcSE"] <- "R2_WtvTox_lfcSE" # log2 fold change standard error

R2 <- dplyr::select(transform[[1]], locus_tag, R2_WtvTox, R2_WtvTox_lfc, R2_WtvTox_lfcSE)

x <- R2 # x is dataframe to be annotated

anno <- function(x){
  mat <- match(x$locus_tag, anno.R2$Region)
  add <- anno.R2[mat, ]
  x$description <- add$product
  x$protID <- add$Protein.id
  x$locus_number <- add$locus_number
  x$gene_name <- add$Name
  x$LocusTag <- add$locus_tag
  x$R2_WtvTox_sig <- cut(x$R2_WtvTox, c(-100,-2,2,100), labels = c("down", "n.s.", "up")) # sig cut-off 0.01
  x
}

R2_full <- anno(x)
