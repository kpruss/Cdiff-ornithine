######################################
# Kali Pruss
# Pruss et al, Nature Metabolism 2021

# Code for in vivo transcriptional profiling differential gene expression for experiments with 630 WT and 630 Tox-
# Figs. 1c-e, Supplementary Table 2

#####################################

## Input: takes raw counts table mapped to C. difficile 630 genome
# rows are genes, with unique row IDs and
# columns are samples, with unique column names that match sample names in metadata provided

# Input: takes metadata dataframe (library_trmts) where first column (sample IDs) matches column headers in counts table and
# additional metadata is provided in subsequent columns, in this case Treatment is the comparison of interest

# Input: annotated genome dataframe where rows = genes that match to RNA-seq mapping

library(tidyverse)
library(DESeq2)

# Organize dataframes -----------------------------------------------------


# dataframe Cd1: reads mapped to 630 genome in the presence of defined community 1
Cd1

# dataframe lib_trmts_01: metadata with sample names in rows, includes column with experimental group (630_Wt or 630_Tox)
lib_trmts_01

# dataframe Cd2: reads mapped to 630 genome in the presence of defined community 1
Cd2

# dataframe lib_trmts_01: metadata with sample names in rows, includes column with experimental group (630_Wt or 630_Tox)
lib_trmts_02

# add defined community to metadata
lib_trmts_01$Def_Comm <- "DC1"
lib_trmts_02$Def_Comm <- "DC2"

# remove control group from experiment 1
lib_trmts_01 <- lib_trmts_01[-c(1:4), ]

# change metadata to match lib_trmts_02 (better to distinguish 630 vs  R20291 anyway)
levels(as.factor(lib_trmts_01$Treatment))

# make factor so that can re-name to match
lib_trmts_01$Treatment <- as.factor(lib_trmts_01$Treatment)
lib_trmts_01$Library <- as.factor(lib_trmts_01$Library)

levels(lib_trmts_01$Treatment)[1] <-  "630_Tox"
levels(lib_trmts_01$Treatment)[2] <- "630_Wt"


# rename libraries that have overlapping names in metadata
levels(lib_trmts_01$Library)
levels(lib_trmts_01$Library)[3] <- "KP120"
levels(lib_trmts_01$Library)[4] <- "KP130"
lib_trmts_01
lib_trmts_12 <- rbind(lib_trmts_01, lib_trmts_02)

# dataframe containing metadata for RNA-seq library sample preparation, Treatment condition, and defined community
lib_trmts_12


# rename samples in count data
Cd1 <- Cd1 %>% dplyr::rename("KP120" = "KP12")
Cd1 <- Cd1 %>% dplyr::rename("KP130" = "KP13")

Cd2 <- Cd2 %>% dplyr::rename("ID" = "Name")
head(Cd1)
head(Cd2)

# merge counts matrix from two independent experiments
Cd12 <- merge(Cd1, Cd2, by = "Identifier")

# view rRNAs
Cd12[grep("_rRNA_", Cd12$ID.x), ]
# remove them 
Cd12_norRNAs <- Cd12[-grep("_rRNA_", Cd12$ID.x), ]
row.names(Cd12_norRNAs) <- Cd12_norRNAs$Identifier

# retain gene IDs as rownames, counts (rows) by sample IDs (columnes) only 
Cd12_norRNAs <- Cd12_norRNAs %>% dplyr::select(-c(Identifier,ID.x,ID.y))


# assign libraries (Sample IDs) to rownames to match column names for DESeq pipeline
row.names(lib_trmts_12) <- lib_trmts_12$Library # need to do this explicity from removing levels, etc rownames get wonky

# ensure concordance
match(rownames(lib_trmts_12), colnames(Cd12_norRNAs))

x <- Cd12_norRNAs # rename for ease


# create new unique treatment_defined community metadata identifier
# Group_Defined Community
lib_trmts_12 <- lib_trmts_12 %>% unite(col = Group, 
                                       Treatment, Def_Comm, sep = "_", remove = FALSE)
lib_trmts <- lib_trmts_12 # rename for ease



# Run differential expression ---------------------------------------------


## by +/- toxin status
exp.all <- x
geneIDs = row.names(exp.all)
colData <- DataFrame(lib_trmts)
ds <- DESeqDataSetFromMatrix(countData = as.data.frame(exp.all),
                             colData = colData, design = ~Treatment) # run by C. diff infection (i.e. combined)
ds 
dds <- estimateSizeFactors(ds)

# filter lowly-expressed genes
dds <- dds[ rowSums(counts(dds)) > 8, ] 

dds <- DESeq(dds) 
res <- results(dds)


WtvTox_both_res <- results(dds, alpha = 0.01, contrast = c("Treatment", "630_Wt", "630_Tox")) 

## by experimental group
ds2 <- DESeqDataSetFromMatrix(countData = as.data.frame(exp.all),
                              colData = colData, design = ~Group) # run by experimental  group
ds2
dds2 <- estimateSizeFactors(ds2)
dds2 <- dds2[ rowSums(counts(dds2)) > 8, ] 
dds2 <- DESeq(dds2)
res2 <- results(dds2)

DC1_WtvTox_res <- results(dds2, alpha = 0.01, contrast = c("Group", "630_Wt_DC1", "630_Tox_DC1"))
DC2_WtvTox_res <- results(dds2, alpha = 0.01, contrast = c("Group", "630_Wt_DC2", "630_Tox_DC2")) 



# Annotate dataframes -----------------------------------------------------



list.res <- list(DC1_WtvTox_res, DC2_WtvTox_res, WtvTox_both_res)

# subset
Cd12_subset_p01 <- lapply(list.res, function(res){
  sub <- subset(res, padj < 0.01)
  df <- as.data.frame(sub)
})

# full dataframes
Cd12_dataframes <- lapply(list.res, function(res){
  df <- as.data.frame(res)
})


# now iterating over the full dataframes
transform <- lapply(Cd12_dataframes, function(df){
  df$padj[is.na(df$padj)] <- 1
  df$neg.log.padj <- -log10(df$padj)
  df$log2FoldChange[is.na(df$log2FoldChange)] <- 0
  df$flip = 0
  idxs = df$log2FoldChange > 0
  df$flip[idxs] = df$neg.log.padj[idxs]
  df$flip[!idxs] = -1*df$neg.log.padj[!idxs]
  df$locus_tag <- row.names(df); df
})


# rename flip columns by condition
names(transform[[1]])[names(transform[[1]])=="flip"] <- "Cd630_v_Tox_DC1" # neg log padj
names(transform[[2]])[names(transform[[2]])=="flip"] <- "Cd630_v_Tox_DC2"
names(transform[[3]])[names(transform[[3]])=="flip"] <- "Cd630_v_Tox_both"

names(transform[[1]])[names(transform[[1]])=="log2FoldChange"] <- "Cd630_v_Tox_DC1_lfc" # log2 fold change
names(transform[[2]])[names(transform[[2]])=="log2FoldChange"] <- "Cd630_v_Tox_DC2_lfc"
names(transform[[3]])[names(transform[[3]])=="log2FoldChange"] <- "Cd630_v_Tox_both_lfc"

names(transform[[1]])[names(transform[[1]])=="lfcSE"] <- "Cd630_v_Tox_DC1_lfcSE" # log2 fold change standard error
names(transform[[2]])[names(transform[[2]])=="lfcSE"] <- "Cd630_v_Tox_DC2_lfcSE"
names(transform[[3]])[names(transform[[3]])=="lfcSE"] <- "Cd630_v_Tox_both_lfcSE"


DC1 <- dplyr::select(transform[[1]], locus_tag, Cd630_v_Tox_DC1, Cd630_v_Tox_DC1_lfc, Cd630_v_Tox_DC1_lfcSE)
DC2 <- dplyr::select(transform[[2]], locus_tag, Cd630_v_Tox_DC2, Cd630_v_Tox_DC2_lfc, Cd630_v_Tox_DC2_lfcSE) 
both <- dplyr::select(transform[[3]], locus_tag, Cd630_v_Tox_both, Cd630_v_Tox_both_lfc, Cd630_v_Tox_both_lfcSE) 


mer <- merge(DC1, DC2, by = "locus_tag")
mer <- merge(mer, both, by = "locus_tag")
head(mer)


# read in genome annotation dataframe 
head(anno.630)


# if add gene name column, is there still complete lack of annotation information?
m <- match(mer$locus_tag, Cd12$Identifier)
anyNA(m) # FALSE
add <- Cd12[m, ]
mer$gene_ID <- add$ID.x
x <- mer # x is the dataframe to be annotated 

# add annotations and significance cutoff (adjusted p val 0.05) 
comb.anno <- function(x){
  mat <- match(x$locus_tag, anno.630$GeneID)
  add <- anno.630[mat, ]
  x$description <- add$Description
  x$protID <- add$Protein.product
  x$gene <- add$Locus
  x$locus_tag <- add$Locus.tag
  x$start <- add$Start
  x$Cd630_v_Tox_DC1_sig <- cut(x$Cd630_v_Tox_DC1, c(-100,-1.3,1.3,100), labels = c("down", "n.s.", "up")) # sig cut-off 0.05
  x$Cd630_v_Tox_DC2_sig <- cut(x$Cd630_v_Tox_DC2, c(-100,-1.3,1.3,100), labels = c("down", "n.s.", "up"))
  x$Cd630_v_Tox_both_sig <- cut(x$Cd630_v_Tox_both, c(-100,-1.3,1.3,100), labels = c("down", "n.s.", "up"))
  x
}

Cd630 <- comb.anno(x)

load("~/Box/PhD_papers/Ornithine/Final_submission_files/Source_Data/630_RNA-seq.RData")
