# Daniel J Fulop
# TCGA & ICGC Immune Cell Data derived from immune deconvolution methods
# in immunedeconv package (https://icbi-lab.github.io/immunedeconv/)
# to infer cell-type proportions from bulk RNA-seq data from TCGA tumor samples.
# Functions to clean tcga and icgc data into required format,
# run deconvolution methods

# load dependent libraries
library(tidyverse)
library(ggpubr)
library(reshape2)
library(knitr)
library(immunedeconv) # recommend using conda to install due to issues with dependent package versions
library(annotables)
library(WGCNA)
library(XLConnect)
library(sqldf)

# load in gzip tsvs or tsv filetypes
# increase buffering size to load large data size
# returns list of data frames
setData <- function(filePaths) {
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 4)
  l <- list()
  for (f in 1:length(filePaths)) {
    # gzip file containing tsv
    if (substr(filePaths[f],nchar(filePaths[f])-1,nchar(filePaths[f]))=="gz" &
        substr(filePaths[f],nchar(filePaths[f])-5,nchar(filePaths[f])-3)=="tsv") {
      l[[f]] <- read_tsv(gzfile(filePaths[f]))
      }
    if (substr(filePaths[f],nchar(filePaths[f])-2,nchar(filePaths[f]))=="tsv") {
      l[[f]] <- read_tsv(filePaths[f])
    }
  }
  return(l)
}

# immunedeconv::deconvolute function accepts matrix/data frame of
# HGNC symbols as row names and samples as columns
tcgaFormat <- function(tcgaDat) {
  # clean HGNC symbols
  tcgaDat$gene_id <- gsub("\\|.*","",tcgaDat$gene_id)
  # remove hypothetical loci/unnamed transcripts coded with "?"
  tcgaDat <- subset(tcgaDat, tcgaDat$gene_id != "?")
  # remove duplicate "SLC35E2" row with smaller values per "https://rdrr.io/github/jefferys/FusionExpressionPlot/man/extractGeneModels.html"
  tcgaDat <- tcgaDat %>% 
    group_by(gene_id) %>% 
    filter(`TCGA-OR-A5J1-01A-11R-A29S-07` == max(`TCGA-OR-A5J1-01A-11R-A29S-07`)) %>% 
    distinct
  # convert gene_id variable to row name
  tcgaDat <- column_to_rownames(tcgaDat, var = 'gene_id')
  return(tcgaDat)
}

# Function below modifies cleaned TCGA data and runs Timer outputting immune cell fractions
# timer method requires indications vector of cancer types per sample
# necessitating a modification to the already cleaned TCGA data
# also requires data file containing mapping of samples to cancer type
# in file barcodeCancerMapping, bcr_patient_barcode = TCGA_barcode, acronym = HGNC codes
tcgaTimerDeconv <- function(tcgaClean,barcodeCancerMapping) {
  barcodeCancerMapping <- barcodeCancerMapping %>%
    select("bcr_patient_barcode", "acronym")
  # create df of TCGA_barcodes in order
  tcgaSampleBarcodes <- tibble(colnames(tcgaClean))
  colnames(tcgaSampleBarcodes) <- c("TCGA_BARCODE")
  tcgaSampleBarcodes$bcr_patient_barcode <- str_trunc(tcgaSampleBarcodes$TCGA_BARCODE, 12, "right", ellipsis = "")
  # join cancer map and barcode table, 10856 samples have cancer designation, 214 do not 
  tcgaBarcodeCancernoNAs <- left_join(tcgaSampleBarcodes, barcodeCancerMapping, by = "bcr_patient_barcode") %>%
    filter(!is.na(acronym))
  # transpose TCGA data for join
  tcgaTransposed <- as.data.frame(t(tcgaClean))
  tcgaTransposed <- rownames_to_column(tcgaTransposed, var = "TCGA_BARCODE")
  # INNER join with TCGA_barcode_cancer_noNAs to remove samples without cancer designation
  tcgaTransposedNoNAsCancer <- inner_join(tcgaTransposed, tcgaBarcodeCancernoNAs, by = "TCGA_BARCODE")
  # remove truncated barcode column
  tcgaTransposedNoNAsCancer <- tcgaTransposedNoNAsCancer[ , !(names(tcgaTransposedNoNAsCancer) %in% c("bcr_patient_barcode"))]
  tcgaTransposedNoNAsCancer <- column_to_rownames(tcgaTransposedNoNAsCancer, var = "TCGA_BARCODE")
  # transpose data to sample barcodes as variables (format to be used by deconvolution methods)
  tcgaNoNAsCancer <- as.data.frame(t(tcgaTransposedNoNAsCancer))
  # create vector of indications (HGNC cancer symbols)
  cancerIndications <- as.character(tcgaNoNAsCancer[20502,])
  # drop acronym row from dataset
  tcgaNoNAsCancer <- tcgaNoNAsCancer[-20502,]
  # run timer
  tcgaTimerFrac <- deconvolute(tcgaNoNAsCancer,"timer",indications=cancerIndications)
  return(tcgaTimerFrac)
}

# get Cibersort TCGA data into same format deconvolution output
# run and published by Thorsson et al.
tcgaCibersortFormat <- function(tcgaCibersort) {
  tcgaCibersort$SampleID <- str_replace_all(tcgaCibersort$SampleID, "\\.", "-")
  tcgaCibersort <- select(tcgaCibersort, -c("CancerType"))
  cibersortTransposed <- as.data.frame(t(tcgaCibersort))
  names(cibersortTransposed) <- as.character(unlist(cibersortTransposed[1,]))
  cibersortTransposed <- cibersortTransposed[-1,]
  cibersortTransposed <- rownames_to_column(cibersortTransposed, var = "cell_type")
  cibersortTransposed$cell_type <- str_replace_all(cibersortTransposed$cell_type, "\\.", " ")
  cibersortTransposed <- cibersortTransposed[!cibersortTransposed$cell_type %in% c("P value","Correlation","RMSE"),]
  newNames <- c("B cell naive", "B cell memory","B cell plasma","T cell CD8+","T cell CD4+ naive",
                +               "T cell CD4+ memory resting","T cell CD4+ memory activated","T cell follicular helper",
                +               "T cell regulatory (Tregs)","T cell gamma delta","NK cell resting","NK cell activated",
                +               "Monocyte","Macrophage M0","Macrophage M1","Macrophage M2","Myeloid dendritic cell resting",
                +               "Myeloid dendritic cell activated","Mast cell resting","Mast cell activated","Eosinophil",
                +               "Neutrophil")
  cibersortTransposed$cell_type <- newNames
  return(cibersortTransposed)
}

# Format ICGC data for use by immunedeconv::deconvolute
icgcFormat <- function(icgcDat) {
  # load ensemble to HGNC mapping, eliminating duplicate rows
  HGNC_ensemble_map <- grch38 %>% 
    select(ensgene, symbol) %>%
    distinct()
  icgcDat$ensgene <- str_trunc(icgcDat$feature, 15, "right", ellipsis = "")
  # join ICGC with HGNC map on contains "ensgene"
  icgcDatJoin <- left_join(icgcDat, HGNC_ensemble_map, by = "ensgene") %>%
    select(-c("ensgene", "feature"))
  # 5249 nas (ie. rows without HGNC symbol designation)
  # remove observations where symbol is na
  icgcDatNoNA <- icgcDatJoin[!is.na(icgcDatJoin$symbol),]
  rowID <- rownames(icgcDatNoNA)
  rowGroup <- icgcDatNoNA$symbol
  # remove symbol
  icgcDatNoNA <- select(icgcDatNoNA, -c("symbol"))
  # 142 HGNC symbols with multiple observations (1359 total) in ICGC 
  # collapse symbol to one observation/HGNC symbol by selecting the row with the fewest missing values (select MaxMean)
  icgcDatCollapsed <- collapseRows(icgcDatNoNA, rowGroup=rowGroup, rowID=rowID, selectFewestMissing = TRUE)
  # ICGC_collapsed to df with HGNC symbols as rownames and only one observation/symbol
  icgcDf <- data.frame(icgcDatCollapsed$group2row, icgcDatCollapsed$datETcollapsed)
  # remove unnecessary columns
  icgcDf <- select(icgcDf, -c("group", "selectedRowID"))
  return(icgcDf)
}

icgcTimerDeconv <- function(icgcClean,cancerSubtypes) {
  # get cancer type from dcc_project_code
  cancerSubtypes$cancer_type <- gsub("-.*","",cancerSubtypes$dcc_project_code)
  # aliquot_ID contains key for joining ICGC_collapsed_df)
  # get rid of unnecessary columns -- only keep aliquot_id and cancer_type
  cancerSubtypes <- cancerSubtypes %>%
    select(aliquot_id, cancer_type)
  # Transpose ICGC data to join cancer subtypes
  icgcTransposed <- as.data.frame(t(icgcClean))
  icgcTransposed <- rownames_to_column(icgcTransposed, var = "aliquot_id")
  # Transform aliquot_id: 1) remove "X" from string; 2) replace "." with "-"
  icgcTransposed$aliquot_id <- sub('.', '', icgcTransposed$aliquot_id)
  icgcTransposed$aliquot_id <- gsub('\\.', '-', icgcTransposed$aliquot_id)
  # Inner join ICGC_transposed and cancer_subtypes by aliquot_id
  # 980 observations
  icgcCancertype <- inner_join(icgcTransposed, cancerSubtypes, by = "aliquot_id")
  # get only cancer types compatible with timer -- 620 samples
  # capitalize timer_available_cancers
  upper_timer_available_cancers <- toupper(timer_available_cancers)
  timer_cancers_only <- as.data.frame(icgcCancertype[icgcCancertype$cancer_type %in% upper_timer_available_cancers, ])
  # reset dataframe rownames to allow for aliquot_id to become rownames for transposition
  rownames(timer_cancers_only) <- NULL
  # Column to rownames
  timer_cancers_only <- column_to_rownames(timer_cancers_only, var = "aliquot_id")
  # transpose data to aliquot_id as variables (format to be used by deconvolution methods)
  icgcCancerTransposed <- as.data.frame(t(timer_cancers_only))
  # create vector of indications (HGNC cancer symbols)
  cancerIndications <- as.character(icgcCancerTransposed[51355,])
  # drop HGNC symbol row from dataset
  icgcCancerTransposed <- icgcCancerTransposed[-51355,]
  # run deconvolution
  icgcTimerFrac <- deconvolute(icgcCancerTransposed,"timer",indications=cancerIndications)
  return(icgcTimerFrac)
}

# run deconvolution methods and output single data frame with immune fractions 
# for all input samples generated by passed deconvolution method
getFractions <- function(dataList,method,celltypes=NULL) {
  listofdfs <- list()
  for (i in 1:length(dataList)) {
      listofdfs[[i]] <- deconvolute(dataList[[i]],method,tumor=TRUE,expected_cell_types=celltypes)
      # drop cell_type column from subsequent deconvolutions for easy df binding
      if (i>1) {
        listofdfs[[i]] <- listofdfs[[i]][,-1]
      }
  }
  # generate one df containing immune fractions for all samples
  immuneFrac <- do.call(cbind,listofdfs)
  return(immuneFrac)
}

# map immune fractions to common cell vocabulary shared by all methods using immunedeconv::map_result_to_celltypes
# requires user generated vector of cell types from immunedeconv::cell_type_list hierarchically maps to common cell types
# accepts list of data frames
commonCells <- function(listofdfs,desiredCells) {
  dfs <- list()
  for (i in 1:length(listofdfs)) {
    dfs[[i]] <- map_result_to_celltypes(listofdfs[[i]],desiredCells)
  }
  return(dfs)
}

# load data
files <- c("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv.gz","tophat_star_fpkm_uq.v2_aliquot_gl.tsv.gz","TCGA.Kallisto.fullIDs.cibersort.relative.tsv")
dfs <- setData(files)
tcgaDat <- dfs[[1]]
icgcDat <- dfs[[2]]
tcgaCibersortUnformatted <- dfs[[3]] # From Thorsson et al. paper

# clean data
tcgaClean <- tcgaFormat(tcgaDat)
icgcClean <- icgcFormat(icgcDat)

# Subset TCGA data to create managable dfs able to be processed by deconvolute function
tcga1 <- tcgaClean[,1:2767]
tcga2 <- tcgaClean[,2768:4152]
tcga3 <- tcgaClean[,4153:5534] 
tcga4 <- tcgaClean[,5535:8301]
tcga5 <- tcgaClean[,8302:9000] 
tcga6 <- tcgaClean[,9001:9685] 
tcga7 <- tcgaClean[,9686:10300]
tcga8 <- tcgaClean[,10301:11069]



### RUN DECONVOLUTIONS

# Different subsets of tcga data used for different deconvolution methods
# due to errors resulting from missingness in samples handled by different
# deconvolution methods. To be used as input to getFractions.
# Quantiseq should be run separately from the other methods, which can be run together.
tcgaQuantiseqDfs <- list(tcga1,tcga2,tcga3,tcga4,tcga5)
# Data as input to getFractions for mcp, epic, and xcell deconvolution methods 
tcgaMcpEpicXcellDfs <- list(tcga1,tcga2,tcga3,tcga4,tcga5,tcga6,tcga7,tcga8)
# icgc deconvolutions can be run without subsetting dataset
icgcDatList <- list()
icgcDatList[[1]] <- icgcClean

# RUN TCGA DECONVOLUTIONS -- Quantiseq, MCP Counter, Epic, Xcell, Timer
# Note that this is a computationally intensive process and may take a long time
# to run on all samples. Consider subsetting or only running on samples of desired cancer type.
tcgaQuantiseq <- getFractions(tcgaQuantiseqDfs,"quantiseq",celltypes=NULL)
tcgaMcp <- getFractions(tcgaMcpEpicXcellDfs,"mcp",celltypes=NULL)
tcgaEpic <- getFractions(tcgaMcpEpicXcellDfs,"epic",celltypes=NULL)
# When running xcell, cell type vector should be used to avoid
# overcompensation by spillover correction as described by Sturm et al.
# set celltypes to cellTypeVect
cellTypeVect = c("B cell", "Myeloid dendritic cell", "T cell CD4+", "T cell CD8+", "Macrophage/Monocyte", "NK cell", "T cell CD4+ (non-regulatory)", "T cell regulatory (Tregs)", "Cancer associated fibroblast", "Endothelial cell")
tcgaXcell <- getFractions(tcgaMcpEpicXcellDfs,"xcell",celltypes=cellTypeVect)
# sample-cancer type mapping required for timer
# use in tcgaTimerDeconv
barcodeCancerMapping <- read_tsv("clinical_PANCAN_patient_with_followup.tsv")
tcgaTimer <- tcgaTimerDeconv(tcgaClean,barcodeCancerMapping)
tcgaCibersort <- tcgaCibersortFormat(tcgaCibersortUnformatted)

# RUN ICGC DECONVOLUTIONS -- Quantiseq, MCP Counter, Epic, Xcell, Timer
icgcQuantiseq <- getFractions(icgcDatList,"quantiseq",celltypes=NULL)
icgcMcp <- getFractions(icgcDatList, "mcp_counter",celltypes=NULL)
icgcEpic <- getFractions(icgcDatList, "epic",celltypes=NULL)
icgcXcell <- getFractions(icgcDatList, "xcell",celltypes=cellTypeVect)
# timer -- requires file containing cancer subtypes
cancerSubtypes <- read_tsv("pcawg_sample_sheet.tsv")
icgcTimer <- icgcTimerDeconv(icgcClean,cancerSubtypes)



### COMMON CELL VOCABULARY

# Based on Sturm et al. benchmarking cell-type quantifications, the cells should be:
# B, DC, Mac/Mono, NK, T CD4+, T CD8+, T CD4+ nr, T reg, CAF, Endo
desiredCells <- c("B cell", "Myeloid dendritic cell", "T cell CD4+", "T cell CD8+", "Macrophage/Monocyte",
                  "NK cell", "T cell CD4+ (non-regulatory)","T cell regulatory (Tregs)","Cancer associated fibroblast",
                  "Endothelial cell")

#tcga
tcgaListOfDfs <- list(tcgaQuantiseq,tcgaMcp,tcgaEpic,tcgaXcell,tcgaTimer,tcgaCibersort)
tcgaCv <- commonCells(tcgaListOfDfs,desiredCells)
# data frame for each method
tcgaQuantiseqMapped <- tcgaCv[[1]]
tcgaMcpMapped <- tcgaCv[[2]][-c(3,7:8),]
tcgaEpicMapped <- tcgaCv[[3]]
tcgaXcellMapped <- tcgaCv[[4]]
tcgaTimerMapped <- tcgaCv[[5]]
tcgaCibersortMapped <- tcgaCv[[6]]


# icgc 
icgcListOfDfs <- list(icgcQuantiseq,icgcMcp,icgcEpic,icgcXcell,icgcTimer)
icgcCv <- commonCells(icgcListOfDfs,desiredCells)
# dataframe for each method
icgcQuantiseqMapped <- icgcCv[[1]]
icgcMcpMapped <- icgcCv[[2]]
icgcEpicMapped <- icgcCv[[3]]
icgcXcellMapped <- icgcCv[[4]]
icgcTimerMapped <- icgcCv[[5]]



  