library(BiocGenerics)
library(readxl)
library(SummarizedExperiment)
library(abind)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[[1]], "download")
processed_dir <- paste0(args[[1]], "processed")

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_BeatAML-snakemake/bhklab_orcestra/snakemake/PSet_BeatAML/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_BeatAML-snakemake/bhklab_orcestra/snakemake/PSet_BeatAM/processed"
 

###########################################################################################################################################################################################################
#######Creating eSet for  rnaseq########

#kallisto updated in Aug 2020 by Anthony
#code - https://github.com/sisiranair/BeatAML/blob/master/R/createPset/rnaseq_kallisto.R
#kallisto rerun in Mar 2021 by Anthony
#load("~/Downloads/beatAML_SE.RData") # backed up on H4H and One drive
# BeatAML_gene_exp_20210304 <- rnaseq_se$rnaseq
# saveRDS(BeatAML_gene_exp_20210304, "data/BeatAML_gene_exp_20210304.rds")
# Updated to include all rnaseq data types on Feb 8, 2022
kallisto_rnaseq_all <- readRDS(file.path(download_dir, "beatAML_SE.rds"))


fixSEAnnotations <- function(kallisto_rnaseqAML){
  
  
  rownames(colData(kallisto_rnaseqAML)) <-
    gsub("Sample_", "A_", rownames(colData(kallisto_rnaseqAML)))
  
  #in id col of pdata
  kallisto_rnaseqAML$sample <-
    gsub("Sample_", "A_", kallisto_rnaseqAML$sample)
  
  #renaming id to cellid as cellid should be present in pdata
  colnames(colData(kallisto_rnaseqAML))[colnames(colData(kallisto_rnaseqAML)) == 'sample'] <-
    'cellid'
  
  #to unlock the object/environment
  #storageMode(kallisto_rnaseqAML) <- "environment"
  colnames(assay(kallisto_rnaseqAML)) <-
    gsub("Sample_", "A_", colnames(assay(kallisto_rnaseqAML)))
  
  #to lock the object/environment
  #storageMode(kallisto_rnaseqAML) <- "lockedEnvironment"
  
  #addind batchid
  #batchid <- NA
  #pData(kallisto_rnaseqAML)<- cbind(pData(kallisto_rnaseqAML),batchid)
  
  #changing colname -gene_name to Symbol
  #colnames(fData(kallisto_rnaseqAML))[12] <- "Symbol"
  colnames(rowData(kallisto_rnaseqAML))[colnames(rowData(kallisto_rnaseqAML)) == 'gene_name'] <-
    'Symbol'
  
  #adding BEST
  rowData(kallisto_rnaseqAML)$BEST <- NA
  
  stopifnot(all(rownames(rowData(kallisto_rnaseqAML)) == rownames(assay(kallisto_rnaseqAML))))
  stopifnot(all(rownames(colData(kallisto_rnaseqAML)) == colnames(assay(kallisto_rnaseqAML))))
  
  return(kallisto_rnaseqAML)
  
}

kallisto_rnaseq_all <- lapply(kallisto_rnaseq_all, fixSEAnnotations)



kallisto_rnaseqAML <- kallisto_rnaseq_all$rnaseq

#convert eset to SE
#new_SE_rnaseq <-SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(kallisto_rnaseqAML)

#stopifnot(all(rownames(colData(new_SE_rnaseq)) == rownames(pData(kallisto_rnaseqAML))))
#stopifnot(all(rownames(rowData(new_SE_rnaseq)) == rownames(fData(kallisto_rnaseqAML))))
###########################################################################################################################################################################################################

########Creating cell object########

mutationData_clinicalSummary <-
  read_excel(file.path(download_dir, "Table-S5-Clinical-Summary.xlsx"), sheet = "Sheet 1")
mutationData_clinicalSummary$LabId <-
  paste("A", mutationData_clinicalSummary$LabId, sep = "_")
cell_AML <- as.data.frame(unique(mutationData_clinicalSummary))
colnames(cell_AML)[colnames(cell_AML) == 'LabId'] <-
  'cellid' #' @check this

#adding the extra cellids from drug experiments
read_drugResponse <-
  read_excel(file.path(download_dir, "Table-S10-Drug-Responses.xlsx"), sheet = "Sheet 1")
read_drugResponse.1 <- read_drugResponse
read_drugResponse.1$lab_id <-
  paste("A", read_drugResponse.1$lab_id, sep = "_")
#to check which labids are absent in  drug data
celllines118 <- setdiff(read_drugResponse.1$lab_id, cell_AML$cellid)
#for i in cellines118, adds a new row wth NA
for (i in 1:length(celllines118)) {
  cell_AML[nrow(cell_AML) + 1, ] = c(celllines118[i], rep(NA, length(colnames(cell_AML)) -
                                                            1))
}

#pull data from kallisto
#adding cellid only present in kallisto rnaseq
diffCellids <- setdiff(kallisto_rnaseqAML$cellid , cell_AML$cellid)

#for i in diffCellids, adds a new row wth NA
for (i in 1:length(diffCellids)) {
  cell_AML[nrow(cell_AML) + 1, ] = c(diffCellids[i], rep(NA, length(colnames(cell_AML)) -
                                                           1))
}
rownames(cell_AML) <- cell_AML$cellid

#########################################################################################################################################################################

#######CREATING CURATIONCELL########
curationCell <-
  data.frame(cell_AML$cellid, row.names = cell_AML$cellid)
curationCell$unique.cellid <- curationCell$cell_AML.cellid

#######CREATING CURATIONDRUG########
drugs_with_ids <-
  read.csv(file.path(download_dir, "drugs_with_ids.csv"), stringsAsFactors = F)
curationDrug <-
  subset(
    drugs_with_ids,
    subset = !is.na(drugs_with_ids$BeatAML.drugid) ,
    select = c("unique.drugid", "BeatAML.drugid")
  )
colnames(curationDrug) <- c("unique.drugid", "AML.drugid")

# remove duplicate drug ids
ind <- duplicated(curationDrug$unique.drugid)
curationDrug <- curationDrug[!ind, ]
rownames(curationDrug) <- curationDrug[, "unique.drugid"]

###########################################################################################################################################################################################################
########Creating drug object########

drug_AML <- as.data.frame(unique(read_drugResponse$inhibitor))
colnames(drug_AML) <- "drug.name"

merge_drugAml <-
  merge(drug_AML, curationDrug, by.x = "drug.name", by.y = "AML.drugid")
drug_AML <- merge_drugAml

rownames(drug_AML) <- drug_AML$unique.drugid
drug_AML$drugid <- paste("drugid", drug_AML$drug.name, sep = "_")

drug_AML <- drug_AML[, -grep("unique.drugid", colnames(drug_AML))]

saveRDS(cell_AML, file.path(processed_dir, 'cell_AML.rds'))
saveRDS(drug_AML, file.path(processed_dir, 'drug_AML.rds'))
saveRDS(kallisto_rnaseq_all, file.path(processed_dir, 'kallisto_rnaseq_all.rds'))
saveRDS(curationDrug, file.path(processed_dir, 'curationDrug.rds'))
saveRDS(curationCell, file.path(processed_dir, 'curationCell.rds'))

