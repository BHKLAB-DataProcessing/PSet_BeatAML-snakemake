library(BiocGenerics)
library(readxl)
library(SummarizedExperiment)
library(abind)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[[1]], "download")
processed_dir <- paste0(args[[1]], "processed")

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_BeatAML-snakemake/bhklab_orcestra/snakemake/PSet_BeatAML/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_BeatAML-snakemake/bhklab_orcestra/snakemake/PSet_BeatAML/processed"

#######Creating eSet for  ALL consequences########
mutationData_all <-
  read_excel(file.path(download_dir, "Table-S7-Variants-for-Analysis.xlsx"), sheet = "Sheet 1")
#to change the labid to non-numeric
mutationData_all_labidchange <-
  paste("A", mutationData_all$labId, sep = "_")
mutationData_all$labId <- mutationData_all_labidchange
mutationData_subset_ALL <-
  as.data.frame(subset(
    x = mutationData_all,
    select = c("labId", "symbol", "ensembl_gene", "all_consequences")
  ))

#creating fdata
f_data_mutation_ALL <-
  unique(mutationData_subset_ALL[, c("symbol", "ensembl_gene"), drop = F])
rownames(f_data_mutation_ALL) <- f_data_mutation_ALL$symbol

#creating assaydata
#concatenate "all_consequences" to one using "///", ultimately mapping it to one gene
exprs_mutationSubset_ALL <-
  (mutationData_subset_ALL[, c("labId", "symbol", "all_consequences"), drop =
                             F])
uniqueSymbol <- unique(exprs_mutationSubset_ALL$symbol)
uniqueLabid <- unique(exprs_mutationSubset_ALL$labId)

#creating the matrix with NA values to fill the mutation data
exprs_mutation_ALL <-
  data.frame(matrix(
    NA,
    ncol = length(uniqueLabid),
    nrow = length(uniqueSymbol)
  ))
rownames(exprs_mutation_ALL) <- uniqueSymbol
colnames(exprs_mutation_ALL) <- uniqueLabid

#for loop to fill the matrix with concatenated values
for (i in 1:nrow(exprs_mutationSubset_ALL)) {
  symbol <- exprs_mutationSubset_ALL[i, "symbol"]
  labID <- exprs_mutationSubset_ALL[i, "labId"]
  consequence <- exprs_mutationSubset_ALL[i, "all_consequences"]
  existing_conseq <- exprs_mutation_ALL[symbol, labID]
  if (is.na(existing_conseq)) {
    exprs_mutation_ALL[symbol, labID] <- consequence
    
  }
  else {
    exprs_mutation_ALL[symbol, labID] <-
      paste(existing_conseq, "///", consequence, sep = "")
  }
}
#replacing NA with wt
###############################################################DO noT fill real NA with wt#############################
exprs_mutation_ALL = as.matrix(exprs_mutation_ALL)
exprs_mutation_ALL[is.na(exprs_mutation_ALL)] <- "wt"
#exprs_mutation_ALL=as.data.frame(exprs_mutation_ALL)
mutationData_ALLconseq <-
  Biobase::ExpressionSet(assayData = exprs_mutation_ALL)
fData(mutationData_ALLconseq) <- f_data_mutation_ALL
pdata_matrix_mutation_ALL <-
  matrix(NA, ncol = 2, nrow = length(unique(mutationData_subset_ALL$labId)))
colnames(pdata_matrix_mutation_ALL) <- c("batchid", "cellid")
pdata_matrix_mutation_ALL[, "cellid"] <-
  unique(mutationData_subset_ALL$labId)
rownames(pdata_matrix_mutation_ALL) <-
  unique(mutationData_subset_ALL$labId)
pData(mutationData_ALLconseq) <-
  as.data.frame(pdata_matrix_mutation_ALL)
stopifnot(all(unique(mutationData_subset_ALL$labId) == rownames(pData(
  mutationData_ALLconseq
))))

stopifnot(all(rownames(pdata_matrix_mutation_ALL) == rownames(pData(
  mutationData_ALLconseq
))))
stopifnot(all(rownames(f_data_mutation_ALL) == rownames(fData(
  mutationData_ALLconseq
))))
stopifnot(all(rownames(fData(
  mutationData_ALLconseq
)) == rownames(exprs(
  mutationData_ALLconseq
))))
stopifnot(all(rownames(pData(
  mutationData_ALLconseq
)) == colnames(exprs(
  mutationData_ALLconseq
))))

BiocGenerics::annotation(mutationData_ALLconseq) <- "mutation"
#convert eset to SE
new_SE_ALL <-
  SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(mutationData_ALLconseq)

stopifnot(all(rownames(colData(new_SE_ALL)) == rownames(pData(
  mutationData_ALLconseq
))))
stopifnot(all(rownames(rowData(new_SE_ALL)) == rownames(fData(
  mutationData_ALLconseq
))))

###########################################################################################################################################################################################################
#######Creating eSet for  Chosen consequences########

mutationData_subset_Chosen <-
  as.data.frame(subset(
    x = mutationData_all,
    select = c("labId", "symbol", "ensembl_gene", "chosen_consequence")
  ))

#creating fdata
f_data_mutation_Chosen <-
  unique(mutationData_subset_Chosen[, c("symbol", "ensembl_gene")])
rownames(f_data_mutation_Chosen) <- f_data_mutation_Chosen$symbol

#creating assaydata
#concatenate "chosen_consequence" to one using "///", ultimately mapping it to one gene

exprs_mutationSubset_Chosen <-
  (mutationData_subset_Chosen[, c("labId", "symbol", "chosen_consequence"), drop =
                                F])
uniqueSymbol_chosen <- unique(exprs_mutationSubset_Chosen$symbol)
uniqueLabid_chosen <- unique(exprs_mutationSubset_Chosen$labId)
#creating the matrix with NA values to fill the mutation data
exprs_mutation_Chosen <-
  data.frame(matrix(
    NA,
    ncol = length(uniqueLabid_chosen),
    nrow = length(uniqueSymbol_chosen)
  ))
rownames(exprs_mutation_Chosen) <- uniqueSymbol_chosen
colnames(exprs_mutation_Chosen) <- uniqueLabid_chosen

#for loop to fill the matrix with concatenated values
for (i in 1:nrow(exprs_mutationSubset_Chosen)) {
  symbol <- exprs_mutationSubset_Chosen[i, "symbol"]
  labID <- exprs_mutationSubset_Chosen[i, "labId"]
  consequence <-
    exprs_mutationSubset_Chosen[i, "chosen_consequence"]
  existing_conseq <- exprs_mutation_Chosen[symbol, labID]
  if (is.na(existing_conseq)) {
    exprs_mutation_Chosen[symbol, labID] <- consequence
    
  }
  else {
    exprs_mutation_Chosen[symbol, labID] <-
      paste(existing_conseq, "///", consequence, sep = "")
  }
}

#replacing NA with wt
exprs_mutation_Chosen = as.matrix(exprs_mutation_Chosen)
exprs_mutation_Chosen[is.na(exprs_mutation_Chosen)] <- "wt"
#exprs_mutation_Chosen=as.data.frame(exprs_mutation_Chosen)
mutationData_ChosenConseq <-
  Biobase::ExpressionSet(assayData = (exprs_mutation_Chosen))
fData(mutationData_ChosenConseq) <- f_data_mutation_Chosen
pdata_matrix_mutation_Chosen <-
  matrix(NA, ncol = 2, nrow = length(unique(mutationData_subset_Chosen$labId)))
colnames(pdata_matrix_mutation_Chosen) <- c("batchid", "cellid")
pdata_matrix_mutation_Chosen[, "cellid"] <-
  unique(mutationData_subset_Chosen$labId)
rownames(pdata_matrix_mutation_Chosen) <-
  unique(mutationData_subset_Chosen$labId)
pData(mutationData_ChosenConseq) <-
  as.data.frame(pdata_matrix_mutation_Chosen)
stopifnot(all(unique(mutationData_ChosenConseq$labId) == rownames(pData(
  mutationData_ChosenConseq
))))

stopifnot(all(rownames(pdata_matrix_mutation_Chosen) == rownames(pData(
  mutationData_ChosenConseq
))))
stopifnot(all(rownames(f_data_mutation_Chosen) == rownames(fData(
  mutationData_ChosenConseq
))))
stopifnot(all(rownames(fData(
  mutationData_ChosenConseq
)) == rownames(exprs(
  mutationData_ChosenConseq
))))
stopifnot(all(rownames(pData(
  mutationData_ChosenConseq
)) == colnames(exprs(
  mutationData_ChosenConseq
))))

BiocGenerics::annotation(mutationData_ChosenConseq) <- "mutation"
#convert eset to SE
new_SE_Chosen <-
  SummarizedExperiment::makeSummarizedExperimentFromExpressionSet(mutationData_ChosenConseq)

stopifnot(all(rownames(colData(new_SE_Chosen)) == rownames(pData(
  mutationData_ChosenConseq
))))
stopifnot(all(rownames(rowData(new_SE_Chosen)) == rownames(fData(
  mutationData_ChosenConseq
))))

saveRDS(new_SE_ALL, file.path(processed_dir, 'new_SE_ALL.rds'))
saveRDS(new_SE_Chosen, file.path(processed_dir, 'new_SE_Chosen.rds'))
