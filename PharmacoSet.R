library(PharmacoGx)
library(Biobase)
library(BiocGenerics)
library(readxl)
library(SummarizedExperiment)
library(abind)

#######Creating eSet for  ALL consequences########
mutationData_all <-
    read_excel("/pfs/beatAML_raw/Table-S7-Variants-for-Analysis.xlsx", sheet = "Sheet 1")
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
###########################################################################################################################################################################################################
#######Creating eSet for  rnaseq########

#kallisto updated in Aug 2020 by Anthony
#code - https://github.com/sisiranair/BeatAML/blob/master/R/createPset/rnaseq_kallisto.R
#kallisto rerun in Mar 2021 by Anthony
#load("~/Downloads/beatAML_SE.RData") # backed up on H4H and One drive
# BeatAML_gene_exp_20210304 <- rnaseq_se$rnaseq
# saveRDS(BeatAML_gene_exp_20210304, "data/BeatAML_gene_exp_20210304.rds")
# Updated to include all rnaseq data types on Feb 8, 2022
kallisto_rnaseq_all <- readRDS("/pfs/beatAML_raw/beatAML_SE.rds")


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
    read_excel("/pfs/beatAML_raw/Table-S5-Clinical-Summary.xlsx", sheet = "Sheet 1")
mutationData_clinicalSummary$LabId <-
    paste("A", mutationData_clinicalSummary$LabId, sep = "_")
cell_AML <- as.data.frame(unique(mutationData_clinicalSummary))
colnames(cell_AML)[colnames(cell_AML) == 'LabId'] <-
    'cellid' #' @check this

#adding the extra cellids from drug experiments
read_drugResponse <-
    read_excel("/pfs/beatAML_raw/Table-S10-Drug-Responses.xlsx", sheet = "Sheet 1")
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
    read.csv("/pfs/beatAML_raw/drugs_with_ids.csv", stringsAsFactors = F)
curationDrug <-
    subset(
        drugs_with_ids,
        subset = !is.na(drugs_with_ids$BeatAML.drugid) ,
        select = c("unique.drugid", "BeatAML.drugid")
    )
colnames(curationDrug) <- c("unique.drugid", "AML.drugid")
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

#########################################################################################################################################################################
########Creating drug object########
#raw_data_pointsAML <- read.delim("data/inhibitor_data_points_2018_10_24.txt",sep="\t", header=TRUE )#old
#csv with conc values that needs to be changed
#concentrationMapping <- read.csv(("data/concentration_mapping.csv"), header = TRUE,fileEncoding="UTF-8-BOM")#old

#PROCESSING RAW DOSE RESPONSE - updated data from Jeffrey Tyner OHSU
raw_data_pointsAML <-
    read.table(
        "/pfs/beatAML_raw/beataml_manuscript_raw_inhib_v2_passed_qc_ceil_100.tsv",
        sep = "\t",
        header = TRUE
    )
#most recent panel version of cleaned up concentrations with consistent rounding from Jeffrey Tyner
clean_conc <-
    read.csv("/pfs/beatAML_raw/report_1614271471952.csv", stringsAsFactors = F)#only for reference

#equating dose values to ref values sent by Jeff
#checks to see if merging values would create duplicated rows
check_val <- function(raw.df, val1, val2) {
    rr <-
        raw.df[raw.df$well_concentration == val1 |
                   raw.df$well_concentration == val2, ]
    if (length(unique(rr$well_concentration)) != "2") {
        print("Not all input values are present in subset")
    } else{
        rr_uni <- unique(rr$inhibitor)
        for (dr in 1:length(rr_uni)) {
            rr_dr <- rr[rr$inhibitor == rr_uni[dr], ]
            if (length(unique(rr_dr$well_concentration)) > 1) {
                print(
                    paste(
                        rr_uni[dr],
                        "has both values. This might create duplicates",
                        sep = " "
                    )
                )
            } else {
                print(paste("Pass - no duplicates found for", rr_uni[dr], sep = " "))
            }
        }
    }
}

#for values 0.001372 &  0.001400
val.pair.1 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "0.001372",
              val2 = "0.0014")
#for values 0.004100 &  0.004115
val.pair.2 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "0.0041",
              val2 = "0.004115")
#for values 0.012300 & 0.012346
val.pair.3 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "0.0123",
              val2 = "0.012346")
#for values 0.013700 & 0.013717
val.pair.4 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "0.0137",
              val2 = "0.013717")
#for values 0.037000 & 0.037037
val.pair.5 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "0.037",
              val2 = "0.037037")
#for values 0.041152 & 0.041200
val.pair.6 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "0.041152",
              val2 = "0.0412")
#for values 0.111100 & 0.111111
val.pair.7 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "0.1111",
              val2 = "0.111111")
#for values 0.123457 & 0.123500
val.pair.8 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "0.123457",
              val2 = "0.1235")
#for values 0.333300 & 0.333333
val.pair.9 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "0.3333",
              val2 = "0.333333")
#for values 0.370370 & 0.370400
val.pair.10 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "0.37037",
              val2 = "0.3704")
#for values 1.111100 & 1.111110
val.pair.11 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "1.1111",
              val2 = "1.11111")
#for values 3.333300 & 3.333330
val.pair.12 <-
    check_val(raw.df = raw_data_pointsAML,
              val1 = "3.3333",
              val2 = "3.33333")


updatedRaw_data_pointAML <- raw_data_pointsAML

#new conversions
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.001372] <-
    0.00137
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.0014] <-
    0.00137

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.0041] <-
    0.00412
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.004115] <-
    0.00412

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.006859] <-
    0.00686

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.0123] <-
    0.0123
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.012346] <-
    0.0123

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.0137] <-
    0.0137
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.013717] <-
    0.0137

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.020576] <-
    0.0206

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.037] <-
    0.037
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.037037] <-
    0.037

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.041152] <-
    0.0412
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.0412] <-
    0.0412

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.061728] <-
    0.0617

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.1111] <-
    0.111
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.111111] <-
    0.111

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.123457] <-
    0.123
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.1235] <-
    0.123

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.185185] <-
    0.185

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.3333] <-
    0.333
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.333333] <-
    0.333

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.37037] <-
    0.37
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.3704] <-
    0.37

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 0.555556] <-
    0.556

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 1.1111] <-
    1.11
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 1.11111] <-
    1.11

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 1.66667] <-
    1.67

updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 3.3333] <-
    3.33
updatedRaw_data_pointAML$well_concentration[updatedRaw_data_pointAML$well_concentration == 3.33333] <-
    3.33
# Doses "1.000000", "5.000000", "10.000000" are the same, hence not changed
##################################################################################################################################################################
###########CREATING 3D ARRAY#############

conc_testedAML <- 21

new_ids_AML <- updatedRaw_data_pointAML
new_ids_AML$lab_id <- paste("A", new_ids_AML$lab_id, sep = "_")
new_ids_AML$merged_id <-
    paste(
        "drugid",
        new_ids_AML$inhibitor,
        new_ids_AML$lab_id,
        new_ids_AML$inhibitor_panel,
        new_ids_AML$run_index,
        new_ids_AML$replicate,
        new_ids_AML$time_of_read ,
        sep = "."
    )

unique_exp_idsAML <- unique(new_ids_AML$merged_id)

#Creating dose Array
doseArrayAML <-
    array(NA, dim = c(length(unique_exp_idsAML), conc_testedAML))
rownames(doseArrayAML) <- unique_exp_idsAML
colnames(doseArrayAML) <-
    c(
        "doses1",
        "doses2",
        "doses3",
        "doses4",
        "doses5",
        "doses6",
        "doses7",
        "doses8",
        "doses9",
        "doses10",
        "doses11",
        "doses12",
        "doses13",
        "doses14",
        "doses15",
        "doses16",
        "doses17",
        "doses18",
        "doses19",
        "doses20",
        "doses21"
    )

#Creating viability Array
viabilityArrayAML <-
    array(NA, dim = c(length(unique_exp_idsAML), conc_testedAML))
rownames(viabilityArrayAML) <- unique_exp_idsAML
colnames(viabilityArrayAML) <-
    c(
        "doses1",
        "doses2",
        "doses3",
        "doses4",
        "doses5",
        "doses6",
        "doses7",
        "doses8",
        "doses9",
        "doses10",
        "doses11",
        "doses12",
        "doses13",
        "doses14",
        "doses15",
        "doses16",
        "doses17",
        "doses18",
        "doses19",
        "doses20",
        "doses21"
    )

dosea <- sort(unique(new_ids_AML$well_concentration))
names(dosea) <- colnames(doseArrayAML)


for (row in 1:nrow(doseArrayAML)) {
    rid <- rownames(doseArrayAML)[row]
    tdf <- new_ids_AML[new_ids_AML$merged_id == rid,]
    for (i in 1:nrow(tdf))
    {
        doname <- names(dosea[dosea == tdf[i, "well_concentration"]])
        doseArrayAML[rid, doname] <- tdf[i, "well_concentration"]
        viabilityArrayAML[rid, doname] <-
            tdf[i, "normalized_viability_ceil_100"]
    }
    
}

doseViabilityArray_AML <-
    abind(doseArrayAML, viabilityArrayAML, along = 3)
dimnames(doseViabilityArray_AML)[[3]] <- c("Dose", "Viability")
#saveRDS(doseViabilityArray_AML, "data/doseViabilityArray_AML.rds")
#load array
#doseViabilityArray_AML <- readRDS("data/doseViabilityArray_AML.rds")
###########CREATING INFO#############
info_new_ids_AML_subset <-
    new_ids_AML[, c(
        "lab_id",
        "inhibitor",
        "inhibitor_panel",
        "run_index",
        "replicate",
        "time_of_read",
        "merged_id"
    )]
info_new_ids_AML_subset <- unique(info_new_ids_AML_subset)

merge_infoAml_drugid <-
    merge(info_new_ids_AML_subset,
          curationDrug,
          by.x = "inhibitor",
          by.y = "AML.drugid")
#merge_infoAml_drugid <- merge_infoAml_drugid[,-1]
#colnames(merge_infoAml_drugid)[3] <- "drugid"
names(merge_infoAml_drugid)[names(merge_infoAml_drugid) == 'unique.drugid'] <-
    'drugid'

merge_infoAml_drugid[, "nbr.conc.tested"] <- 21
#colnames(merge_infoAml_drugid)[2] <- "cellid"
names(merge_infoAml_drugid)[names(merge_infoAml_drugid) == 'lab_id'] <-
    'cellid'

rownames(merge_infoAml_drugid) <- merge_infoAml_drugid$merged_id
merge_infoAml_drugid <-
    merge_infoAml_drugid[, -grep("merged_id", colnames(merge_infoAml_drugid))]
info_AML <- merge_infoAml_drugid

#load("data/doseViabilityArray_AML.RData")
info_AML <- info_AML[dimnames(doseViabilityArray_AML)[[1]], ]

##################################################################################################################################################################
###########CREATING PROFILES#############

#coded on 18-jan for second version of pset, with replicates
#AUC_IC50_AML <- PharmacoGx:::.calculateFromRaw(doseViabilityArray_AML, cap = NA, nthread = 2, family = c("normal","Cauchy"), scale = 0.07, n=1 )
#AUC_IC50_AML <- PharmacoGx:::.calculateFromRaw(doseViabilityArray_AML, cap = NA, nthread = 2, n=1)
AUC_IC50_AML <-
    PharmacoGx:::.calculateFromRaw(
        doseViabilityArray_AML,
        cap = 100,
        nthread = 2,
        n = 1
    )

#saveRDS(AUC_IC50_AML, "data/AUC_IC50_AML.rds")
#load pre-run auc and ic50
#AUC_IC50_AML <- readRDS("data/AUC_IC50_AML.rds")
#common id col
xauc <- as.data.frame(AUC_IC50_AML$AUC)
xauc$id <- rownames(xauc)
xic50 <- as.data.frame(AUC_IC50_AML$IC50)
xic50$id <- rownames(xic50)

read_drugResponse <-
    read_excel("/pfs/beatAML_raw/Table-S10-Drug-Responses.xlsx", sheet = "Sheet 1")
read_drugResponse.3 <- read_drugResponse
read_drugResponse.3$lab_id <-
    paste("A", read_drugResponse.3$lab_id, sep = "_")
merge_auc_ic50_aml <-
    merge(
        new_ids_AML,
        read_drugResponse.3,
        by.x = c("lab_id", "inhibitor"),
        by.y = c("lab_id", "inhibitor"),
        all.x = TRUE
    )

#subset_list_aml <- c("merged_id", "ic50","auc")
subset_merge_auc_ic50_aml <-
    merge_auc_ic50_aml[, c("merged_id", "ic50", "auc"), drop = F]
unique_subset_merge_auc_ic50_aml <-
    unique(subset_merge_auc_ic50_aml)
merge_pub_recom_AML <-
    merge(unique_subset_merge_auc_ic50_aml,
          xauc,
          by.x = "merged_id",
          by.y = "id")
merge_pub_recom_AML <-
    merge(merge_pub_recom_AML, xic50, by.x = "merged_id", by.y = "id")

profiles_AML <- merge_pub_recom_AML
rownames(profiles_AML) <- profiles_AML$merged_id
profiles_AML <- profiles_AML[dimnames(doseViabilityArray_AML)[[1]], ]
profiles_AML <-
    profiles_AML[, -grep("merged_id", colnames(profiles_AML))]

colnames(profiles_AML) <-
    c("ic50_published",
      "auc_published",
      "aac_recomputed",
      "ic50_recomputed")

###########CREATING PSet#############

BeatAML <- PharmacoGx::PharmacoSet(
    "BeatAML",
    molecularProfiles <-
        c(kallisto_rnaseq_all,list(
            "mutationall" = new_SE_ALL,
            "mutationchosen" = new_SE_Chosen
        )),
    cell <- cell_AML,
    drug <- drug_AML,
    sensitivityInfo <- info_AML,
    sensitivityRaw <-
        doseViabilityArray_AML,
    sensitivityProfiles <-
        profiles_AML,
    curationDrug = curationDrug,
    curationCell = curationCell,
    datasetType <- "sensitivity",
    verify = TRUE
)
BeatAML@annotation$notes <-
    "1. All cellid(labIDs) in the PSet have a prefix 'A_'. 2. Mutation status from mutect and varscan is concatenated by '///' in the expression matrix of each object.3. All cell and drug metadata can be found in 'cell' and 'drug' objects respectively. 4. Dose values are equated as per the additional inputs from authors. 5. Throughout the 'sensitivity' object, a unique identifier has been created by concatenating drug-sample-inhibitor_panel-run_index replicate-time for representing each experiment; these are the rownames of each sub object. 6. All raw dose and viability values are in the 'sensitivity-raw' object. 7. 'sensitivity-profiles have measures like AUC,AAC,IC50. The new data has more unique replicate sets that are summarized in published data from paper. Published IC50 and AUC values will be repeated to best accomodate the new dose response data from authors. This structure might change in the future to accomodate these data more efficiently. 8. rnaseq and isoform values are log2(TPM+0.001) values. Counts are log2(count+1) transformed."
BeatAML@annotation$version <- 2
saveRDS(BeatAML, file = "/pfs/out/BeatAML.rds")
# write.csv(BeatAML@cell$cellid, "output/cell.csv")
# write.csv(rownames(BeatAML@drug), "output/lab.drug.csv")
# write.csv(BeatAML@drug$drug.name, "output/dataset.drug.csv")
