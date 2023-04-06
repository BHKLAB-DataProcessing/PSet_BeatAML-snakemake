library(PharmacoGx)
library(Biobase)
library(BiocGenerics)
library(readxl)
library(SummarizedExperiment)
library(abind)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- paste0(args[[1]], "download")
processed_dir <- paste0(args[[1]], "processed")

curationDrug <- readRDS(file.path(processed_dir, "curationDrug.rds"))

# download_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_BeatAML-snakemake/bhklab_orcestra/snakemake/PSet_BeatAML/download"
# processed_dir <- "/Users/minoru/Code/bhklab/DataProcessing/PSet/PSet_BeatAML-snakemake/bhklab_orcestra/snakemake/PSet_BeatAML/processed"

#########################################################################################################################################################################
######## Creating drug object########
# raw_data_pointsAML <- read.delim("data/inhibitor_data_points_2018_10_24.txt",sep="\t", header=TRUE )#old
# csv with conc values that needs to be changed
# concentrationMapping <- read.csv(("data/concentration_mapping.csv"), header = TRUE,fileEncoding="UTF-8-BOM")#old

# PROCESSING RAW DOSE RESPONSE - updated data from Jeffrey Tyner OHSU
raw_data_pointsAML <-
  read.table(
    file.path(download_dir, "beataml_manuscript_raw_inhib_v2_passed_qc_ceil_100.tsv"),
    sep = "\t",
    header = TRUE
  )
# most recent panel version of cleaned up concentrations with consistent rounding from Jeffrey Tyner
clean_conc <-
  read.csv(file.path(download_dir, "report_1614271471952.csv"), stringsAsFactors = F) # only for reference

# equating dose values to ref values sent by Jeff
# checks to see if merging values would create duplicated rows
check_val <- function(raw.df, val1, val2) {
  rr <-
    raw.df[raw.df$well_concentration == val1 |
      raw.df$well_concentration == val2, ]
  if (length(unique(rr$well_concentration)) != "2") {
    print("Not all input values are present in subset")
  } else {
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

# for values 0.001372 &  0.001400
val.pair.1 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "0.001372",
    val2 = "0.0014"
  )
# for values 0.004100 &  0.004115
val.pair.2 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "0.0041",
    val2 = "0.004115"
  )
# for values 0.012300 & 0.012346
val.pair.3 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "0.0123",
    val2 = "0.012346"
  )
# for values 0.013700 & 0.013717
val.pair.4 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "0.0137",
    val2 = "0.013717"
  )
# for values 0.037000 & 0.037037
val.pair.5 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "0.037",
    val2 = "0.037037"
  )
# for values 0.041152 & 0.041200
val.pair.6 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "0.041152",
    val2 = "0.0412"
  )
# for values 0.111100 & 0.111111
val.pair.7 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "0.1111",
    val2 = "0.111111"
  )
# for values 0.123457 & 0.123500
val.pair.8 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "0.123457",
    val2 = "0.1235"
  )
# for values 0.333300 & 0.333333
val.pair.9 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "0.3333",
    val2 = "0.333333"
  )
# for values 0.370370 & 0.370400
val.pair.10 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "0.37037",
    val2 = "0.3704"
  )
# for values 1.111100 & 1.111110
val.pair.11 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "1.1111",
    val2 = "1.11111"
  )
# for values 3.333300 & 3.333330
val.pair.12 <-
  check_val(
    raw.df = raw_data_pointsAML,
    val1 = "3.3333",
    val2 = "3.33333"
  )


updatedRaw_data_pointAML <- raw_data_pointsAML

# new conversions
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
########### CREATING 3D ARRAY#############

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
    new_ids_AML$time_of_read,
    sep = "."
  )

unique_exp_idsAML <- unique(new_ids_AML$merged_id)

# Creating dose Array
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

# Creating viability Array
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
  tdf <- new_ids_AML[new_ids_AML$merged_id == rid, ]
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
# saveRDS(doseViabilityArray_AML, "data/doseViabilityArray_AML.rds")
# load array
# doseViabilityArray_AML <- readRDS("data/doseViabilityArray_AML.rds")
########### CREATING INFO#############
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
    by.y = "AML.drugid"
  )
# merge_infoAml_drugid <- merge_infoAml_drugid[,-1]
# colnames(merge_infoAml_drugid)[3] <- "drugid"
names(merge_infoAml_drugid)[names(merge_infoAml_drugid) == "unique.drugid"] <-
  "drugid"

merge_infoAml_drugid[, "nbr.conc.tested"] <- 21
# colnames(merge_infoAml_drugid)[2] <- "cellid"
names(merge_infoAml_drugid)[names(merge_infoAml_drugid) == "lab_id"] <-
  "cellid"

rownames(merge_infoAml_drugid) <- merge_infoAml_drugid$merged_id
merge_infoAml_drugid <-
  merge_infoAml_drugid[, -grep("merged_id", colnames(merge_infoAml_drugid))]
info_AML <- merge_infoAml_drugid

# load("data/doseViabilityArray_AML.RData")
info_AML <- info_AML[dimnames(doseViabilityArray_AML)[[1]], ]

##################################################################################################################################################################
########### CREATING PROFILES#############

# coded on 18-jan for second version of pset, with replicates
# AUC_IC50_AML <- PharmacoGx:::.calculateFromRaw(doseViabilityArray_AML, cap = NA, nthread = 2, family = c("normal","Cauchy"), scale = 0.07, n=1 )
# AUC_IC50_AML <- PharmacoGx:::.calculateFromRaw(doseViabilityArray_AML, cap = NA, nthread = 2, n=1)
AUC_IC50_AML <-
  PharmacoGx:::.calculateFromRaw(
    doseViabilityArray_AML,
    cap = 100,
    nthread = 2,
    n = 1
  )

# saveRDS(AUC_IC50_AML, "data/AUC_IC50_AML.rds")
# load pre-run auc and ic50
# AUC_IC50_AML <- readRDS("data/AUC_IC50_AML.rds")
# common id col
xauc <- as.data.frame(AUC_IC50_AML$AUC)
xauc$id <- rownames(xauc)
xic50 <- as.data.frame(AUC_IC50_AML$IC50)
xic50$id <- rownames(xic50)

read_drugResponse <-
  read_excel(file.path(download_dir, "Table-S10-Drug-Responses.xlsx"), sheet = "Sheet 1")
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

# subset_list_aml <- c("merged_id", "ic50","auc")
subset_merge_auc_ic50_aml <-
  merge_auc_ic50_aml[, c("merged_id", "ic50", "auc"), drop = F]
unique_subset_merge_auc_ic50_aml <-
  unique(subset_merge_auc_ic50_aml)
merge_pub_recom_AML <-
  merge(unique_subset_merge_auc_ic50_aml,
    xauc,
    by.x = "merged_id",
    by.y = "id"
  )
merge_pub_recom_AML <-
  merge(merge_pub_recom_AML, xic50, by.x = "merged_id", by.y = "id")

profiles_AML <- merge_pub_recom_AML
rownames(profiles_AML) <- profiles_AML$merged_id
profiles_AML <- profiles_AML[dimnames(doseViabilityArray_AML)[[1]], ]
profiles_AML <-
  profiles_AML[, -grep("merged_id", colnames(profiles_AML))]

colnames(profiles_AML) <-
  c(
    "ic50_published",
    "auc_published",
    "aac_recomputed",
    "ic50_recomputed"
  )

saveRDS(info_AML, file.path(processed_dir, "info_AML.rds"))
saveRDS(doseViabilityArray_AML, file.path(processed_dir, "doseViabilityArray_AML.rds"))
saveRDS(profiles_AML, file.path(processed_dir, "profiles_AML.rds"))
