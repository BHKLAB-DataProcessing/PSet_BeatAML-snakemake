require(downloader)
library(curl)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- paste0(args[1], "download")

root_url <- "https://orcestradata.blob.core.windows.net/beataml/"

curl_download(paste0(root_url, "Table-S5-Clinical-Summary.xlsx"), destfile = file.path(out_dir, "Table-S5-Clinical-Summary.xlsx"))
curl_download(paste0(root_url, "Table-S7-Variants-for-Analysis.xlsx"), destfile = file.path(out_dir, "Table-S7-Variants-for-Analysis.xlsx"))
curl_download(paste0(root_url, "Table-S10-Drug-Responses.xlsx"), destfile = file.path(out_dir, "Table-S10-Drug-Responses.xlsx"))
curl_download(paste0(root_url, "beatAML_SE.rds"), destfile = file.path(out_dir, "beatAML_SE.rds"))
curl_download(paste0(root_url, "beataml_manuscript_raw_inhib_v2_passed_qc_ceil_100.tsv"), destfile = file.path(out_dir, "beataml_manuscript_raw_inhib_v2_passed_qc_ceil_100.tsv"))
curl_download(paste0(root_url, "report_1614271471952.csv"), destfile = file.path(out_dir, "report_1614271471952.csv"))
curl_download("https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/drugs_with_ids.csv", destfile = file.path(out_dir, "drugs_with_ids.csv"))
