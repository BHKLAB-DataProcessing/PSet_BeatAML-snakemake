from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"],
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]
filename = config["filename"]

rule get_pset:
    input:
        S3.remote(prefix + "processed/new_SE_ALL.rds"),
        S3.remote(prefix + "processed/new_SE_Chosen.rds"),
        S3.remote(prefix + "processed/cell_AML.rds"),
        S3.remote(prefix + "processed/drug_AML.rds"),
        S3.remote(prefix + "processed/kallisto_rnaseq_all.rds"),
        S3.remote(prefix + "processed/curationDrug.rds"),
        S3.remote(prefix + "processed/curationCell.rds"),
        S3.remote(prefix + "processed/info_AML.rds"),
        S3.remote(prefix + "processed/doseViabilityArray_AML.rds"),
        S3.remote(prefix + "processed/profiles_AML.rds")
    output:
        S3.remote(prefix + filename)
    shell:
        """
        Rscript scripts/getBeatAML.R {prefix} {filename}
        """

rule process_sensitivity:
    input:
        S3.remote(prefix + "processed/curationDrug.rds"),
        S3.remote(
            prefix + "download/beataml_manuscript_raw_inhib_v2_passed_qc_ceil_100.tsv"),
        S3.remote(prefix + "download/report_1614271471952.csv"),
        S3.remote(prefix + "download/Table-S10-Drug-Responses.xlsx")
    output:
        S3.remote(prefix + "processed/info_AML.rds"),
        S3.remote(prefix + "processed/doseViabilityArray_AML.rds"),
        S3.remote(prefix + "processed/profiles_AML.rds")
    shell:
        """
        Rscript scripts/process_sensitivity.R {prefix}
        """

rule process_seq_and_curation:
    input:
        S3.remote(prefix + "download/beatAML_SE.rds"),
        S3.remote(prefix + "download/Table-S5-Clinical-Summary.xlsx"),
        S3.remote(prefix + "download/Table-S10-Drug-Responses.xlsx"),
        S3.remote(prefix + "download/drugs_with_ids.csv")
    output:
        S3.remote(prefix + "processed/cell_AML.rds"),
        S3.remote(prefix + "processed/drug_AML.rds"),
        S3.remote(prefix + "processed/kallisto_rnaseq_all.rds"),
        S3.remote(prefix + "processed/curationDrug.rds"),
        S3.remote(prefix + "processed/curationCell.rds")
    shell:
        """
        Rscript scripts/process_seq_and_curation.R {prefix}
        """

rule process_se:
    input:
        S3.remote(prefix + "download/Table-S7-Variants-for-Analysis.xlsx")
    output:
        S3.remote(prefix + "processed/new_SE_ALL.rds"),
        S3.remote(prefix + "processed/new_SE_Chosen.rds")
    shell:
        """
        Rscript scripts/process_SE.R {prefix}
        """

rule download_data:
    output:
        S3.remote(prefix + "download/Table-S5-Clinical-Summary.xlsx"),
        S3.remote(prefix + "download/Table-S7-Variants-for-Analysis.xlsx"),
        S3.remote(prefix + "download/Table-S10-Drug-Responses.xlsx"),
        S3.remote(prefix + "download/beatAML_SE.rds"),
        S3.remote(
            prefix + "download/beataml_manuscript_raw_inhib_v2_passed_qc_ceil_100.tsv"),
        S3.remote(prefix + "download/report_1614271471952.csv"),
        S3.remote(prefix + "download/drugs_with_ids.csv")
    shell:
        """
        Rscript scripts/downloadData.R {prefix}
        """
