## Plotting RXE after removing specific XCI status genes.

library(here)
library(dplyr)
library(ggpubr)


save_ggplots <- function(p, fn, w=7, h=7){
    for(ext in c('.svg', '.pdf')){
        ggsave(p, filename=paste0(fn, ext), width=w, height=h)
    }
}

get_pheno <- function(){
    fn <- here("input/phenotypes/_m/phenotypes.csv")
    return(data.table::fread(fn) |>
           filter(Region == "HIPPO", Sex == "M") |>
           distinct() |> select(c("RNum", "Sex", "Dx", "Age")))
}

get_degs <- function(){
    fn <- here("differential_expression/hippocampus",
               "_m/genes/diffExpr_maleVfemale_full.txt")
    return(data.table::fread(fn) |>
           filter(`adj.P.Val` < 1) |>
           rename("feature_id"="V1",
                  "gene_name"="Symbol",
                  "gencode_id"="gencodeID") |>
           select(c("feature_id", "gencode_id", "gene_name",
                    "logFC", "adj.P.Val", "SE")))
}

get_xci <- function(){
    fn <- here("xci_dosage_analysis/xci_enrichment/_h",
               "xci_status_hg19.txt")
    return(data.table::fread(fn) |>
           select(c("Gene name", "Combined XCI status")) |>
           rename("gene_name"="Gene name",
                  "xci_status"="Combined XCI status") |>
           inner_join(get_degs(), by="gene_name"))
}

xci_annotation <- function(){
    fn <- here("input/counts/text_files_counts/_m",
               "hippocampus/gene_annotation.txt")
    return(data.table::fread(fn) |>
           select(c("name", "seqnames")) |>
           rename("feature_id"="name", "chrom"="seqnames") |>
           inner_join(get_xci(), by="feature_id"))
}


get_annotation <- function(){
    fn <- here("input/counts/text_files_counts/_m",
               "hippocampus/gene_annotation.txt")
    return(data.table::fread(fn) |>
           select(c("name", "gene_name", "seqnames")) |>
           rename("feature_id"="name", "chrom"="seqnames"))
}

get_logTPM <- function(){
    fn <- here("input/counts/text_files_counts/tpm/_m",
               "hippocampus/gene.log2tpm.csv")
    return(data.table::fread(fn) |>
           rename("feature_id"="name") |>
           inner_join(get_annotation(), by="feature_id") |>
           mutate(chrom_type=case_when(
                      chrom == "chrX" ~ "X",
                      stringr::str_detect(chrom, "chr\\d+") ~ "Autosome",
                      TRUE ~ "Other")
                  ))
}
memTPM <- memoise::memoise(get_logTPM)

subset_data <- function() {
    logTPM  <- memTPM()
    samples <- logTPM |> select(starts_with("R")) |> names()
    return(mutate(logTPM,
                  sum=rowSums(select(logTPM, all_of(samples)))) |>
           filter(sum >= (0.2 * length(samples)) &
                  chrom_type %in% c("X", "Autosome")) |>
           select(-gene_name, -chrom, -sum))
}
memDF <- memoise::memoise(subset_data)

cal_rxe <- function(XCI, STATUS){
    if (XCI) {
        if (!is.null(STATUS)) {
            xci <- xci_annotation() |> filter(xci_status == STATUS)
            df  <- anti_join(memDF(), xci, by="feature_id")
        } else {
            df <- anti_join(memDF(), xci_annotation(),
                            by="feature_id")
        }
    } else {
        df <- memDF()
    }
    return(group_by(df, chrom_type) |>
           summarise(across(where(is.numeric), \(x) mean(x, na.rm=TRUE))) |>
           tidyr::pivot_longer(cols=-chrom_type, names_to="RNum",
                               values_to="value") |>
           tidyr::pivot_wider(names_from="chrom_type",
                              values_from="value") |>
           mutate(RXE= X - Autosome))
}

plotting_rxe <- function(XCI, STATUS){
    dx  <- inner_join(cal_rxe(XCI, STATUS), get_pheno(),
                     by="RNum") |> distinct()
    bxp <- ggboxplot(dx, x="Dx", y="RXE", add="jitter", fill="Dx",
                     palette=c("#999995ff", "#e69f00ff"),
                     outlier.shape=NA, ylab="Relative X Expression",
                     xlab="Diagnosis", legend="bottom",
                     ggtheme=theme_pubr(base_size=15)) +
        stat_compare_means(comparison=list(c("Control", "SCZD"))) +
        font("xy.title", size=16, face="bold")
    if(XCI){
        fname <- paste0("boxplot_rxe.no_", STATUS, ".xci_genes")
    } else {
        fname <- paste0("boxplot_rxe.all_genes")
    }
    save_ggplots(bxp, fname, 5, 5)
}

#### MAIN

plotting_rxe(FALSE, NULL)
plotting_rxe(TRUE, "escape")
plotting_rxe(TRUE, "variable")
plotting_rxe(TRUE, "inactive")

#### Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
