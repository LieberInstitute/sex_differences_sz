## Implementation of RRHO for sex-specific expression

library(RRHO2)
library(dplyr)

get_de <- function(){
    err <- 1e-320
    fn  <- paste0("../../summary_table/_m/",
                  "differential_expression_schizophrenia_by_sex_4features.txt.gz")
    return(data.table::fread(fn) |>
           filter(Type == "Gene") |>
           mutate(Score = -log10(`P.Value` + err) * sign(logFC)) |>
           select(Tissue, Sex, gencodeID, Score))
}
memDE <- memoise::memoise(get_de)

subset_sex <- function(sex){
    return(filter(memDE(), Sex == sex) |> select(-Sex))
}

subset_regions <- function(tissue1, tissue2, sex){
    df1 <- subset_sex(sex) |> filter(Tissue == tissue1) |>
        select(-Tissue)
    df2 <- subset_sex(sex) |> filter(Tissue == tissue2) |>
        select(-Tissue)
    genes <- intersect(df1$gencodeID, df2$gencodeID)
    return(list(tissue1=filter(df1, gencodeID %in% genes),
                tissue2=filter(df2, gencodeID %in% genes)))
}

run_rrho <- function(tissue1, tissue2, sex){
    lt       <- subset_regions(tissue1, tissue2, sex)
    rrho_obj <- RRHO2_initialize(lt[["tissue1"]], lt[["tissue2"]],
                                 labels=c(tissue1, tissue2),
                                 log10.ind=TRUE)
    outfile  <- paste("rrho_heatmap", sex, tissue1, tissue2,
                      "pdf", sep=".")
    pdf(file = tolower(outfile), width = 10, height = 10)
    RRHO2_heatmap(rrho_obj, maximum=600, minimum=0,
                  colorGradient=viridis::viridis(100))
    dev.off()
}

within_tissue_rrho <- function(tissue){
    df1 <- filter(subset_sex("Female"), Tissue == tissue) |>
        select(-Tissue)
    df2 <- filter(subset_sex("Male"), Tissue == tissue) |>
        select(-Tissue)
    genes <- intersect(df1$gencodeID, df2$gencodeID)
    rrho_obj <- RRHO2_initialize(filter(df1, gencodeID %in% genes),
                                 filter(df2, gencodeID %in% genes),
                                 labels=c("Female", "Male"),
                                 log10.ind=TRUE)
    outfile  <- paste0("rrho_sex_comparison.",tissue,".pdf")
    pdf(file = tolower(outfile), width = 10, height = 10)
    RRHO2_heatmap(rrho_obj, maximum=600, minimum=0,
                  colorGradient=viridis::viridis(100))
    dev.off()
}

####### MAIN
for(sex in c("Female", "Male")){
    run_rrho("Caudate", "DLPFC", sex)
    run_rrho("Caudate", "Hippocampus", sex)
    run_rrho("DLPFC", "Hippocampus", sex)
}

within_tissue_rrho("Caudate")
within_tissue_rrho("DLPFC")
within_tissue_rrho("Hippocampus")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
