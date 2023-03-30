## Implementation of RRHO for sex-specific expression

library(RRHO2)
library(dplyr)

get_de <- function(){
    err <- 1e-320
    fn  <- paste0("../../summary_table/_m/",
                  "differential_expression_analysis_4features_sex.txt.gz")
    return(data.table::fread(fn) |>
           filter(Type == "Gene") |>
           mutate(Score = -log10(`P.Value` + err) * sign(logFC)) |>
           select(Tissue, gencodeID, Score))
}
memDE <- memoise::memoise(get_de)

subset_regions <- function(tissue1, tissue2){
    df1 <- memDE() |> filter(Tissue == tissue1) |>
        select(-Tissue)
    df2 <- memDE() |> filter(Tissue == tissue2) |>
        select(-Tissue)
    genes <- intersect(df1$gencodeID, df2$gencodeID)
    return(list(tissue1=filter(df1, gencodeID %in% genes),
                tissue2=filter(df2, gencodeID %in% genes)))
}

run_rrho <- function(tissue1, tissue2){
    lt       <- subset_regions(tissue1, tissue2)
    rrho_obj <- RRHO2_initialize(lt[["tissue1"]], lt[["tissue2"]],
                                 labels=c(tissue1, tissue2),
                                 log10.ind=TRUE)
    outfile  <- paste("rrho_heatmap", tissue1, tissue2,
                      "pdf", sep=".")
    pdf(file = tolower(outfile), width = 10, height = 10)
    RRHO2_heatmap(rrho_obj, maximum=600, minimum=0,
                  colorGradient=viridis::viridis(100))
    dev.off()
}
    
####### MAIN
run_rrho("Caudate", "DLPFC")
run_rrho("Caudate", "Hippocampus")
run_rrho("DLPFC", "Hippocampus")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
