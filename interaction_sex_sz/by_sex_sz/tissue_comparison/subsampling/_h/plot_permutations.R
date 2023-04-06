#### This script generates boxplots for permutation analysis
#### (male downsampling to female levels)

suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
})

save_ggplots <- function(p, fn, w=6, h=6){
    for(ext in c('.svg', '.png', '.pdf')){
        ggsave(p, filename=paste0(fn, ext), width=w, height=h)
    }
}

read_data <- function(region){
    fn <- paste0("../../../", tolower(region),
                 "/subsampling_male/deg_summary/_m/permutations.csv")
    return(data.table::fread(fn) |>
           select(Symbol, gencodeID, logFC, P.Value, adj.P.Val, Permutation) |>
           mutate(Tissue = region))
}

merge_data <- function(){
    lt <- c("Caudate", "DLPFC", "Hippocampus")
    return(purrr::map(lt, read_data) |> bind_rows())
}
memDF <- memoise::memoise(merge_data)

get_perm_data <- function(region){
    return(memDF() |> group_by(Permutation, Tissue) |>
           summarize(Size = n()) |>
           tidyr::replace_na(list(DLPFC=0, Hippocampus=0, Caudate=0)) |>
           filter(Tissue == region) |> as.data.frame())
}

get_female_degs <- function(region){
    fn <- paste0("../../../", tolower(region),
                 "/female_analysis/_m/genes/diffExpr_szVctl_full.txt")
    return(data.table::fread(fn) |> filter(`adj.P.Val` < 0.05))
}

cal_zscore <- function(region){
    x     <- dim(get_female_degs(region))[1]
    mu    <- mean(get_perm_data(region)$Size)
    sigma <- sd(get_perm_data(region)$Size)
    return((x - mu) / sigma)
}

cal_significant <- function(){
    z_scores = c(); two_tail = c(); tissues = c();
    for(region in c("Caudate", "DLPFC", "Hippocampus")){
        tissues = c(tissues, region)
                                        # Z-score
        q <- cal_zscore(region); z_scores <- c(z_scores, q)
                                        # Convert to p-value
        two_tail <- c(two_tail, 2*pnorm(q, mean=0, sd=1, lower.tail=FALSE))
    }
    return(data.frame("Tissue"=tissues, "Z_score"=z_scores,
                      "P_Value"=two_tail))
}

plot_hist <- function(){
    df <- memDF() |> group_by(Permutation, Tissue) |>
        summarize(Size = n()) |> as.data.frame() |>
        tidyr::pivot_wider(names_from=Tissue, values_from=Size) |>
        tidyr::replace_na(list(DLPFC=0, Hippocampus=0, Caudate=0)) |>
        tidyr::pivot_longer(-Permutation, names_to="Tissue",
                            values_to="DEGs") |>
        mutate_if(is.character, as.factor)
    hist <- gghistogram(df, x="DEGs", fill="lightgray", bins=30,
                        rug=TRUE, facet.by="Tissue", ncol=1,
                        ylab="Count in Perumtation",
                        xlab="Number of SZ DEGs\n(Subsampled Male Only)",
                        panel.labs.font=list(face="bold", size=18),
                        ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        font("xy.title", face="bold", size=18)
    return(hist)
}

#### Main
                                        # Calculate permutation p-values
cal_significant() |>
    data.table::fwrite("permutation_pvalues.tsv", sep='\t')
                                        # Plot histogram
hist <- plot_hist()
save_ggplots(hist, "permutation_histogram", 6, 7)
                                        # Print summary
dx <- memDF() |> group_by(Permutation, Tissue) |>
    summarize(Size = n()) |> as.data.frame() |>
    tidyr::pivot_wider(names_from=Tissue, values_from=Size) |>
    tidyr::replace_na(list(DLPFC=0, Hippocampus=0, Caudate=0)) |> 
    tidyr::pivot_longer(-Permutation, names_to="Tissue", values_to="DEGs") |>
    mutate_if(is.character, as.factor) |> group_by(Tissue) |>
    summarize(Mean = mean(DEGs), Median = median(DEGs), Std = sd(DEGs))
print(dx)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
