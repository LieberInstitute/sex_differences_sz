## This script generates plots for GSEA analysis.
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

save_plot <- function(p, fn, w, h){
    for(ext in c('.svg', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

get_top_GO <- function(tissue, sex){
    new_tissue <- tolower(tissue)
    filenames  <- list.files(path=tolower(sex), full.names=TRUE,
                             pattern=paste0(new_tissue, "_gsea*"))
    filenames  <- filenames[!grepl("DGN", filenames)]
    df_list <- list()
    for(i in seq_along(filenames)){
        df_list[[i]] <- data.table::fread(filenames[i])
    }
    return(bind_rows(df_list) |> arrange(pvalue) |> head(10) |>
           mutate(`Log10`=-log10(qvalue), Tissue=tissue, Sex=sex))
}

generate_dataframe <- function(){
    ##                                     # Female
    ## female_df <- get_top_GO("Caudate", "Female")
    ##                                     # Male
    df_list1 <- list(); df_list2 <- list()
    tissues  <- c("Caudate", "DLPFC", "Hippocampus")
    for(jj in seq_along(tissues)){
        df_list1[[jj]] <- get_top_GO(tissues[jj], "Male")
        df_list2[[jj]] <- get_top_GO(tissues[jj], "Female")
    }
    male_df   <- bind_rows(df_list1)
    female_df <- bind_rows(df_list2)
    return( bind_rows(female_df, male_df) )
}

plot_GO <- function(sex, CAUDATE=FALSE){
    dt <- generate_dataframe() |> filter(Sex == sex)
    if(CAUDATE) { dt <- filter(dt, Tissue == "Caudate") }
    cbPalette <- ggpubr::get_palette(palette = "npg", 4)
    gg1 = ggplot(dt, aes(x=NES, y=Description, color=Tissue, size=Log10)) +
        geom_point(shape=18, alpha=0.8) + xlim(-3, 3) +
        labs(y='', size='-Log10 (qvalue)') +
        scale_colour_manual(name="Brain Region", values=cbPalette,
                            labels=c("Caudate","DLPFC","Hippocampus")) +
        scale_size_continuous(range = c(2, 10)) +
        theme_bw(base_size=15) +
        theme(axis.title=element_text(face='bold'),
              strip.text=element_text(face='bold'))
    return(gg1)
}

#### MAIN
                                        # Females
gg = plot_GO("Female", TRUE)
save_plot(gg, "sex_sz_females.caudate.GSEA_top10_stacked", 9, 5)
gg = plot_GO("Female")
save_plot(gg, "sex_sz_females.GSEA_top10_stacked", 9, 5)

                                        # Males
gg = plot_GO("Male")
save_plot(gg, "sex_sz_males.GSEA_top10_stacked", 14, 7)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
