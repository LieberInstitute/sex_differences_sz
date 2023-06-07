#### This script explores sex steroid hormone levels
#### in the brain

suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
})

get_expression <- function(region){
    fn <- here::here("differential_expression",
                     tolower(region), "_m/genes/voomSVA.RData")
    load(fn)
    return(v)
}
memVOOM <- memoise::memoise(get_expression)

get_covs <- function(region){
    v      <- memVOOM(region)
    sample <- select(as.data.frame(v$targets), "Sex", "Dx") |>
        tibble::rownames_to_column("RNum")
    pheno  <- as.data.frame(v$design) |> select(-Intercept)
    return(pheno)
}

extract_sex_hormones <- function(region){
    v      <- memVOOM(region)
    return(v$genes |> as.data.frame() |>
           filter(Symbol %in% c("AR", "ESR1", "ESR2", "PGR")) |>
           select(gencodeID, Symbol) |>
           tibble::rownames_to_column("feature_id"))
}

get_norm_expr <- function(region){
    v     <- memVOOM(region)
    genes <- extract_sex_hormones(region)
    return(as.data.frame(t(v$E[genes$feature_id, ])))
}

get_res_expr <- function(region, protect_vars, FEMALE=TRUE){
    sex   <- ifelse(FEMALE, 0, 1)
    genes <- extract_sex_hormones(region)
    covs  <- get_covs(region) |> select(-any_of(protect_vars))
    norm  <- get_norm_expr(region)
    nsize <- dim(select(filter(covs, Male == sex), -Male))[1]
    res   <- norm[1:nsize, ]
    for(ii in 1:length(genes$feature_id)){
        d       <- as.data.frame(cbind(y=norm[,ii],covs)) |>
            filter(Male == sex) |> select(-Male) ## Select female
	model   <- lm(y ~ .,data=d)
	res[,ii] <- resid(model)
    }
    return(res)
}

merge_data <- function(region, protect_vars, FEMALE=TRUE){
    genes <- extract_sex_hormones(region)
    res   <- get_res_expr(region, protect_vars, FEMALE) |>
        tibble::rownames_to_column("RNum") |>
        tidyr::pivot_longer(-RNum, names_to="feature_id",
                            values_to="expression")
    covs  <- get_covs(region) |>
        select(any_of(protect_vars)) |>
        tibble::rownames_to_column("RNum")
    return(inner_join(res, covs, by="RNum") |>
           inner_join(genes, by="feature_id"))
}

loop_regions <- function(region){
    protect_vars <- c("Age", "SCZD")
    data1 <- merge_data(region, protect_vars, TRUE) |>
        mutate(sex="Female")
    data2 <- merge_data(region, protect_vars, FALSE) |>
        mutate(sex="Male")
    return( bind_rows(data1, data2) |> mutate(tissue = region) )
}

plot_corr <- function(df, SEX){
    fn <- paste(tolower(SEX), "age.hormones.pdf", sep=".")
    sca <- df |> filter(sex == SEX) |>
        mutate(Dx=ifelse(SCZD == 1, "SZ", "CTL")) |>
        ggscatter(x="Age", y="expression", color="Dx", add="reg.line", size=1,
                  ylab="Residualized Expression", alpha=0.3, palette="npg", 
                  panel.labs.font=list(face="bold"), ncol=4, cor.coef.size=3,
                  facet.by=c("tissue", "Symbol"), conf.int=FALSE, cor.coef=FALSE, 
                  ggtheme=theme_pubr(base_size=15, border=TRUE))
    ggsave(file=fn, plot=sca, width=10, height=6)
}

plot_corr_annot <- function(df, SEX){
    fn <- paste(tolower(SEX), "age.hormones.annot.pdf", sep=".")
    sca <- df |> filter(sex == SEX) |>
        mutate(Dx=ifelse(SCZD == 1, "SZ", "CTL")) |>
        ggscatter(x="Age", y="expression", color="Dx", add="reg.line", size=1,
                  ylab="Residualized Expression", alpha=0.3, palette="npg", 
                  panel.labs.font=list(face="bold"), ncol=4, cor.coef.size=3,
                  facet.by=c("tissue", "Symbol"), conf.int=TRUE, cor.coef=FALSE, 
                  ggtheme=theme_pubr(base_size=15, border=TRUE))
    ggsave(file=fn, plot=sca, width=10, height=6)
}

plot_expr <- function(df, SEX){
    fn <- paste(tolower(SEX), "age.hormones.boxplot.pdf", sep=".")
    bxp <- df |> filter(sex == SEX) |>
        mutate(Dx=ifelse(SCZD == 1, "SZ", "CTL")) |>
        ggboxplot(x="Dx", y="expression", add="jitter", palette="npg",
                  xlab="Diagnosis", ylab="Residualized Expression",
                  add.params=list(size=1, alpha=0.3), color="Dx",
                  panel.labs.font=list(face="bold"),
                  facet.by=c("tissue", "Symbol"),
                  ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        geom_pwc(aes(group=Dx))
    ggsave(file=fn, plot=bxp, width=10, height=8)
}

#### Main
                                        # Combine Data
df <- purrr::map(c("Caudate", "DLPFC", "Hippocampus"),
                 loop_regions) |> bind_rows()
plot_corr(df, "Female"); plot_corr(df, "Male")
plot_corr_annot(df, "Female"); plot_corr_annot(df, "Male")
plot_expr(df, "Female"); plot_expr(df, "Male")

#### Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
