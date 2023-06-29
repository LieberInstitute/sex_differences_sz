suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(argparse)
})

## Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character",
                    help="feature name used in eQTL analysis")
parser$add_argument("-p", "--perm_pval", type="double", default=0.0001,
                    help="eQTL p-value threshold [default: %default]")
args <- parser$parse_args()

## Functions
save_ggplots <- function(p, fn, w=6, h=6){
    for(ext in c('.pdf', '.png', '.svg')){
        ggsave(p, filename=paste0(fn, ext), width=w, height=h)
    }
}

get_gwas <- function(feature, sex){
    fn <- paste0(sex, '/gwas_pgc3_', feature, ".txt")
    return(data.table::fread(fn) %>% mutate(Sex=sex))
}

get_eqtl <- function(feature, sex){
    fn <- paste0(sex, '/eqtl_', feature, '.txt')
    return(data.table::fread(fn) %>% mutate(Sex=sex))
}

merge_data <- function(feature, perm_pval){
    gwas_df <- bind_rows(get_gwas(feature, "female"),
                         get_gwas(feature, "male"))
    eqtl_df <- bind_rows(get_eqtl(feature, "female"),
                         get_eqtl(feature, "male"))
    return(gwas_df %>%
           inner_join(eqtl_df, by=c("SNP"="SNP.Id", "Sex")) %>%
           mutate("eqtl"=-log10(P.Value), "gwas"=-log10(P)) %>%
           filter(P.Value < perm_pval) %>%
           mutate_if(is.character, as.factor) %>%
           mutate(Sex=recode_factor(Sex, female="Female",
                                    male="Male")))
}

plotNsave_PP <- function(feature, perm_pval){
    sca <- merge_data(feature, perm_pval) %>%        
        ggscatter(x="eqtl", y="gwas", facet.by="Sex",
                  ylab="-log10(P[GWAS])", xlab="-log10(P[eQTL])",
                  panel.labs.font=list(face="bold", size=16),
                  add="reg.line", conf.int=TRUE, cor.coef=TRUE,
                  add.params=list(fill="lightgray"), scales="free_x",
                  cor.coeff.args=list(label.sep="\n", label.y=8), 
                  ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        font("xy.title", face="bold", size=16)
    save_ggplots(sca, paste0(feature, "_PP_plot"), 6, 4)
}

## Main
plotNsave_PP(args$feature, args$perm_pval)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
