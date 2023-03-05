## Explore eigengene values and correlate with phenotype
suppressPackageStartupMessages({
    library(dplyr)
})

fn      <- here::here("input/phenotypes/_m/phenotypes.csv")
pheno   <- data.table::fread(fn)
eigen   <- data.table::fread("../../_m/eigengenes.csv")
modules <- eigen |> select(-V1) |> colnames()
dt      <- left_join(eigen, pheno, by=c("V1"="RNum"), multiple="all") |>
    mutate_at("Sex", as.factor) |> mutate_at("Sex", as.numeric)

pvals = c(); rsq = c();
for(mod in modules){
    model = paste0("Sex ~ ", mod)
    res = summary(lm(model, data=dt))
    rsq = c(rsq, res$r.squared)
    pvals = c(pvals, res$coefficients[mod, "Pr(>|t|)"])
}
fdr <- p.adjust(pvals, method="fdr")
df1 <- data.frame("Modules"=modules, "Rsq"=rsq,
                  "Pvalue"=pvals, "FDR"=fdr)
print(df1 |> filter(Pvalue < 0.05))
df1 |> mutate(Tissue="Hippocampus") |>
    data.table::fwrite("eigen_correlation_sex.tsv", sep='\t')

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
