## Explore eigengene values and correlate with phenotype
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(tidyverse)
    library(ggpubr)
})

pheno = data.table::fread(paste0("/ceph/projects/v4_phase3_paper/inputs/",
                                 "phenotypes/_m/merged_phenotypes.csv"))
eigen = data.table::fread("../../_m/eigengenes.csv")
modules = eigen %>% select(-V1) %>% colnames
dt = eigen %>% left_join(pheno, by=c("V1"="RNum")) %>%
    mutate_at("Sex", as.factor) %>% mutate_at("Sex", as.numeric)

pvals = c()
for(mod in modules){
    model = paste0("Sex ~ ", mod)
    res = anova(lm(model, data=dt))
    pvals = c(pvals, res[mod, "Pr(>F)"])
}
fdr <- p.adjust(pvals, method="fdr")
df1 = data.frame("Modules"=modules, "Pvalue"=pvals, "FDR"=fdr)
df1 %>% mutate(Tissue="Caudate") %>%
    data.table::fwrite("eigen_correlation_sex.tsv", sep='\t')

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
