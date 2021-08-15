## Prepare covariates for FastQTL analysis
suppressPackageStartupMessages({
    library(tidyverse)
    library(optparse)
})

                                        # Create parser object
option_list <- list(
    make_option(c("-f", "--feature"), type="character", default="genes",
                help="Feature to extract covariates [default: %(default)]",
                metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

get_covs <- function(feature){
    load("covariates.rda", verbose=TRUE)
    if(feature == "genes"){
        covs0 = covsGene0
    }else if(feature == "transcripts"){
        covs0 = covsTx0
    } else if (feature == "exons"){
        covs0 = covsExon0
    } else {
        covs0 = covsJxn0
    }
    return(covs0)
}

extract_covs <- function(feature){
                                        # Load data
    sample_df = data.table::fread("../../_m/sample_id_to_brnum.tsv")
                                        # Extract covariates
    covs = get_covs(feature) %>% as.data.frame %>% select(sample_df$RNum) %>%
        t %>% as.data.frame %>% rownames_to_column("RNum") %>%
        inner_join(sample_df, by=c("RNum")) %>% select(-"RNum") %>%
        rename("ID"="BrNum") %>% column_to_rownames("ID") %>% t %>%
        as.data.frame %>% rownames_to_column("ID")
    return(covs)
}


if(length(opt$feature) == 0){
    print("Missing feature!")
} else {
    fn = paste(opt$feature, "combined_covariates.txt", sep='.')
    covs = extract_covs(opt$feature)
    covs %>% data.table::fwrite(fn, sep='\t')
}
