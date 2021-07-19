library(dplyr)

annotated_phe_position <- function(){
    data.table::fread(paste0("/ceph/projects/v4_phase3_paper/inputs/",
                             "counts/text_files_counts/_m/hippocampus/gene.bed")) %>%
        mutate(chr=gsub("chr", "", seqnames)) %>% select(gene_id, chr, start, end) %>%
        data.table::fwrite("phe.position.txt", sep='\t')
}

## Get annotations for gene positions
annotated_phe_position()
