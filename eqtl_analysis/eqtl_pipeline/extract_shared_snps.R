                                        # Select common SNPs between the sexes
library(dplyr)

ff = data.table::fread("../../_m/LIBD_TOPMed_female.pvar")
mm = data.table::fread("../../_m/LIBD_TOPMed_male.pvar")
shared = ff %>% inner_join(mm, by=c("#CHROM", "POS", "ID", "REF", "ALT")) %>%
    select(ID) %>% data.table::fwrite("shared_snps.tsv", sep='\t', col.names=F)
