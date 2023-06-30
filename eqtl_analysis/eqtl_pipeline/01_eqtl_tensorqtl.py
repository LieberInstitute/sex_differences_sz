"""
This script runs tensorQTL in python.
"""
import pandas as pd
import argparse, session_info
from functools import lru_cache
from tensorqtl import cis, read_phenotype_bed
from tensorqtl import calculate_qvalues, genotypeio

@lru_cache()
def get_genotypes():
    plink_prefix_path = "../../plink_format/_m/genotypes"
    pr = genotypeio.PlinkReader(plink_prefix_path)
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
    variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
    return pr.load_genotypes(fam_id="fid"), variant_df


@lru_cache()
def get_interaction():
    f_inter = "../../_m/sex_interaction.txt"
    return pd.read_csv(f_inter, sep='\t', index_col=0).loc[:,["Sex"]]


@lru_cache()
def get_covars(feature):
    covar_file = f"../../covariates/_m/{feature}.combined_covariates.txt"
    return pd.read_csv(covar_file, sep='\t', index_col=0)\
             .transpose()


@lru_cache()
def get_phenotype(feature):
    expr_bed = f"../../normalize_expression/_m/{feature}.expression.bed.gz"
    return read_phenotype_bed(expr_bed)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    feature = args.feature
    
    # Load data
    phenotype_df, phenotype_pos_df = get_phenotype(feature)
    genotype_df, variant_df = get_genotypes()
    prefix = "BrainSEQ_TOPMed"

    # Interaction analysis
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                    prefix, covariates_df=get_covars(feature),
                    interaction_df=get_interaction(), maf_threshold=0.05,
                    window=500000, maf_threshold_interaction=0.05,
                    output_dir=".", run_eigenmt=True,
                    write_top=True, write_stats=True)
    
    # Reproducibility information
    session_info.show()

     
if __name__ == "__main__":
    main()
