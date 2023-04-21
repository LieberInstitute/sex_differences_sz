"""
This script prepares the eQTL results (tensorQTL) for mash modeling.
"""
import argparse
import pandas as pd
from pyhere import here

def load_eqtl(filename):
    df = pd.read_csv(filename, sep='\t', nrows=100,
                     usecols=["phenotype_id","variant_id","pval_gi"])
    return pd.read_csv(filename, sep='\t', dtype=df.dtypes.to_dict(),
                       usecols=["phenotype_id","variant_id","pval_gi"],
                       compression="gzip")


def extract_eqtls(feature, fn):
    ## Load eQTLs for mashr
    ### Caudate
    print("Loading Caudate data!")
    cc_file = here(f"prep_eqtl_analysis/caudate/{feature}",
                   f"interaction_model/_m/{fn}")
    caudate = load_eqtl(cc_file)
    ### DLPFC
    print("Loading DLPFC data!")
    dd_file = here(f"prep_eqtl_analysis/dlpfc/{feature}",
                   f"interaction_model/_m/{fn}")
    dlpfc = load_eqtl(dd_file)
    ### Hippocampus
    print("Loading hippocampus data!")
    hh_file = here(f"prep_eqtl_analysis/hippocampus/{feature}",
                   f"interaction_model/_m/{fn}")
    hippo = load_eqtl(hh_file)
    return caudate, dlpfc, hippo


def extract_dataframe(feature):
    fn1 = "BrainSEQ_TOPMed.cis_qtl_top_assoc.txt.gz"
    print("Load data!")
    ## Top eQTL pairs
    cc, dd, hh = extract_eqtls(feature, fn1)
    # Combine eQTL
    print("Extract strong set!")
    df = pd.concat([cc, dd, hh])
    df["effect"] = df.phenotype_id + "_" + df.variant_id
    # Save files
    df.drop_duplicates(subset="effect")\
      .to_csv("top_eqtls_interaction.txt.gz",sep='\t',index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Main
    extract_dataframe(args.feature)


if __name__=='__main__':
    main()
