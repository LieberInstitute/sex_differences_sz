"""
This script prepares the eQTL results (tensorQTL) for mash modeling.
"""
import argparse
import pandas as pd
from pyhere import here

def load_eqtl(filename):
    df = pd.read_csv(filename, sep='\t', nrows=100,
                     usecols=["phenotype_id","variant_id","b_gi","b_gi_se"])
    return pd.read_csv(filename, sep='\t', dtype=df.dtypes.to_dict(),
                       usecols=["phenotype_id","variant_id","b_gi","b_gi_se"],
                       compression="gzip")


def extract_eqtls(feature):
    ## Load eQTLs for mashr
    ### Caudate
    cc_file = here(f"prep_eqtl_analysis/caudate/{feature}",
                   "interaction_model/_m",
                   "BrainSEQ_TOPMed.interaction.txt.gz")
    caudate = load_eqtl(cc_file)
    ### DLPFC
    dd_file = here(f"prep_eqtl_analysis/dlpfc/{feature}",
                   "interaction_model/_m",
                   "BrainSEQ_TOPMed.interaction.txt.gz")
    dlpfc = load_eqtl(dd_file)
    ### Hippocampus
    hh_file = here(f"prep_eqtl_analysis/hippocampus/{feature}",
                   "interaction_model/_m",
                   "BrainSEQ_TOPMed.interaction.txt.gz")
    hippo = load_eqtl(hh_file)
    return caudate, dlpfc, hippo


def extract_dataframe(caudate, dlpfc, hippo, variable, label):
    ## Caudate
    dfc = caudate.loc[:, ["phenotype_id","variant_id",variable]]\
                 .rename(columns={variable: "Caudate"})
    # DLPFC
    dfd = dlpfc.loc[:, ["phenotype_id","variant_id",variable]]\
               .rename(columns={variable: "DLPFC"})
    ## Hippocampus
    dfh = hippo.loc[:, ["phenotype_id","variant_id",variable]]\
               .rename(columns={variable: "Hippocampus"})
    df = dfc.merge(dfd, on=["phenotype_id", "variant_id"])\
            .merge(dfh, on=["phenotype_id", "variant_id"])
    df.to_csv(f"{label}_interaction_3regions.txt.gz",
              sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Main
    cc, dd, hh = extract_eqtls(args.feature)
    extract_dataframe(cc, dd, hh, "b_gi", "bhat")
    extract_dataframe(cc, dd, hh, "b_gi_se", "shat")


if __name__=='__main__':
    main()
