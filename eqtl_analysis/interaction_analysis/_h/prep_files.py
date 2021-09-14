import argparse
import pandas as pd


def load_eqtl(filename):
    df = pd.read_csv(filename, sep='\t', nrows=100)
    return pd.read_csv(filename, sep='\t', dtype=df.dtypes.to_dict())


def extract_eqtls(feature):
    ## Load eQTLs for mashr
    ### Caudate
    cc_file = "../../../../prep_eqtl_analysis/caudate/%s/" % feature+\
        "prepare_expression/fastqtl_nominal/_m/Brainseq_LIBD.allpairs.txt.gz"
    caudate = load_eqtl(cc_file)
    ### DLPFC
    dlpfc_file = "../../../../prep_eqtl_analysis/dlpfc/%s/" % feature+\
        "prepare_expression/fastqtl_nominal/_m/Brainseq_LIBD.allpairs.txt.gz"
    dlpfc = load_eqtl(dlpfc_file)
    ### Hippocampus
    hippo_file = "../../../../prep_eqtl_analysis/hippocampus/%s/" % feature+\
        "prepare_expression/fastqtl_nominal/_m/Brainseq_LIBD.allpairs.txt.gz"
    hippo = load_eqtl(hippo_file)
    return caudate, dlpfc, hippo


def extract_dataframe(caudate, dlpfc, hippo, variable, label, feature):
    ## Caudate
    df1 = caudate.loc[:, ["gene_id","variant_id",variable]]\
                 .rename(columns={variable: "Caudate"})
    # DLPFC
    df2 = dlpfc.loc[:, ["gene_id","variant_id",variable]]\
               .rename(columns={variable: "DLPFC"})
    ## Hippocampus
    df3 = hippo.loc[:, ["gene_id","variant_id",variable]]\
               .rename(columns={variable: "Hippocampus"})
    df = df1.merge(df2, on=["gene_id", "variant_id"])\
            .merge(df3, on=["gene_id", "variant_id"])
    df.to_csv("%s/%s_fastqtl_3tissues.tsv" % (feature, label),
              sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Main
    cc1, dd1, hh1 = extract_eqtls(args.feature)
    extract_dataframe(cc1, dd1, hh1, "pval_nominal", "pvalue", args.feature)
    extract_dataframe(cc1, dd1, hh1, "slope", "bhat", args.feature)
    extract_dataframe(cc1, dd1, hh1, "slope_se", "shat", args.feature)


if __name__=='__main__':
    main()
