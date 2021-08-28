"""
This script combines the chromsomes and corrects for multiple
testing with BH.
"""
import pandas as pd
from glob import iglob
from statsmodels.stats.multitest import multipletests

def get_eFeatures():
    lt = []
    for filename in iglob("*tsv"):
        lt.append(pd.read_csv(filename, sep='\t'))
    df = pd.concat(lt, axis=0)
    _, fdr, _, _ = multipletests(df.BF, method="fdr_bh")
    df["eigenMT_BH"] = fdr
    return df.rename(columns={"p-value": "pval_nominal"})


def main():
    get_eFeatures().to_csv("Brainseq_LIBD.txt.gz", sep='\t', index=False)


if __name__ == '__main__':
    main()
