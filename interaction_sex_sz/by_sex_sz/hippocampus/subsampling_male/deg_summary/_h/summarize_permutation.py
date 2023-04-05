## Examine permutation based on female sample size
import re
import pandas as pd
from glob import iglob

def load_permutation():
    df = pd.DataFrame()
    for filename in iglob("../../_m/permutation_*"):
        m = re.search("\d+", filename)
        dt = pd.read_csv(filename, sep='\t', index_col=0)
        deg = dt[(dt["adj.P.Val"] < 0.05)].copy()
        deg["Permutation"] = m.group(0)
        df = pd.concat([df, deg], axis=0)
    df.to_csv("permutations.csv")
    return df


def get_female_deg():
    fn = "../../../female_analysis/_m/genes/diffExpr_szVctl_FDR05.txt"
    return pd.read_csv(fn, sep='\t', index_col=0)


def summarize_permutation(df):
    xx = df.groupby("Permutation").size()\
           .reset_index().rename(columns={0:"DEGs"})\
           .merge(pd.DataFrame({"Permutation":
                                [str(x).zfill(2) for x in range(1,11)]}), 
                  on="Permutation", how="outer")\
           .fillna(0).sort_values("Permutation")
    print(f"Median: {xx.DEGs.median()}")
    print(xx.DEGs.describe())


def main():
    df = load_permutation()
    summarize_permutation(df)
    print("There are %d DEGs with females!" % get_female_deg().shape[0])


if __name__ == "__main__":
    main()
