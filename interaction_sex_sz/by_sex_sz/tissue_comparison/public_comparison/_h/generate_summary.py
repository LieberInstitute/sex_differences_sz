"""
This script summarized DE results.
"""
import pandas as pd
import session_info
from pyhere import here
from functools import lru_cache

@lru_cache()
def get_specific_deg():
    fn = "../../summary_table/_m/"+\
        "differential_expression_schizophrenia_by_sex_4features.sig.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == "Gene")].copy()


@lru_cache()
def get_deg():
    fn = "../../summary_table/_m/"+\
        "differential_expression_schizophrenia_by_sex_4features.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == "Gene") & 
              (df["adj.P.Val"] < 0.05)].copy()


@lru_cache()
def get_qin():
    fn = here("input/public_results/_m/qin/qin_results_probesets.csv")
    df = pd.read_csv(fn)
    df["Symbol"] = df.loc[:, 'Gene symbol ']\
                     .str.replace(" ", "", regex=True)
    df["Direction"] = df.loc[:, "Fold difference "]\
                        .str.replace(" ", "", regex=True)
    return df


def print_summary(tissue):
    df1 = get_specific_deg()[(get_specific_deg()["Tissue"] == tissue)].copy()
    df2 = get_deg()[(get_deg()["Tissue"] == tissue)].copy()
    qin = get_qin().loc[:, ["Symbol", "Direction"]]
    tot = set(qin.Symbol)
    if tissue == "Caudate":
        w_mode = "w"
    else:
        w_mode = "a"
    statement = f"{tissue} overlap with Qin et al."
    with open("overlap_summary.log", mode=w_mode) as f:
        print(statement, file=f)
        for sex in ["Female", "Male"]:
            print(sex, file=f)
            specific = df1.loc[(df1["Sex"] == sex),
                               ["Symbol", "Direction"]].copy()
            xx = set(qin.Symbol) & set(specific.Symbol)
            print(f"\nSpecific:\t{len(xx)} "+\
                  f"({len(xx) / len(tot):.1%};{list(xx)})", file=f)
            print(pd.merge(qin, specific, on="Symbol",
                           suffixes=["_qin", "_bs"]), file=f)
            general  = df2.loc[(df2["Sex"] == sex),
                               ["Symbol", "Direction"]].copy()
            yy = set(qin.Symbol) & set(general.Symbol)
            print(f"\nGeneral:\t{len(yy)} "+\
                  f"({len(yy) / len(tot):.1%};{list(yy)})", file=f)
            print(pd.merge(qin, general, on="Symbol",
                           suffixes=["_qin", "_bs"]), file=f)


def main():
    # Summarize overlap with Qin et al
    for region in ["Caudate", "DLPFC", "Hippocampus"]:
        print_summary(region)
        
    # Session infomation
    session_info.show()


if __name__ == '__main__':
    main()
