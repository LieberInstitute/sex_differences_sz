
# Examine enrichment in psychiatric disorders TWAS and DEGs
import pandas as pd
import session_info
from os import environ
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

environ['NUMEXPR_MAX_THREADS'] = '10'

@lru_cache()
def get_eqtl(tissue):
    new_tissue = tissue.replace(" ", "_")
    fn = "../../_m/lfsr.sex_interaction.txt.gz"
    df0 = pd.read_csv(fn, sep='\t', nrows=100,
                      usecols=["gene_id", new_tissue])
    return pd.read_csv(fn, sep='\t', dtype=df0.dtypes.to_dict(),
                       usecols=["gene_id", new_tissue],
                       compression="gzip")


@lru_cache()
def get_de(tissue):
    fn = here("differential_expression/tissue_comparison",
              "summary_table/_m",
              "differential_expression_analysis_4features_sex.txt.gz")
    de_df = pd.read_csv(fn, sep='\t')
    return de_df[(de_df["Type"] == "Gene") &
                 (de_df["Tissue"] == tissue)].copy()


@lru_cache()
def get_background(tissue):
    return list(set(get_de(tissue).gencodeID) |
                set(get_eqtl(tissue).gene_id))


@lru_cache()
def get_eGene(tissue):
    return get_eqtl(tissue)[(get_eqtl(tissue)[tissue] < 0.05)]\
        .drop_duplicates(subset="gene_id")


@lru_cache()
def get_deg(tissue, direction):
    df = get_de(tissue)
    if direction=="all":
        return df[(df["adj.P.Val"] < 0.05)].copy()
    elif direction=="up":
        return df[(df["adj.P.Val"] < 0.05) &
                  (df["logFC"] > 0)].copy()
    else:
        return df[(df["adj.P.Val"] < 0.05) &
                  (df["logFC"] < 0)].copy()


def cal_enrichment(tissue, direction):
    """
    Calculates Fisher's Exact test.
    Inputs: brainr region and DE direction of effect.
    """
    universe = set(get_background(tissue))
    de_set = set(get_deg(tissue, direction).gencodeID);
    eGene_set = set(get_eGene(tissue).gene_id)
    yes_de = universe.intersection(de_set)
    yes_eGene = universe.intersection(eGene_set)
    no_de = universe - de_set
    no_eGene = universe - eGene_set
    m = [[len(yes_de.intersection(yes_eGene)),
          len(no_de.intersection(yes_eGene))],
         [len(yes_de.intersection(no_eGene)),
          len(no_de.intersection(no_eGene))]]
    print(m)
    return fisher_exact(m)


def enrichment_loop():
    dir_dict = {"all": "All", "up": "Upregulated in Males",
                "down": "Upregulated in Females"}
    dt = pd.DataFrame()
    for direction in ["all", "up", "down"]:
        print(direction)
        or_lt = []; pval_lt = []; tissue_lt = []
        for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
            print(tissue)
            oddratio, pvals = cal_enrichment(tissue, direction)
            or_lt.append(oddratio); pval_lt.append(pvals);
            tissue_lt.append(tissue)
        fdr = multipletests(pval_lt, method='fdr_bh')[1]
        dtx = pd.DataFrame({"Tissue": tissue_lt, "OR": or_lt,
                            "P-value": pval_lt, "FDR": fdr,
                            "Direction": dir_dict[direction]})
        dt = pd.concat([dt, dtx], axis=0)
    return dt


def main():
    ## Enrichment analysis
    dt = enrichment_loop()
    ## Save enrichment
    dt.to_csv("eGene_enrichment_fishers.tsv", sep='\t', index=False)
    session_info.show()


if __name__ == '__main__':
    main()
