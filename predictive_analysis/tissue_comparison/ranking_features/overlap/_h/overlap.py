## Examine the overlap between predictive features across brain regions

import functools
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

@functools.lru_cache()
def get_rank():
    fn = "../../_m/BrainSeq_sex_prediction_median_rank_genes.txt.gz"
    return pd.read_csv(fn, sep='\t')


def extract_predictive():
    cc = get_rank()[(get_rank()["Tissue"] == "Caudate") &
                    (get_rank()["Rank"] <= 48)].copy()
    dd = get_rank()[(get_rank()["Tissue"] == "DLPFC") &
                    (get_rank()["Rank"] <= 48)].copy()
    hh = get_rank()[(get_rank()["Tissue"] == "Hippocampus") &
                    (get_rank()["Rank"] <= 48.5)].copy()
    return cc, dd, hh


def print_overlap():
    cc, dd, hh = extract_predictive()
    overlapping = set(cc.Geneid) & set(dd.Geneid) & set(hh.Geneid)
    total = set(cc.Geneid) | set(dd.Geneid) | set(hh.Geneid)
    print("There are {} ({:.1%}) overlapping genes!"\
          .format(len(overlapping),len(overlapping)/len(total)))
    print(overlapping)


def get_n_save_venn():
    cc, dd, hh = extract_predictive()
    setC = set(cc.Geneid)
    setD = set(dd.Geneid)
    setH = set(hh.Geneid)
    venn3([setC, setD, setH], ("Caudate", "DLPFC", "Hippocampus"))
    plt.savefig("venn_diagram_geneOverlap_predictive.pdf")


def main():
    print_overlap()
    get_n_save_venn()


if __name__ == "__main__":
    main()
