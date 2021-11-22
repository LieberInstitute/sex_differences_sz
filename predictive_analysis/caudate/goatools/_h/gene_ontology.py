"""
Conducts gene term enrichment analysis with GO database for
max predictive and redundant subsets from dRFE.
"""

import numpy as np
import pandas as pd
import functools, os
import collections as cx
from pybiomart import Dataset
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# GO analysis
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

os.environ["NUMEXPR_MAX_THREADS"] = '32'

@functools.lru_cache()
def get_summary():
    fn = "../../_m/genes/dRFEtools_10folds.txt"
    return pd.read_csv(fn, sep='\t').median()


@functools.lru_cache()
def get_pred_features():
    fn = "../../_m/genes/rank_features.txt"
    ml_df = pd.read_csv(fn, sep='\t', header=None,
                        names=["Geneid", "Fold", "rank"])\
              .pivot(index="Geneid", columns="Fold")\
              .median(axis=1).reset_index().rename(columns={0: "Rank"})\
              .sort_values("Rank")
    ml_df["Rank"] = ml_df["Rank"].rank()
    return ml_df[(ml_df["Rank"] <= get_summary()["n_features"])]


@functools.lru_cache()
def get_redundant_features():
    fn = "../../_m/genes/rank_features.txt"
    ml_df = pd.read_csv(fn, sep='\t', header=None,
                        names=["Geneid", "Fold", "rank"])\
              .pivot(index="Geneid", columns="Fold")\
              .median(axis=1).reset_index().rename(columns={0: "Rank"})\
              .sort_values("Rank")
    ml_df["Rank"] = ml_df["Rank"].rank()
    return ml_df[(ml_df["Rank"] <= get_summary()["n_redundant"])]


@functools.lru_cache()
def get_database():
    dataset = Dataset(name="hsapiens_gene_ensembl",
                      host="http://www.ensembl.org",
                      use_cache=True)
    db = dataset.query(attributes=["ensembl_gene_id",
                                   "external_gene_name",
                                   "entrezgene_id"],
                       use_attr_names=True).dropna(subset=['entrezgene_id'])
    return db


def convert2entrez(fnc):
    df = fnc()
    df["ensemblID"] = df.Geneid.str.replace("\\..*", "", regex=True)
    return df.merge(get_database(), left_on='ensemblID',
                    right_on='ensembl_gene_id')


def obo_annotation(alpha=0.05):
    # database annotation
    fn_obo = download_go_basic_obo()
    fn_gene2go = download_ncbi_associations() # must be gunzip to work
    obodag = GODag(fn_obo) # downloads most up-to-date
    anno_hs = Gene2GoReader(fn_gene2go, taxids=[9606])
    # get associations
    ns2assoc = anno_hs.get_ns2assc()
    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated human genes".format(NS=nspc, N=len(id2gos)))
    goeaobj = GOEnrichmentStudyNS(
        get_database()['entrezgene_id'], # List of human genes with entrez IDs
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = alpha, # default significance cut-off
        methods = ['fdr_bh'])
    return goeaobj


def run_goea(fnc, label):
    df = convert2entrez(fnc)
    geneids_study = {z[0]:z[1] for z in zip(df['entrezgene_id'],
                                            df['external_gene_name'])}
    goeaobj = obo_annotation()
    goea_results_all = goeaobj.run_study(geneids_study)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    ctr = cx.Counter([r.NS for r in goea_results_sig])
    print('Significant results[{TOTAL}] = {BP} BP + {MF} MF + {CC} CC'.format(
        TOTAL=len(goea_results_sig),
        BP=ctr['BP'],  # biological_process
        MF=ctr['MF'],  # molecular_function
        CC=ctr['CC'])) # cellular_component
    goeaobj.wr_xlsx("GO_analysis_%s.xlsx" % label, goea_results_sig)
    goeaobj.wr_txt("GO_analysis_%s.txt" % label, goea_results_sig)


def main():
    # Run GO enrichment with GOATOOLS
    run_goea(get_pred_features, "predictive")
    run_goea(get_redundant_features, "redundant")

if __name__ == '__main__':
    main()
