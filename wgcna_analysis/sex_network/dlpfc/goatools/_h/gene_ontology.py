"""
Conducts gene term enrichment analysis with GO database for
each module generated from WGCNA.
"""

import functools
import numpy as np
import pandas as pd
import collections as cx
from pybiomart import Dataset
from gtfparse import read_gtf
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# GO analysis
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

@functools.lru_cache()
def get_gtf_genes_df():
    gtf_df = read_gtf("/ceph/genome/human/gencode25/gtf.CHR/_m/gencode.v25.annotation.gtf")
    return gtf_df[gtf_df["feature"] == "gene"][['gene_id', 'gene_name']]


@functools.lru_cache()
def get_wgcna_modules():
    return pd.read_csv("../../_m/modules.csv", index_col=0)


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


def convert2entrez(mod):
    df = get_wgcna_modules()[(get_wgcna_modules().module) == mod].copy()
    df["ensemblID"] = df.index.str.replace("\\..*", "", regex=True)
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


def run_goea(mod):
    df = convert2entrez(mod)
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
    goeaobj.wr_xlsx("GO_analysis_module_%s.xlsx" % mod, goea_results_sig)
    goeaobj.wr_txt("GO_analysis_module_%s.txt" % mod, goea_results_sig)


def main():
    # Load gene annotation
    gtf = get_gtf_genes_df()
    # Load WGCNA modules and annotate
    wgcna_df = get_wgcna_modules()\
        .merge(gtf, left_index=True, right_on="gene_id", how="left")
    wgcna_df.to_csv("module_annotated.csv", index=False)
    # Run GO enrichment with GOATOOLS
    for mod in get_wgcna_modules().module.unique():
        run_goea(mod)

if __name__ == '__main__':
    main()
