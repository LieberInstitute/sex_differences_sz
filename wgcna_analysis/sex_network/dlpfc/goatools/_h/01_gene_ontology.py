"""
Preform gene term enrichment analysis with GOATools.
"""
import session_info
import pandas as pd
import collections as cx
from pybiomart import Dataset
from functools import lru_cache
from statsmodels.stats.multitest import multipletests

# GO analysis
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

@lru_cache()
def get_wgcna_modules():
    return pd.read_csv("../../_m/modules.csv")\
             .rename(columns={"Unnamed: 0": "feature_id"})


@lru_cache()
def get_database():
    dataset = Dataset(name="hsapiens_gene_ensembl",
                      host="http://www.ensembl.org",
                      use_cache=True)
    db = dataset.query(attributes=["ensembl_gene_id",
                                   "external_gene_name",
                                   "entrezgene_id"],
                       use_attr_names=True).dropna(subset=['entrezgene_id'])
    return db


@lru_cache()
def convert2entrez():
    bg = get_wgcna_modules()
    gnames = bg.feature_id.str.split("|", expand=True)\
                              .rename(columns={0: "gene_name",
                                               1: "gene_id"})
    gnames["ensemblID"] = gnames.gene_id.str.replace("\\..*", "",
                                                     regex=True)
    return pd.concat([gnames, bg], axis=1)\
             .set_index("feature_id")\
             .merge(get_database(), left_on='ensemblID',
                    right_on='ensembl_gene_id')


@lru_cache()
def get_module_genes(mod):
    return convert2entrez()[(convert2entrez()["module"] == mod)].copy()


def obo_annotation(alpha=0.05):
    # database annotation
    bg = convert2entrez()
    fn_obo = download_go_basic_obo()
    fn_gene2go = download_ncbi_associations() # must be gunzip to work
    obodag = GODag(fn_obo) # downloads most up-to-date
    anno_hs = Gene2GoReader(fn_gene2go, taxids=[9606])
    # get associations
    ns2assoc = anno_hs.get_ns2assc()
    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated human genes".format(NS=nspc, N=len(id2gos)))
    goeaobj = GOEnrichmentStudyNS(
        bg['entrezgene_id'], # List of human genes with entrez IDs
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = alpha, # default significance cut-off
        methods = ['fdr_bh'])
    return goeaobj


def run_goea(mod):
    df = get_module_genes(mod)
    label = mod
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
    goeaobj.wr_xlsx(f"GO_analysis_mash_{label}.xlsx", goea_results_sig)
    goeaobj.wr_txt(f"GO_analysis_mash_{label}.txt", goea_results_sig)


def main():
    # Run GO enrichment with GOATOOLS
    for mod in get_wgcna_modules().module.unique():
        run_goea(mod)
    ## Session info
    session_info.show()


if __name__ == '__main__':
    main()
