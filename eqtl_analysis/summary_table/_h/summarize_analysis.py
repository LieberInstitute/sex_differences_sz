"""
This script summarizes eQTL analysis.
"""
import functools
import pandas as pd

config = {
    "genes": "/ceph/projects/v4_phase3_paper/inputs/counts/"+\
    "text_files_counts/_m/caudate/gene_annotation.tsv",
    "transcripts": "/ceph/projects/v4_phase3_paper/inputs/counts/"+\
    "text_files_counts/_m/caudate/tx_annotation.tsv",
    "exons": "/ceph/projects/v4_phase3_paper/inputs/counts/"+\
    "text_files_counts/_m/caudate/exon_annotation.tsv",
    "junctions": "/ceph/projects/v4_phase3_paper/inputs/counts/"+\
    "text_files_counts/_m/caudate/jxn_annotation.tsv"
}

@functools.lru_cache()
def get_eFeatures(path):
    return pd.read_csv(path, sep='\t')


@functools.lru_cache()
def get_sig_eFeatures(feature, tissue):
    path = "../../../prep_eqtl_analysis/%s/%s/" % (tissue, feature)+\
        "prepare_expression/fastqtl_nominal/multiple_corrections/_m/"+\
        "Brainseq_LIBD.txt.gz"
    return get_eFeatures(path)[(get_eFeatures(path)["BF"] < 0.05)]


@functools.lru_cache()
def annotate_eqtls(feature, tissue):
    annot = pd.read_csv(config[feature], sep='\t').loc[:,["names", "gencodeID"]]
    return pd.merge(get_sig_eFeatures(feature, tissue), annot,
                    left_on="gene_id", right_on="names")\
             .drop(["names"], axis=1)


@functools.lru_cache()
def load_pgc2():
    pgc2_file = '/ceph/projects/v4_phase3_paper/inputs/sz_gwas/'+\
               'pgc2_clozuk/map_phase3/_m/libd_hg38_pgc2sz_snps_p5e_minus8.tsv'
    return pd.read_csv(pgc2_file, sep='\t', low_memory=False)


@functools.lru_cache()
def merge_pgc2_N_eqtl(feature, tissue):
    return load_pgc2().merge(annotate_eqtls(feature, tissue), how='inner',
                             left_on='our_snp_id', right_on='variant_id',
                             suffixes=['_PGC2', '_eqtl'])


def tissue_map(tissue):
    return {"caudate": "Caudate", "dlpfc": "DLPFC",
            "hippocampus": "Hippocampus"}[tissue]


def load_data(tissue, fnc):
    ## Significant eFeatures (eigenMT corrected)
    genes = fnc("genes", tissue)
    genes["Type"] = "Gene"
    trans = fnc("transcripts", tissue)
    trans["Type"] = "Transcript"
    exons = fnc("exons", tissue)
    exons["Type"] = "Exon"
    juncs = fnc("junctions", tissue)
    juncs["Type"] = "Junction"
    return genes, trans, exons, juncs


def print_summary(genes, trans, exons, juncs):
    ## Total eFeatures
    gg = genes.shape[0]
    tt = trans.shape[0]
    ee = exons.shape[0]
    jj = juncs.shape[0]
    print("eFeatures\neGene:\t\t%d\neTranscript:\t%d\neExon:\t\t%d\neJunction:\t%d" %
          (gg, tt, ee, jj))
    gg = len(set(genes['gencodeID']))
    tt = len(set(trans['gencodeID']))
    ee = len(set(exons['gencodeID']))
    jj = len(set(juncs['gencodeID']))
    print("\neGenes\nGene:\t\t%d\nTranscript:\t%d\nExon:\t\t%d\nJunction:\t%d" %
          (gg, tt, ee, jj))


def combine_data(tissue, PGC2=False):
    print("\n"+tissue_map(tissue))
    if PGC2:
        genes, trans, exons, juncs = load_data(tissue, merge_pgc2_N_eqtl)
        print_summary(genes, trans, exons, juncs)
        df = pd.concat([genes, trans, exons, juncs])\
               .sort_values(["Type", "gene_id", "pval_nominal", "P"])\
               .loc[:, ["variant_id", "rsid", "hg38chrc", "gene_id","gencodeID",
                        "Freq.A1", "A1", "slope", "statistic", "OR", "SE", "P",
                        "pval_nominal", "BF", "eigenMT_BH", "TESTS",
                        "pgc2_a1_same_as_our_counted", "is_index_snp", "Type"]]
    else:
        genes, trans, exons, juncs = load_data(tissue, annotate_eqtls)
        print_summary(genes, trans, exons, juncs)
        df = pd.concat([genes, trans, exons, juncs])\
               .sort_values(["Type", "gene_id", "pval_nominal"])\
               .loc[:, ["variant_id", "gene_id","gencodeID","slope","statistic",
                        "pval_nominal", "BF", "eigenMT_BH", "TESTS", "Type"]]
    df["Type"] = df.Type.astype("category")\
                        .cat.reorder_categories(["Gene", "Transcript",
                                                 "Exon", "Junction"])
    df["Tissue"] = tissue_map(tissue)
    return df


def main():
    ## eFeatures
    with open("summary.log", "w") as f:
        dt = pd.DataFrame()
        print("Interacting features (BF < 0.05):", file=f)
        for tissue in ["caudate", "dlpfc", "hippocampus"]:
            dt = pd.concat([dt, combine_data(tissue, False)])
        print("\neignMT (BF < 0.01):", file=f)
        print(dt[(dt["BF"] < 0.01)].groupby(["Tissue", "Type"]).size(), file=f)
        print("\neignMT FDR corrected (q-value < 0.25):", file=f)
        print(dt[(dt["eigenMT_BH"] < 0.25)].groupby(["Tissue", "Type"]).size(),
              file=f)
        dt.to_csv("Brainseq_sex_interacting_4features_3regions.eFeatures.txt.gz",
                  sep='\t', index=False)
        ## Overlapping schizophrenia risk variants
        """
        Overlap with PGC2+CLOZUK would be very rare as we are only look at top
        variants for each features. We will want to also look at fine mapped
        SNPs and colocalization.
        """
        dt = pd.DataFrame()
        print("\nOverlapping with PGC2+CLOZUK (BF < 0.05):", file=f)
        for tissue in ["caudate", "dlpfc", "hippocampus"]:
            dt = pd.concat([dt, combine_data(tissue, True)])
        print("\neignMT (BF < 0.01):", file=f)
        print(dt[(dt["BF"] < 0.01)].groupby(["Tissue", "Type"]).size(), file=f)
        print("\neignMT FDR corrected (q-value < 0.25):", file=f)
        print(dt[(dt["eigenMT_BH"] < 0.25)].groupby(["Tissue", "Type"]).size(),
              file=f)
        dt.to_csv("Brainseq_sex_interacting_4features_3regions_PGC2.eFeatures.txt.gz",
                  sep='\t', index=False)


if __name__ == '__main__':
    main()
