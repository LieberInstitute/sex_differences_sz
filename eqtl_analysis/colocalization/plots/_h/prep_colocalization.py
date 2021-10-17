"""
This script is used to prep files for plotting with
eQTpLot for eQTL-GWAS colocalization plots.

It is modified from a jupyter-notebook Apua development.
Changes are for module plotting.
"""
import numpy as np
import pandas as pd
import subprocess, argparse, functools
from rpy2.robjects import r, pandas2ri, globalenv

@functools.lru_cache()
def get_gwas():
    gwas_fn = '/ceph/projects/v4_phase3_paper/inputs/sz_gwas/'+\
              'pgc2_clozuk/map_phase3/zscore/_m/libd_hg38_pgc2sz_snps.tsv'
    return pd.read_csv(gwas_fn, sep="\t", dtype={'chrN':str}, index_col=0)


@functools.lru_cache()
def subset_gwas(chrom, pos, window):
    gwas_df = get_gwas().loc[(get_gwas()['chrN'] == chrom) &
                             (get_gwas()['pos'] > pos - window) &
                             (get_gwas()['pos'] < pos + window),
                             ['chrN', 'pos', 'our_snp_id', 'P']]\
                        .rename(columns={'chrN':'CHR', 'pos':'BP', 'our_snp_id':'SNP'})
    ## Flip direction of OR based on alleles matching
    gwas_df['BETA'] = np.log(get_gwas()[["OR"]])
    gwas_df['PHE'] = 'SCZD'
    gwas_df['CHR'] = gwas_df['CHR'].astype(int)
    gwas_df['pgc2_a1_same_as_our_counted'] = get_gwas()[["pgc2_a1_same_as_our_counted"]]
    return gwas_df


@functools.lru_cache()
def get_eqtl(fn, feature):
    cmd = '''
    zcat %s | head -1; zcat %s | awk '$1 == "%s" {print}'
    ''' % (fn,fn,feature)
    #print(cmd)
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) as p:
        df = pd.read_csv(p.stdout, sep='\t')
    return df


@functools.lru_cache()
def annotate_eqtls(fn, feature, tissue):
    df = get_eqtl(fn, feature)
    eqtl_df = pd.DataFrame({'SNP.Id': df['variant_id'],
                            'Gene.Symbol': df['gene_id'],
                            'P.Value': df['pval_nominal'],
                            'NES': df['slope'], 'Tissue': tissue},
                           index=df.index)
    return eqtl_df


def flip_slope_by_allele(row):
    return [-1, 1][bool(row["pgc2_a1_same_as_our_counted"])] * row["NES"]


def get_ld(fn, eqtl_dfx, gwas_dfx, label):
    shared_df = gwas_dfx.merge(eqtl_dfx, left_on='SNP', right_on='SNP.Id')\
                        .sort_values('P', ascending=True)
    shared_df[['SNP.Id']].to_csv('snps_%s.txt' % label, index=None, header=None)
    cmd = '''/ceph/opt/plink-ng/1.9/plink \
                --bfile /ceph/users/jbenja13/github_projects/sex_differences_sz/input/genotypes/subset_by_sex/shared_snps/_m/LIBD_Brain_TopMed \
                --extract snps_%s.txt \
                --keep-fam %s --r2 inter-chr \
                --write-snplist --ld-window-r2 0 \
                --out shared_snps_%s;
            sed -i 's/ \+//; s/ \+/\t/g' shared_snps_%s.ld
      ''' % (label,fn,label,label)
    subprocess.run(cmd, shell=True)
    return pd.read_csv("shared_snps_%s.ld" % label, sep='\t', usecols=[*range(7)])


def get_eqtl_by_genes(sex, tissue, gene):
    fn = '/ceph/users/jbenja13/github_projects/sex_differences_sz/'+\
        'prep_eqtl_analysis/by_sex/%s/%s/' % (tissue, sex) +\
        'prepare_expression/fastqtl_nominal/_m/Brainseq_LIBD.allpairs.txt.gz'
    return annotate_eqtls(fn, gene, tissue)


def merge_gwas(eqtl_df, gwas_df):
    eqtl_df = pd.merge(eqtl_df, gwas_df, left_on="SNP.Id", right_on="SNP", how="left")\
                .drop(["CHR", "SNP", "BP", "P", "BETA", "PHE"], axis=1).fillna(True)
    eqtl_df.loc[:,'NES'] = eqtl_df.apply(flip_slope_by_allele, axis=1)
    return eqtl_df.drop(["pgc2_a1_same_as_our_counted"], axis=1)


def get_ld_by_tissue(eqtl_df, gwas_df, tissue, label, sex):
    fn_fam = "/ceph/users/jbenja13/github_projects/sex_differences_sz/"+\
        "prep_eqtl_analysis/by_sex/%s/%s/_m/keepFam.txt" % (tissue, sex)
    return get_ld(fn_fam, eqtl_df, gwas_df, "%s" % (label))


def plot_coloc(gwas_df, genes_df, perm_pval, ld_df, eqtl_df, tissue):
    pandas2ri.activate()
    globalenv["gwas_df"] = gwas_df
    globalenv["genes_df"] = genes_df
    globalenv["ld_df"] = ld_df
    globalenv["eqtl_df"] = eqtl_df
    globalenv["perm_pval"] = perm_pval
    globalenv["tissue"] = tissue
    r("""
    library(eQTpLot)
    pval = perm_pval$perm_pval[1]
    gene = perm_pval$Gene[1]
    p = eQTpLot(GWAS.df=gwas_df, eQTL.df=eqtl_df, Genes.df=genes_df,
                getplot=FALSE, LD.df=ld_df, LDmin=10, R2min=0.25,
                LDcolor='black', gene=gene, trait='SCZD', gbuild='hg38',
                tissue=tissue, sigpvalue_eQTL=pval, CollapseMethod="min",
                congruence=FALSE)
    """)


def main(args):
    gwas_df = subset_gwas("%d" % args.chrom, args.start, args.window)
    perm_pval = pd.DataFrame({"Gene": [args.feature],
                              "perm_pval": [args.perm_pval]})
    genes_df = pd.DataFrame({'CHR':[args.chrom], 'Start':[args.start],
                             'Stop':[args.end], 'Gene':[args.feature],
                             'Build': ['hg38']})
    eqtl_df = get_eqtl_by_genes(args.sex, args.tissue, genes_df.Gene[0])
    eqtl_df = merge_gwas(eqtl_df, gwas_df)
    gwas_df.drop(["pgc2_a1_same_as_our_counted"], axis=1, inplace=True)
    ld_df = get_ld_by_tissue(eqtl_df, gwas_df, args.tissue, args.feature,
                             args.sex)
    plot_coloc(gwas_df, genes_df, perm_pval, ld_df, eqtl_df, args.tissue)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate colocalization prep files')
    parser.add_argument('--chrom', type=np.int32, help='feature chromosome')
    parser.add_argument("--start", type=np.int32, help="feature start position")
    parser.add_argument("--end", type=np.int32, help="feature end position")
    parser.add_argument("--window", default=2e5, type=np.int32,
                        help="window for GWAS extraction")
    parser.add_argument("--feature", type=str, help="feature name used in eQTL analysis")
    parser.add_argument("--perm_pval", type=np.double, help="feature permutation p-value")
    parser.add_argument("--tissue", default="caudate", type=str,
                        help="brain region")
    parser.add_argument("--sex", default="female", help="sex to analyze")
    args = parser.parse_args()
    main(args)
