#### Examining enrichment in psychiatric disorders (TWAS and DEGs)

import session_info
import pandas as pd
from os import environ
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

environ['NUMEXPR_MAX_THREADS'] = '10'
config = {
    "caud8_file": here("input/neuropsychiatric_results/_m",
                       "BrainSeq_Phase3_Caudate_DifferentialExpression_DxSZ_all.txt.gz"),
    "dlpfc_file": here("input/neuropsychiatric_results/_m",
                       "Brainseq_phase2_qsvs_age17_noHGold_sz_DEG.tsv"),
    "hippo_file": here("input/neuropsychiatric_results/_m",
                       "Brainseq_phase2_qsvs_age17_noHGold_sz_DEG.tsv"),
    'cmc_file': here("input/neuropsychiatric_results/_m",
                     'CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500'+\
                     '_gene-adjustedSVA-differentialExpression-includeAncestry-DxSCZ-DE.tsv'),
    'PW_parkinson_file': here("input/neuropsychiatric_results/_m/Nido2020_PD_pwCohort.csv"),
    'alzheimers_file': here("input/neuropsychiatric_results/_m/MarquesCoelho2021_AD.csv"),
    "libd_bipolar_file": here("input/neuropsychiatric_results/_m",
                              "bipolarControl_deStats_byRegion_"+\
                              "qSVAjoint_withAnnotation_all.csv"),
    'gandal_de_file': here("input/neuropsychiatric_results/_m",
                           "gandal2018_psychENCODE_DE_results.xlsx"),
    'twas_asd_file': here("input/neuropsychiatric_results/_m/psychENCODE_twas_asd.csv"),
    'twas_sz_file': here("input/neuropsychiatric_results/_m/psychENCODE_twas_sz.csv"),
    'twas_bd_file': here("input/neuropsychiatric_results/_m/psychENCODE_twas_bd.csv"),
    'twas_bs_file': here("input/neuropsychiatric_results/_m/TWAS_gene_tissue_summary.csv"),
    'twas_ad_file': here("input/neuropsychiatric_results/_m/Gockley2021_AD.tsv"),
    'twas_bd_libd_file': here("input/neuropsychiatric_results/_m",
                              'twas_2regions_bpd_output_allTested.csv'),
}

@lru_cache()
def get_deg(fn, REGION):
    df = pd.read_csv(fn, sep='\t')
    if REGION.lower() != "cmc":
        df = df[(df["Type"] == "Gene")].copy()
    if (REGION.lower() == "caudate") | (REGION.lower() == "cmc"):
        return df
    else:
        return df[(df["Tissue"] == REGION)].copy()


@lru_cache()
def get_sig_deg(fn, REGION):
    return get_deg(fn, REGION)[(get_deg(fn, REGION)["adj.P.Val"] < 0.05)].copy()


@lru_cache()
def get_bd_deg():
    df = pd.read_csv(config["libd_bipolar_file"], index_col=0)
    df["ensemblID"] = df.gencodeID.str.replace("\..*", "", regex=True)
    return df[(df["Type"] == "Gene")].copy()


@lru_cache()
def get_bd_sig_deg(label):
    return get_bd_deg()[(get_bd_deg()[label] < 0.05)].copy()


@lru_cache()
def get_gandal_deg():
    return pd.read_excel(config["gandal_de_file"], sheet_name="DGE")


@lru_cache()
def get_parkinson_deg():
    return pd.read_csv(config['PW_parkinson_file'])


@lru_cache()
def get_alzheimers_deg(dataset):
    df = pd.read_csv(config['alzheimers_file'], skiprows=1)\
           .drop_duplicates(subset="GENEID")
    df["ensemblID"] = df.GENEID.str.replace("\..*", "", regex=True)
    return df[(df["gene.padj"] < 0.05) & (df["dataset"] == dataset)].copy()


@lru_cache()
def load_eQTL():
    fn = "../../summary_data/_m/BrainSeq_sexGenotypes_genes_3regions.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    df["ensemblID"] = df.gencode_id.str.replace("\\..*", "", regex=True)
    return df.drop_duplicates("ensemblID")


@lru_cache()
def get_eGenes(tissue, direction):
    fn = "../../../feature_summary/_m/BrainSeq_sexGenotypes_4features_3regions.txt.gz"
    df = pd.read_csv(fn, sep="\t")
    df["ensemblID"] = df.gencode_id.str.replace("\\..*", "", regex=True)
    df = df[(df["region"] == tissue) & (df["feature_type"] == "Gene") &
            (df["lfsr"] < 0.05)].copy()
    if direction=="all":
        return df
    elif direction=="up":
        return df[(df["posterior_mean"] > 0)].copy()
    else:
        return df[(df["posterior_mean"] < 0)].copy()


@lru_cache()
def get_bs_twas():
    return pd.read_csv(config["twas_bs_file"])


@lru_cache()
def get_libd_bd_twas():
    return pd.read_csv(config['twas_bd_libd_file'])


@lru_cache()
def get_gandal_twas(fn):
    df = pd.read_csv(fn)
    return df[(df["TWAS.Bonferroni"] < 0.05)].copy()


@lru_cache()
def get_ad_twas():
    df = pd.read_csv(config['twas_ad_file'], sep='\t', header=None,
                     names=["chr", "set", "feature", "start", "end", "na1",
                            "strand", "na2", "geneid", "gene_name"])
    df["ensemblID"] = df.geneid.str.replace("\..*", "", regex=True)
    return df


def fet(a, b):
    """
    Calculates Fisher's Exact test (fet) with sets a and b in universe u.
    Inputs are sets.
    """
    u = set(load_eQTL().ensemblID)
    a = set(a); b = set(b)
    yes_a = u.intersection(a)
    yes_b = u.intersection(b)
    no_a = u - a
    no_b = u - b
    m = [[len(yes_a.intersection(yes_b)), len(no_a.intersection(yes_b))],
         [len(yes_a.intersection(no_b)), len(no_a.intersection(no_b))]]
    return fisher_exact(m)


def load_eGenes(direction="all"):
    caudate = set(get_eGenes("Caudate", direction).ensemblID)
    dlpfc = set(get_eGenes("DLPFC", direction).ensemblID)
    hippocampus = set(get_eGenes("Hippocampus", direction).ensemblID)
    return {"Caudate": caudate, "DLPFC": dlpfc, "Hippocampus": hippocampus}


def get_public_degs():
    # DEGs
    ## BrainSeq SZ
    bs_caudate_degs = set(get_sig_deg(config["caud8_file"], "Caudate").ensemblID)
    bs_dlpfc_degs = set(get_sig_deg(config["dlpfc_file"], "DLPFC").ensemblID)
    bs_hippo_degs = set(get_sig_deg(config["hippo_file"], "Hippocampus").ensemblID)
    ## CommonMind SZ
    cmc_dlpfc_degs = set(get_sig_deg(config["cmc_file"], "CMC").genes)
    ## LIBD Bipolar
    bs_sacc_degs = set(get_bd_sig_deg("adj.P.Val_sACC").ensemblID)
    bs_amyg_degs = set(get_bd_sig_deg("adj.P.Val_Amyg").ensemblID)
    ## PsychENCODE (Gandal)
    psy_sz = set(get_gandal_deg()[(get_gandal_deg()["SCZ.fdr"] < 0.05)].ensembl_gene_id)
    psy_asd = set(get_gandal_deg()[(get_gandal_deg()["ASD.fdr"] < 0.05)].ensembl_gene_id)
    psy_bd = set(get_gandal_deg()[(get_gandal_deg()["BD.fdr"] < 0.05)].ensembl_gene_id)
    ## Parkinson's
    pd_degs = set(get_parkinson_deg()[(get_parkinson_deg()["padj"] < 0.05)].EnsemblID)
    ## Alzheimer's
    ad_mayo = set(get_alzheimers_deg("MAYO").ensemblID)
    ad_rosmap = set(get_alzheimers_deg("ROSMAP").ensemblID)
    ad_msbb_mb36 = set(get_alzheimers_deg("MSBB BM36").ensemblID)
    ad_msbb_mb44 = set(get_alzheimers_deg("MSBB BM44").ensemblID)
    ad_msbb_mb22 = set(get_alzheimers_deg("MSBB BM22").ensemblID)
    ad_msbb_mb10 = set(get_alzheimers_deg("MSBB BM10").ensemblID)
    return {"BS_Caudate_SZ": bs_caudate_degs, "BS_DLPFC_SZ": bs_dlpfc_degs,
            "BS_Hippocampus_SZ": bs_hippo_degs, "CMC_DLPFC_SZ": cmc_dlpfc_degs,
            "BS_sACC_BD": bs_sacc_degs, "BS_Amyg_BD": bs_amyg_degs,
            "PSY_SZ": psy_sz, "PSY_ASD": psy_asd, "PSY_BD": psy_bd,
            "Parkinsons": pd_degs, "MAYO_AD": ad_mayo, "ROSMAP_AD": ad_rosmap,
            "MSBB_MD36_AD": ad_msbb_mb36, "MSBB_MD44_AD": ad_msbb_mb44,
            "MSBB_MD22_AD": ad_msbb_mb22, "MSBB_MD10_AD": ad_msbb_mb10}


def get_public_twas():
    # TWAS
    ## BrainSeq SZ
    bs_caudate_twas = set(get_bs_twas()[(get_bs_twas()["Caudate_FDR"] < 0.05)].Geneid)
    bs_dlpfc_twas = set(get_bs_twas()[(get_bs_twas()["DLPFC_FDR"] < 0.05)].Geneid)
    bs_hippo_twas = set(get_bs_twas()[(get_bs_twas()["HIPPO_FDR"] < 0.05)].Geneid)
    ## LIBD Bipolar
    bs_sacc_twas = set(get_libd_bd_twas()[(get_libd_bd_twas()["Region"] == "sACC") &
                                          (get_libd_bd_twas()["TWAS.fdrP"]<0.05)].ID)
    bs_amyg_twas = set(get_libd_bd_twas()[(get_libd_bd_twas()["Region"] == "Amygdala") &
                                          (get_libd_bd_twas()["TWAS.fdrP"]<0.05)].ID)
    ## PsychENCODE
    psy_asd_twas = set(get_gandal_twas(config["twas_asd_file"]).GeneID)
    psy_sz_twas = set(get_gandal_twas(config["twas_sz_file"]).GeneID)
    psy_bd_twas = set(get_gandal_twas(config["twas_bd_file"]).ID)
    ## Alzheimer's
    ad_twas = set(get_ad_twas().ensemblID)
    return {"BS_Caudate_SZ": bs_caudate_twas, "BS_DLPFC_SZ": bs_dlpfc_twas,
            "BS_Hippocampus_SZ": bs_hippo_twas, "PSY_SZ": psy_sz_twas,
            "PSY_ASD": psy_asd_twas, "PSY_BD": psy_bd_twas,
            "BS_sACC_BD": bs_sacc_twas, "BS_Amyg_BD": bs_amyg_twas,
            "Alzheimers": ad_twas}


def enrichment_loop(comp_list, comp_dict, label):
    dir_dict = {"all": "All", "up": "Up", "down": "Down"}
    dt = pd.DataFrame()
    for direction in ["all", "up", "down"]:
        ## Load eGenes
        egenes_dict = load_eGenes(direction)
        or_lt = []; pval_lt = []; tissue_lt = []; comp_lt = []
        for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
            for comp in comp_list:
                oddratio, pvals = fet(egenes_dict[tissue], comp_dict[comp])
                or_lt.append(oddratio); pval_lt.append(pvals);
                tissue_lt.append(tissue); comp_lt.append(comp)
        fdr = multipletests(pval_lt, method='fdr_bh')[1]
        dtx = pd.DataFrame({"Tissue": tissue_lt, "Dataset": comp_lt,
                            "OR": or_lt, "P-value": pval_lt, "FDR": fdr,
                            "Direction":dir_dict[direction], "Method":label})
        dt = pd.concat([dt, dtx], axis=0)
    return dt


def main():
    ## Load public DEGs
    degs_dict = get_public_degs()
    ## Load public TWAS
    twas_dict = get_public_twas()
    ## Labels
    degs_list = ["BS_Caudate_SZ", "BS_DLPFC_SZ", "BS_Hippocampus_SZ",
                 "BS_sACC_BD", "BS_Amyg_BD", "CMC_DLPFC_SZ", "PSY_SZ",
                 "PSY_ASD", "PSY_BD", "Parkinsons", "MAYO_AD", "ROSMAP_AD",
                 "MSBB_MD36_AD", "MSBB_MD44_AD", "MSBB_MD22_AD", "MSBB_MD10_AD"]
    twas_list = ["BS_Caudate_SZ", "BS_DLPFC_SZ", "BS_Hippocampus_SZ",
                 "BS_sACC_BD", "BS_Amyg_BD", "PSY_SZ", "PSY_ASD", "PSY_BD",
                 "Alzheimers"]
    ## Direction dictionary
    dt_degs = enrichment_loop(degs_list, degs_dict, "DEG")
    dt_twas = enrichment_loop(twas_list, twas_dict, "TWAS")
    ## Save enrichment
    dt_degs.to_csv("clincial_phenotypes_enrichment_analysis_DEGs.tsv",
                   sep='\t', index=False)
    dt_twas.to_csv("clincial_phenotypes_enrichment_analysis_TWAS.tsv",
                   sep='\t', index=False)


if __name__ == '__main__':
    main()
