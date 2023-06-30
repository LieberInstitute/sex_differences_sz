#!/usr/bin/env python
# Author: Francois Aguet
## Edited by Kynon J Benjamin (05/05/2021)

import gzip, os
import subprocess
import rnaseqnorm
import numpy as np
import pandas as pd
import argparse, feather
import scipy.stats as stats


def prepare_bed(df, bed_template_df, chr_subset=None):
    bed_df = pd.merge(bed_template_df, df, left_on="gene_id", right_index=True)
    # sort by start position
    bed_df = bed_df.groupby('chr', sort=False, group_keys=False)\
                   .apply(lambda x: x.sort_values('start'))
    if chr_subset is not None:
        # subset chrs from VCF
        bed_df = bed_df[bed_df.chr.isin(chr_subset)]
    return bed_df


def write_bed(bed_df, output_name):
    """
    Write DataFrame to BED format
    """
    assert bed_df.columns[0]=='chr' and bed_df.columns[1]=='start' and bed_df.columns[2]=='end'
    # header must be commented in BED format
    header = bed_df.columns.values.copy()
    header[0] = '#chr'
    bed_df.to_csv(output_name, sep='\t', index=False, header=header)
    subprocess.check_call(f'bgzip -f {output_name}', shell=True,
                          executable='/bin/bash')
    subprocess.check_call(f'tabix -f {output_name}.gz', shell=True,
                          executable='/bin/bash')


def read_gct(gct_file, sample_ids=None, dtype=None):
    """
    Load GCT as DataFrame. The first two columns must be 'Name' and 'Description'.
    """
    if sample_ids is not None:
        sample_ids = ['Name']+list(sample_ids)
    if gct_file.endswith('.gct.gz') or gct_file.endswith('.gct'):
        if dtype is not None:
            with gzip.open(gct_file, 'rt') as gct:
                gct.readline()
                gct.readline()
                sample_ids = gct.readline().strip().split()
            dtypes = {i:dtype for i in sample_ids[2:]}
            dtypes['Name'] = str
            dtypes['Description'] = str
            df = pd.read_csv(gct_file, sep='\t', skiprows=2, usecols=sample_ids,
                             index_col=0, dtype=dtypes)
        else:
            df = pd.read_csv(gct_file, sep='\t', skiprows=2, usecols=sample_ids,
                             index_col=0)
    elif gct_file.endswith('.parquet'):
        df = pd.read_parquet(gct_file, columns=sample_ids)
    elif gct_file.endswith('.ft'):  # feather format
        df = feather.read_dataframe(gct_file, columns=sample_ids)
        df = df.set_index('Name')
    else:
        raise ValueError('Unsupported input format.')
    df.index.name = 'gene_id'
    if 'Description' in df.columns:
        df = df.drop('Description', axis=1)
    return df


def prepare_expression(counts_df, tpm_df, vcf_lookup_s,
                       sample_frac_threshold=0.2, count_threshold=6,
                       tpm_threshold=0.1, mode='tmm', sort=True):
    """
    Genes are thresholded based on the following expression rules:
      TPM >= tpm_threshold in >= sample_frac_threshold*samples
      read counts >= count_threshold in sample_frac_threshold*samples

    vcf_lookup: lookup table mapping sample IDs to VCF IDs

    Between-sample normalization modes:
      tmm: TMM from edgeR
      qn:  quantile normalization
    """
    if sort:
        ix = np.intersect1d(counts_df.columns, vcf_lookup_s.index)
        tpm_df = tpm_df[ix]
        counts_df = counts_df[ix]
    ns = tpm_df.shape[1]
    # expression thresholds
    mask = (
        (np.sum(tpm_df>=tpm_threshold,axis=1)>=sample_frac_threshold*ns) &
        (np.sum(counts_df>=count_threshold,axis=1)>=sample_frac_threshold*ns)
    ).values
    # apply normalization
    if mode.lower()=='tmm':
        tmm_counts_df = rnaseqnorm.edgeR_cpm(counts_df, normalized_lib_sizes=True)
        norm_df = rnaseqnorm.inverse_normal_transform(tmm_counts_df[mask])
    elif mode.lower()=='qn':
        qn_df = rnaseqnorm.normalize_quantiles(tpm_df.loc[mask])
        norm_df = rnaseqnorm.inverse_normal_transform(qn_df)
    else:
        raise ValueError(f'Unsupported mode {mode}')
    return norm_df


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Generate normalized expression BED files for eQTL analyses')
    parser.add_argument('tpm_gct',
                        help='GCT file with expression in normalized units, e.g., TPM or FPKM')
    parser.add_argument('counts_gct', help='GCT file with read counts')
    parser.add_argument('sample_participant_lookup',
                        help='Lookup table linking samples to participants')
    parser.add_argument('vcf_chr_list', help='List of chromosomes in VCF')
    parser.add_argument('prefix', help='Prefix for output file names')
    parser.add_argument('-o', '--output_dir', default='.',
                        help='Output directory')
    parser.add_argument('--sample_id_list', default=None,
                        help='File listing sample IDs to include')
    parser.add_argument('--convert_tpm', action='store_true',
                        help='Convert to TPM (in case input is in RPKM/FPKM)')
    parser.add_argument('--legacy_mode', action='store_true',
                        help='Run in legacy mode (generates separate output for PEER factor calculation)')
    parser.add_argument('--tpm_threshold', type=np.double, default=0.1,
                        help='Selects genes with > expression_threshold expression in at least sample_frac_threshold')
    parser.add_argument('--count_threshold', type=np.int32, default=6,
                        help='Selects genes with >= count_threshold reads in at least sample_frac_threshold samples')
    parser.add_argument('--sample_frac_threshold', type=np.double, default=0.2,
                        help='Minimum fraction of samples that must satisfy thresholds')
    parser.add_argument('--normalization_method', default='tmm',
                        help='Normalization method: TMM or quantile normalization (qn)')
    parser.add_argument('--feature', default='gene', help='gene, transcript or exon')
    parser.add_argument('--bed_file', help='this is the bed file annotation')
    args = parser.parse_args()

    print('Loading expression data', flush=True)
    sample_ids = None

    if args.sample_id_list is not None:
        with open(args.sample_id_list) as f:
            sample_ids = f.read().strip().split('\n')
            print(f'  * Loading {len(sample_ids)} samples', flush=True)

    counts_df = read_gct(args.counts_gct, sample_ids)
    tpm_df = read_gct(args.tpm_gct, sample_ids)

    if args.convert_tpm:
        print('  * Converting to TPM', flush=True)
        tpm_df = tpm_df/tpm_df.sum(0)*1e6

    print(f'Normalizing data ({args.normalization_method})', flush=True)
    sample_participant_lookup_s = pd.read_csv(args.sample_participant_lookup,
                                              sep='\t', index_col=0, dtype=str)\
                                    .squeeze("columns")
    norm_df = prepare_expression(counts_df, tpm_df, sample_participant_lookup_s,
                                 sample_frac_threshold=args.sample_frac_threshold,
                                 count_threshold=args.count_threshold,
                                 tpm_threshold=args.tpm_threshold,
                                 mode=args.normalization_method, sort=False)
    print(f'  * {counts_df.shape[0]} genes in input tables.', flush=True)
    print(f'  * {norm_df.shape[0]} genes remain after thresholding.', flush=True)

    # change sample IDs to participant IDs
    norm_df.rename(columns=sample_participant_lookup_s.to_dict(), inplace=True)

    bed_template_df = pd.read_csv(args.bed_file, sep='\t', index_col=0)\
                        .rename(columns={'seqnames':'chr'})\
                        .loc[:, ["chr", "start", "end", "gencode_id"]]\
                        .rename(columns={"gencode_id": "gene_id"})
    bed_template_df.loc[:, "end"] = bed_template_df.start + 1
    print('bed_template_df.shape', bed_template_df.shape, flush=True)
    with open(args.vcf_chr_list) as f:
        chr_list = f.read().strip().split('\n')
    norm_bed_df = prepare_bed(norm_df, bed_template_df, chr_subset=chr_list)
    print(f'  * {norm_bed_df.shape[0]} genes remain after removing contigs absent from VCF.', flush=True)
    print('Writing BED file', flush=True)
    write_bed(norm_bed_df, os.path.join(args.output_dir,
                                        args.prefix+'.expression.bed'))

    if args.legacy_mode:
        # for consistency with v6/v6p pipeline results, write unsorted
        # expression file for PEER factor calculation
        norm_df.to_csv(os.path.join(args.output_dir,
                                    args.prefix+'.expression.txt'), sep='\t')
