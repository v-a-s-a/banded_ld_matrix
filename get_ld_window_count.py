#!/usr/bin/env python

import pandas as pd
import allel as ska
import argparse as arg
import os
import zarr
import numpy as np

def __main__():

    parser = arg.ArgumentParser()
    parser.add_argument('--chr', dest='chrom')
    args = parser.parse_args()

    # read in extra data
    bed = pd.read_csv('/psych/ripke/vasa/reference_data/ldetect-data/EUR/fourier_ls-chr{}.bed'.format(args.chrom),
        sep = '\s+')
    eur_samples = pd.read_csv('/psych/ripke/1000Genomes_reference/1KG_Oct14/1000GP_Phase3_sr_0517d/integrated_call_samples_v3.20130502.ALL.panel.fam.EUR',
        sep='\t', names=['fid', 'iid', 'mid', 'pid', 'sex', 'pheno'], header=None)

    # read in genotype data
    zarr_path = '/psych/ripke/vasa/reference_data/1000G/loc.ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.zarr'.format(args.chrom)
    callset = zarr.open_group(zarr_path, mode='r')
    pos = ska.SortedIndex(callset['variants/POS'])

    callset_samples = list(callset['samples'][:])
    eur_samples['callset_index'] = [callset_samples.index(s) for s in eur_samples['iid']]

    gt = callset['calldata/GT']
    gt_da = ska.GenotypeDaskArray(gt)

    print('Subsetting to europeans')
    eur_da = gt_da.take(eur_samples['callset_index'].values, axis=1)
    eur_ac = eur_da.count_alleles()

    print('Filtering european singletons and invariants')
    flt = (eur_ac.max_allele() == 1) & (eur_ac[:, :2].min(axis=1) > 1)
    flt_mask = flt.compute()
    flt_da = eur_da.compress(flt_mask, axis=0).compute()

    # update variant index
    pos = pos[flt_mask]

    #import ipdb
    #ipdb.set_trace()

    print('Counting region window sizes: ')
    bed['num_variants'] = np.nan
    for i, region in bed.iterrows():
        print('\t{} of {}'.format(i, bed.shape[0]))
        loc_region = pos.locate_range(region['start'], region['stop'])    
        bed.loc[i, ['num_variants']] = flt_da[loc_region, :, :].n_variants

    bed.to_csv('data/1000G_eur_chr{}_region_variant_counts.tsv'.format(args.chrom), sep='\t')

if __name__ == '__main__':
    __main__()
