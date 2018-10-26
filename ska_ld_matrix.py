#!/bin/env python

import allel
import pandas as pd
import numpy as np
import pickle as cp

vcf_file = '/psych/ripke/1000Genomes_reference/1KG_Oct14/1000GP_Phase3_sr_0517d/ALL_v5a.20130502.chr22_1KG_0517.impute.vcf.gz'
#vcf_file = 'data/test.vcf'
bed_file = 'data/fourier_ld-chr22.bed'

bed = pd.read_csv(bed_file, sep='\t')

g = allel.read_vcf(vcf_file)
gt = allel.GenotypeArray(g['calldata/GT'])

X = gt.to_n_alt(fill=-1)

R = allel.windowed_r_squared(g['variants/POS'], X, windows=bed.loc[:, [' start ', ' stop']])

np.dump('data/1000G_eur_R2.pickle', R)

