#!/bin/env python

import allel
import numpy as np

#vcf_file = '/psych/ripke/1000Genomes_reference/1KG_Oct14/1000GP_Phase3_sr_0517d/ALL_v5a.20130502.chr22_1KG_0517.impute.vcf.gz'
vcf_file = 'test.vcf'

g = allel.read_vcf(vcf_file)
gt = allel.GenotypeArray(g['calldata/GT'])

X = gt.to_n_alt(fill=-1)

del(g, gt)

X = np.ma.masked_array(X, X == -1)
means = np.mean(X, axis=1, keepdims=True)
sds = np.std(X, axis=1, keepdims=True)

X = (X - means) / sds

np.ma.fix_invalid(X, fill_value = 0.0, copy=False)

R = np.matmul(X, X.T) * (1/X.shape[0])

np.save('test_R', np.asarray(R))

