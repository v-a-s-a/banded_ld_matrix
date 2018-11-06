#!/usr/bin/env python

import allel as ska
import numcodecs

vcf_file = 'data/test.vcf.gz'
# vcf_file = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
test_fields = ['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT', 'variants/QUAL', 'variants/FILTER', 'variants/INFO']
# test_fields += ['calldata/GT',
# 'variants/CIEND',
# 'variants/CIPOS',
# 'variants/CS',
# 'variants/END',
# 'variants/IMPRECISE',
# 'variants/MC',
# 'variants/MEINFO',
# 'variants/MEND',
# 'variants/MLEN',
# 'variants/MSTART',
# 'variants/SVLEN',
# 'variants/SVTYPE',
# 'variants/TSD',
# 'variants/AC',
# 'variants/AF',
# 'variants/NS',
# 'variants/AN',
# 'variants/EAS_AF',
# 'variants/EUR_AF',
# 'variants/AFR_AF',
# 'variants/AMR_AF',
# 'variants/SAS_AF',
# 'variants/DP',
# 'variants/AA',
# 'variants/VT',
# 'variants/EX_TARGET',
# 'variants/MULTI_ALLELIC']

# test_fields += ['variants/numalt', 'variants/svlen', 'variants/is_snp']
test_fields += ['variants/numalt','variants/is_snp', 'variants/svlen']


# test_fields = ['variants/*']

ska.vcf_to_zarr(vcf_file, vcf_file.replace('.vcf.gz', '.zarr'),
                  fields=test_fields, alt_number=8, overwrite=True,
                  compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False))
