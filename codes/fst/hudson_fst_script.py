import sys
import os
import argparse
import subprocess
import numpy as np
import scipy
import h5py
import allel
from allel import chunked
import pandas as pd
import zarr
import itertools
print('scikit-allel', allel.__version__)

def fst_analysis(zarrfile, meta)
  callset = zarr.open(zarrfile, mode='r')
  pos_index = allel.ChromPosIndex(callset['variants/CHROM'][:], callset['variants/POS'][:])
  genotype_all = allel.GenotypeChunkedArray(callset['calldata']['GT'])
  df_samples = pd.read_csv(meta)
  geo_r = df_samples['region'].unique()
  for comb in itertools.combinations(geo_r,2): #check combs vs permutations
    pop1 = comb[0]
    pop2 = comb[1]
    n_samples_pop1 = np.count_nonzero(df_samples.region == pop1)
    n_samples_pop2 = np.count_nonzero(df_samples.region == pop2)
    print(pop1, n_samples_pop1, pop2, n_samples_pop2)
    
    subpops = {
        pop1: df_samples[df_samples.region == pop1].index,
        pop2: df_samples[df_samples.region == pop2].index,
    }
    # allele counts
    acs = genotype_all.count_alleles_subpops(subpops)
    acs
    
    acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
    flt = acu.is_segregating() & (acu.max_allele() == 1)
    print('retaining', np.count_nonzero(flt), 'SNPs')
    
    ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
    ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
    genotype = genotype_all.compress(flt, axis=0)
    genotype
    # sample indices for population 1
    pop1_idx = subpops[pop1]
    # sample indices for population 2
    pop2_idx = subpops[pop2]
    num, den = allel.hudson_fst(ac1, ac2)
    snp_fst_hudson = num / den
    snp_fst_hudson
    fst_hudson, se_hudson, vb_hudson, _ = allel.average_hudson_fst(ac1, ac2, blen=10000)
    print('%.04f +/- %.04f (Hudson)' % (fst_hudson, se_hudson))
    
    ## now make some dfs!!
    regions = pop1 + "_vs_" + pop2
    f_name1= regions + ".snp_fst_df.csv"
    f_name2= regions + ".hudson_fst_df.csv"
    snp_fst_df = pd.DataFrame({'CHR': callset['variants/CHROM'][:], 'POS': callset['variants/POS'][:] , 'retained':flt})
    snp_fst_df_flt = snp_fst_df[(snp_fst_df['retained']==True)]
    snp_fst_df_flt['Fst'] = snp_fst_hudson
    snp_fst_df_flt.to_csv(f_name1)
    fst_hudson_df = pd.DataFrame({'Region':regions, 'Hudson_Fst': fst_hudson, 'SE': se_hudson, 'n_retained': [np.count_nonzero(flt)]})
    fst_hudson_df.to_csv(f_name2)

def main(args):
  fst_analysis(args.zarr, args.metadata)

parser = argparse.ArgumentParser(description='Pairwise diversity statistics analysis',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata',help='csv file containing metadata',required=True)
parser.add_argument('--zarr',help='zarr file',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)