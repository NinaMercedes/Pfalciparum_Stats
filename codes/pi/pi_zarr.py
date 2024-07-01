import numpy as np
import scipy
import h5py
import allel
from allel import chunked
import pandas as pd
import zarr
import itertools
import sys
import os
import argparse
import subprocess
print('scikit-allel', allel.__version__)

def nucleotide_diversity(zarrfile, metadata)
  callset = zarr.open(zarrfile, mode='r')
  pos_index = allel.ChromPosIndex(callset['variants/CHROM'][:], callset['variants/POS'][:])
  genotype_all = allel.GenotypeChunkedArray(callset['calldata']['GT'])
  df_samples = pd.read_csv(metadata)
  df_samples = df_samples.dropna(subset=['region'])
  df_samples = df_samples.dropna(subset=['country'])
  geo_r = df_samples['region'].unique()
  geo_c = df_samples['country'].unique()
  for pop in geo_r: #check combs vs permutations
    print(pop)
    subpops = {
    pop: df_samples[df_samples.region == pop].index
    }
    # allele counts
    acs = genotype_all.count_alleles_subpops(subpops)
    acu = allel.AlleleCountsArray(acs[pop])
    pi = allel.mean_pairwise_difference(acu)
    pi_df = pd.DataFrame({'pi':pi})
    f_name = pop + "_pi.csv"
    pi_df.to_csv(f_name)   
  for pop in geo_c: #check combs vs permutations
    print(pop)
    subpops = {
    pop: df_samples[df_samples.country == pop].index
    }
    # allele counts
    acs = genotype_all.count_alleles_subpops(subpops)
    acu = allel.AlleleCountsArray(acs[pop])
    pi = allel.mean_pairwise_difference(acu)
    pi_df = pd.DataFrame({'pi':pi})
    f_name = pop + "_pi.csv"
    pi_df.to_csv(f_name)

def main(args):
  nucleotide_diversity(args.zarr, args.metadata)

parser = argparse.ArgumentParser(description='Pairwise diversity statistics analysis',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata',help='csv file containing metadata',required=True)
parser.add_argument('--zarr',help='zarr file',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)