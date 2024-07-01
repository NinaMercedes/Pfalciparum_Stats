#! /usr/bin/env python

import sys
import os
import argparse
import subprocess
import itertools
import pandas as pd
import fastq2matrix as fm
import numpy as np


### Functions
def clean_data(metadata):
  meta = pd.read_csv(metadata)
  meta = meta.dropna(subset=['region'])
  meta = meta.dropna(subset=['country'])
  meta['country'] = meta['country'].replace("\s+", "_", regex= True) 
  meta['country'] = meta['country'].replace("'", "_", regex= True) 
  meta['region'] = meta['region'].replace("\s+", "_", regex= True) 
  meta['region'] = meta['region'].replace("'", "_", regex= True) 
  return(meta)


def region_fst_pair_analysis(meta):
  geo_r = meta['region'].unique()
  high_fst_snps =[]
  fst_region =[]
  for comb in itertools.combinations(geo_r,2): #check combs vs permutations
    pop1 = comb[0]
    pop2 = comb[1]
    filename = comb[0] + "_vs_" + comb[1] + ".snp_fst_df.csv"
    region_fst = pd.read_csv(filename,header=0)
    region_fst[region_fst['Fst']<0] = 0
    mean_val = np.mean(region_fst['Fst'])
    fst_80 = region_fst[region_fst['Fst']>0.80]
    fst_80["Population 1"] = pop1
    fst_80["Population 2"] = pop2
    num_80 = len(fst_80)
    fst_region.append(pd.DataFrame({'Population 1': pop1, 'Population 2': pop2, 'Fst>0.80':num_80}, index=[0]))
    high_fst_snps.append(fst_80)
  fst_region_df = pd.concat(fst_region, ignore_index=True)
  high_fst_snps_df = pd.concat(high_fst_snps, ignore_index=True)
  return(fst_region_df, high_fst_snps_df)
  
def region_fst_one_analysis(meta):
  geo_r = meta['region'].unique()
  high_fst_snps =[]
  fst_region =[]  
  for region in geo_r:
    filename = "single/" + region + "_vs_all.snp_fst_df_single.csv"
    region_fst = pd.read_csv(filename,header=0)
    fst_80 = region_fst[region_fst['Fst']>0.85]
    fst_80["Region"] = region
    num_80 = len(fst_80)
    fst_region.append(pd.DataFrame({'Region': region, 'Fst>0.80':num_80, 'Fst>0.80':num_80}, index=[0]))
    high_fst_snps.append(fst_80)
  fst_region_df = pd.concat(fst_region, ignore_index=True)
  high_fst_snps_df = pd.concat(high_fst_snps, ignore_index=True)
  return(fst_region_df, high_fst_snps_df)


#UP TO THIS BIT
def main(args):
  #params = {"metadata": args.metadata}
  meta_clean = clean_data(args.metadata)
  fst_region_df1, high_fst_snps_df1 = region_fst_pair_analysis(meta_clean)
  fst_region_df1.to_csv("fst_region_pairs.csv")
  high_fst_snps_df1.to_csv("high_fst_snps_region_pairs_80.csv")
  fst_region_df2, high_fst_snps_df2 = region_fst_one_analysis(meta_clean)
  fst_region_df2.to_csv("fst_region_single.csv")
  high_fst_snps_df2.to_csv("high_fst_snps_region_single_80.csv")

parser = argparse.ArgumentParser(description='Pairwise diversity statistics analysis',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata',help='csv file containing metadata',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
