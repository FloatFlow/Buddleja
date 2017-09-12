"""
Written By Wolfgang Rahfeldt
For John Chau, Evolutionary Analyses in Buddleja
Last edited 4 September 2017

This is a simple script to collect statistics based on designated loci

"""

#Import libraries
import pandas as pd
import numpy as np
import os
import sys
import functools
import argparse

#available arguments
parser = argparse.ArgumentParser()
parser.add_argument('-tf', 
                    '--txtfile', 
                    type=str, 
                    help='Txt file with list of loci used for analyses. Single column of loci names, no header or indices.')
parser.add_argument('-d', 
                    '--directory', 
                    type=str, 
                    help= 'Local directory containing fasta files. Fasta files must be .FNA format.')
parser.add_argument('-q', 
                    '--query', 
                    type=str, 
                    default='noquery', 
                    help= 'Loci to query. Can be a single locus or multiple (e.g. "CAL LFY WAXY"). For multiple queries, each query must be seperated by a space.')
parser.add_argument('-qt', 
                    '--querytype',
                    type=str, 
                    default='notype',
                    help= 'Can be "enrichment" or "depletion". "enrichment" will select loci whose names start with query. "depletion" will select loci whose names do NOT start with query. Default None. ')

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()


#read in target loci
loci_df = pd.read_csv(args.txtfile, sep=' ', header=None)
loci_df.columns = ['loci']

#real jankily parse a fasta file
def fasta_parser(locus):
    t_add = os.path.join(args.directory, str(locus + '.FNA'))

    #read in fasta file
    test_lo = pd.read_csv(t_add, sep='>', header=None)
    test_lo.columns = ['sequence', 'species']

    #parse fasta file
    #get indices of nan values
    inds = list(pd.isnull(test_lo[['sequence']]).any(1).nonzero()[0])
    inds.append(len(test_lo) + 1)
    test_lst = []

    #get each sequence and join fragments
    for i in range(len(inds)-1):
        test_lst.append(test_lo['sequence'].tolist()[inds[i]:inds[i+1]])
    nan_strip = [[x for x in sublist if str(x) != 'nan'] for sublist in test_lst]
    nan_strip = [''.join(sublst) for sublst in nan_strip]

    #construct final df
    wdf = test_lo[['species']].dropna()
    wdf[str(locus)] = nan_strip
    return wdf

#store each locus in a df, in a dictionary of dfs
df_dict = {}
for filename in os.listdir(args.directory):
    for locus in loci_df['loci'].values:
        if filename in str(locus + '.FNA'):
            df_dict[locus] = fasta_parser(locus)


#conditional definitions from query inputs
if args.query != 'noquery':
    query_dict = {}
    for q in args.query.split():
        for k, v in df_dict.items():
            if args.querytype == 'enrichment':
                if k.startswith(q):
                    query_dict[k] = v
            if args.querytype == 'depletion':
                if not k.startswith(q):
                    query_dict[k] = v

if args.query == 'noquery':
    query_dict = df_dict


#create list of dfs via dictionary values, then reduce to merge
dval_lst = list(query_dict.values())
df_final = functools.reduce(lambda left,right: pd.merge(left,right,on='species'), dval_lst)

#Slice sequence data for calculations
ival_df = df_final.ix[:, 1:]
def calc_column(col):
    return col.str.len()
val_df = ival_df.apply(calc_column)
calc_df = df_final[['species']]
pd.options.mode.chained_assignment = None
calc_df['mean'] = val_df.mean(axis=1)
calc_df['max'] = val_df.max(axis=1)
calc_df['min'] = val_df.min(axis=1)
calc_df['sum'] = val_df.sum(axis=1)
calc_df['median'] = val_df.median(axis=1)

#Write to file
destination_csv = os.path.join(args.directory, 'speciesstats_{}_{}.csv'.format(str(args.query), str(args.querytype)))
calc_df.to_csv(destination_csv, index=None)

print('Calculations Complete')
print('Data written to %s' % destination_csv)