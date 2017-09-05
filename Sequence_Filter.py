"""
Written by Wolfgang Rahfeldt
For John Chau, Buddleja Project
Written 12 April 2017
Last Edited 4 September 2017

Two (2) required arguments:
    --directory
        folder with output folders for each sample from read_first.py from HybPiper
    --tabfile 
        get_seq_lengths.py output from HybPiper

Three (3) optional arguments:
    --function
        This argument can be either 'paralogs' , 'zeroes', or ‘both’. 
        Default it is set to ‘both’, which will remove genes according to both of the below filters.
        ‘paralogs’ will generate .csv files with a list of loci for Genes with No Paralogs #(Genes_NoParalogs.csv), 
        a list of true paralogs across all samples (Total_Paralogs_Corrected.csv), and a list of true paralogs for each sample inside each 
        sample folder (Paralogs_Corrected.csv)
        'zeroes' will generate a .csv file in the main directory for genes that have less than the threshold of missing data (Genes_Have_Data.csv)

    --missingdata (for zeroes option)
        By default, missingdata=0. This is the number of times a gene can have missing data #before it is removed from the list of candidates

    --contig_cutoff (for paralogs option)
        By default, --contig_cutoff=1. This is the threshold of contigs for which a paralog is considered 'false'. 

On naming scheme of files:
    All files will have two numbers appended to the end (e.g. Genes_NoParalogs_0_1.csv). 
    The first number indicates the value of --missingdata, 
    The second number indicates the value of --contig_cutoff.

Examples: 
    Run the script with defaults:
        $ python Sequence_Filter.py --directory c:/folder --tabfile c:/tab.txt 
    Run script using the 'zeroes' setting, allowing a gene to have missing data for 2 taxa.
        $ python Sequence_Filter.py --directory c:/folder --tabfile c:/tab.txt --function zeroes --missingdata 2

"""

#imports
import pandas as pd
import sys
import os
import subprocess
import itertools
import argparse
from collections import Counter


#available arguments
parser = argparse.ArgumentParser()
parser.add_argument('--tabfile', 
                    type=str, 
                    help='Output from HybPiper, 2D tabfile of gene length and taxa')
parser.add_argument('--directory', 
                    type=str, 
                    help= 'Directory for location of sequencing files')
parser.add_argument('--function', 
                    nargs='?', 
                    const=1, 
                    type=str, 
                    default='both', 
                    help='Determines output. "Zeroes" will remove any missing data, "Paralogs" will remove paralogs, and "Both" will do both. Default = "Both"')
parser.add_argument('--missingdata', 
                    nargs='?', 
                    const=1, 
                    type=int, 
                    default=0, 
                    help= 'Determines how many taxa can be missing data before gene is removed from list of gene candidates. Default = 0')
parser.add_argument('--paralog_cutoff', 
                    nargs='?', 
                    const=1, 
                    type=int, 
                    default=1, 
                    help = 'Determines minimum number taxa that can have paralogs before a particular gene before it is removed from list of gene candidates. Default = 1')
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()


#Arguments conversions
masterlist_address = args.tabfile
folder_address = args.directory
zero_threshold = args.missingdata
function = args.function
nallowed_paralogs = args.paralog_cutoff


#write a list to csv
def write_list_to_csv(loci_list, destination, new_filename):
    loci_df = pd.DataFrame(loci_list, index=None)
    loci_destination = os.path.join(destination, str(new_filename + '_' + str(zero_threshold) + '_' + str(nallowed_paralogs) + '.csv'))
    loci_df.to_csv(loci_destination, index=False)

#Read in pre-generated loci list
masterlist = pd.read_csv(args.tabfile, delim_whitespace=True, header=0)
totalgenes = masterlist.columns[1:].tolist()

#Count number of times a gene is missing data for a taxon
count_zeroes = masterlist[masterlist == 0].count()
candidate_loci = count_zeroes.loc[count_zeroes <= zero_threshold].index.tolist()[1:]
not_candidate_loci = count_zeroes.loc[count_zeroes > zero_threshold].index.tolist()[1:]
if function == 'zeroes' or function == 'both':
    write_list_to_csv(candidate_loci, folder_address, 'Genes_NoMissingData')
    write_list_to_csv(not_candidate_loci, folder_address, 'Genes_MissingData')
    print('Loci without Missing Data: %s') % len(candidate_loci)
    print('Loci with Missing Data: %s') % len(not_candidate_loci)

if function == 'paralogs' or function == 'both':
    species_list = masterlist['Species'].values[1:]

    #get our full true list of paralogs across all species
    full_paralog_list = []
    for species in species_list:

        #define the local directory and retrieve the pregenerated list of potential paralogs
        species_address = os.path.join(folder_address, species)
        paralog_warn = os.path.join(species_address, 'genes_with_paralog_warnings.txt')
        paralog_list = pd.read_csv(paralog_warn, header=None).ix[:,0].tolist()

        #find genes that were flagged as paralogs but only have 1 contig
        false_paralogs = []
        for f in paralog_list:
            if str(f) in totalgenes:
                contig_address = os.path.join(species_address, f, str(f + '_contigs' + '.fasta'))
                contig_count = open(contig_address, 'r').read().count('>')
                if contig_count == 1:
                    false_paralogs.append(str(f))

        #remove false paralogs from the loci that were flagged as paralogs, append to generate full paralog list
        true_paralogs = [locus for locus in paralog_list if locus not in false_paralogs]
        full_paralog_list.append(true_paralogs)

        #write true/false paralogs to local taxon folder
        write_list_to_csv(true_paralogs, species_address, 'Paralogs_Corrected')
        write_list_to_csv(false_paralogs, species_address, 'False_Paralogs')

    #generate final list of true paralogs, write to file
    paralog_dict = Counter(itertools.chain.from_iterable(full_paralog_list))
    paralog_count_df = pd.DataFrame.from_dict(paralog_dict, orient='index').reset_index()
    paralog_count_df.columns = ['gene', 'count']
    nparalog_filtered = paralog_count_df.loc[paralog_count_df['count'] <= nallowed_paralogs]['gene'].tolist()
    write_list_to_csv(nparalog_filtered, folder_address, 'Total_Paralogs_Corrected')

    #generate final list of Loci without paralogs or missing data, write to file
    if function == 'both':
        final_candidate_loci = [locus for locus in candidate_loci if locus not in nparalog_filtered]
        write_list_to_csv(final_candidate_loci, folder_address, 'Genes_NoMissingData_NoParalogs')
        print('Length of Final Candidate Loci: %s') % len(final_candidate_loci)

    #generate final list of Loci without paralogs or missing data, write to file
    if function == 'paralogs' or function == 'both':
        no_paralogs_loci = [locus for locus in totalgenes if locus not in nparalog_filtered]
        write_list_to_csv(no_paralogs_loci, folder_address, 'Genes_NoParalogs')
        print('Length of Loci without Paralogs: %s') % len(no_paralogs_loci)
