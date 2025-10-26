# Author: Isabel Houtkamp

import argparse
import sys
import glob
import os
import pandas as pd
import re

# Create argument parser, to read in arguments specified on command line
parser = argparse.ArgumentParser()
parser.add_argument("--inputfolder", type=str, nargs='?', help="folder containing peptide counts per round")
parser.add_argument("--out", type=str, nargs='?', help="name of outputfile")
parser.add_argument("--minlen", type=int, default=0, nargs='?', help="minimum peptide length")
parser.add_argument("--maxlen", type=int, default=100, nargs='?', help="maximum peptide length")
parser.add_argument("--filepattern", type=str, default="TranslatedAdapted*", nargs='?', help="file pattern for count tables of each sequencing round")
parser.add_argument("--sort_on", type=str, default="Rd1", nargs='?', help="Specify which selection round to sort counts on (for example: --sort on Rd4)")
# arguments are stored in args
args = parser.parse_args()

def read_files(inputfolder, filepattern):

    """ Read file in the input folder, matching the specified filename pattern. This should be a   
    Write the resulting peptide sequences to a fasta file, presevering the read identifiers. 
    """    

    filenames = glob.glob(f"{inputfolder}/{filepattern}.tsv")
    df_list = []
    for fname in filenames:
        basename = os.path.basename(fname)
        print(basename)
        round_spec = basename.find("Rd")
        round = basename[round_spec: round_spec + 3]
        data = pd.read_csv(fname, delimiter = "\t", header = None, dtype = {'count': int, 'peptide': str})
        data.columns = [round, "peptide"]
        data = data.set_index('peptide')
        print(data)
        df_list.append(data)
        print(data)
        print(f"The round nr is: {round}" )
        

    return df_list

def join_dfs(df_list, min_len, max_len, out):
        
    """ Joins a list of pandas dataframes containing the peptide sequence counts per sequencing file into a combined
    pandas dataframe. 
    Peptides sequences are filtered for the specified min and max length. 
    Returns the joined dataframe, and writes it to a csv file. 
    """    

    # Join dataframes into one and replace NA values with 0
    joined_data = df_list[0].join(df_list[1:], how='outer')
    joined_data = joined_data.fillna(0.0).astype(int)

    # Drop peptide sequences that are too short or to loo long from the final output
    joined_data = joined_data[joined_data.index.str.len() > min_len]
    joined_data = joined_data[joined_data.index.str.len() < max_len]
    joined_data = joined_data.sort_values(by=args.sort_on, ascending=False)
    joined_data.sort_index(axis=1, inplace=True)

    joined_data.to_csv(out, sep = "\t")
    return joined_data

def main():
    df_list = read_files(args.inputfolder, args.filepattern)
    joined_data = join_dfs(df_list, args.minlen, args.maxlen, args.out)

if __name__ == "__main__":
    main()
 


