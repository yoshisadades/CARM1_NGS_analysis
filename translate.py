#Author: Isabel Houtkamp

# Required imports
import argparse
import sys
import Bio
import Bio.SeqIO
from Bio.Seq import Seq

"""Specify path to fastq file with DNA sequence to translate on command line, as well
as the path to the output file to create, where protein sequences will be stored.

EXAMPLE to run from terminal:
python3 translate.py myDNAsequences.fastq myPROTsequences.fastq

Customize the codon table if needed! In the current codon table, start codon 'ATG' is reprogrammed to Y. 
"""

# input and output file are specified on the command line
input_fastq = sys.argv[1]
output_fasta = sys.argv[2]

codon_table = {
"TTT": "F", "TCT": "S", "TAT": "Y", "TGT": "C",
"TTC": "F", "TCC": "S", "TAC": "Y", "TGC": "C",
"TTA": "L", "TCA": "S", "TAA": "*", "TGA": "*",
"TTG": "L", "TCG": "S", "TAG": "*", "TGG": "W",

"CTT": "L", "CCT": "P", "CAT": "H", "CGT": "R",
"CTC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
"CTA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
"CTG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",

"ATT": "I", "ACT": "T", "AAT": "N", "AGT": "S",
"ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
"ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
"ATG": "Y", "ACG": "T", "AAG": "K", "AGG": "R",

"GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
"GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
"GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
"GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G"}

start_codon = 'ATG'
stop_codons = ["TAA", "TAG", "TGA"]

# Function to translate DNA seq to protein seq
def translate(seq, codon_table):
    """"Translate a DNA sequence to a protein sequence. 

    Parameters:

    seq: a string containing a DNA sequence made up of letters A,C,G,T or N   
    codon_table: a dictionary containing a codon table of choiche. 

    """
    peptide = ""
    start = seq.find(start_codon)
    for i in range(start, len(seq),3):
        codon = seq[i:i+3]
        #if codon is of length 3, continue translate
        if len(codon) != 3:
            #print(f"Codon encountered of less then 3 amino acids:{codon}. Translation terminated., frameshift may have occured.")
            break   
        if codon in codon_table:
            peptide+= codon_table[codon]
            if codon in stop_codons:
                break
        else:
            peptide+= "X" #unknown amino acid

    return peptide

def main():

    """ Read each record in the fastq file, apply the translate() function.  
    Write the resulting peptide sequences to a fasta file, presevering the read identifiers. 
    """

    with open(input_fastq) as handle:
        record_list = []
        fastq_seq = Bio.SeqIO.parse(handle, "fastq")
        for record in fastq_seq:
            record.letter_annotations = {}
            record.seq = Seq(translate(record.seq, codon_table))
            record_list.append(record)

        Bio.SeqIO.write(record_list, output_fasta, "fasta")
    handle.close()
    
if __name__ == "__main__":
    main()
 
