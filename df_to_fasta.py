# use classify dataset,and convert to fasta format, for compativility with cd-hit. 
#from pos_neg import *
import pandas as pd

fasta_folder = "/home/rubin/proj_CARM1/L/pep_fasta" #CHANGE ME to folder name where output fasta files should be stored. 
file_path = "/home/rubin/proj_CARM1/combined_counts_CARM1_L_e0.1.tsv" #CHANGE ME to name of the files containing peptides and their counts per round. 
def list_to_fasta(seq_list, filename):
    fastafile = open(filename, "w")


    for i in range(len(seq_list)):
        fastafile.write(">" + f"sequence_[{i+1}]_{seq_list[i]}" + "\n" +seq_list[i] + "\n")

    #close the created file
    fastafile.close()

#fulldf: trim of flag and count again
df_full = pd.read_csv(file_path, sep = "\t")
df_full['peptide'] = df_full['peptide'].str[:-6] # trim of flag
df_full= df_full.groupby(['peptide']).sum().reset_index()
print(df_full)
fullpeplist = df_full['peptide'].tolist()

#write all peptides to a fasta file
list_to_fasta(fullpeplist, f"{fasta_folder}/allseq.fasta")

#Rd4 df to fasta
def eachround_to_fasta(n):
    RoundNum = n
    df_round= df_full[df_full[f"Rd{RoundNum}"] != 0]
    df_round= df_round.groupby(['peptide']).sum()
    df_round = df_round.sort_values(f"Rd{RoundNum}", ascending=False).reset_index()
    print(df_round)
    rdpeplist = df_round['peptide'].tolist()
    list_to_fasta(rdpeplist, f"{fasta_folder}/Rd{RoundNum}seq.fasta")
    return

for i in range (1,6):
    eachround_to_fasta(i)
    


