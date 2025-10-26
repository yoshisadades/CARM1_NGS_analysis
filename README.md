# CARM1_NGS_analysis
An NGS analysis package published in xxx
Here is how to use/run these codes.
The project was developed and tested with the following setup:
 - OS: Ubuntu 20.04.3 LTS
 - Anaconda 2022.10
 - Python 3.10.9
 - conda 22.11.1
For reproducibility of each packages used, the same environments cab be used with preprocessing.yml and clustering.yml


1) preprocessing

- Edit the preproc.config
  
    line 1: experiment_name
    line 2: working directory
    line 4: primer1 as a DNA sequnece before peptide encoding region
    line 6: primer2 as a DNA sequence after peptide encoding region
    line 8: amino acid linker sequences used for spacer after randamized region

    line 12: raw data location
    line 16: name pattern of raw data files
    line 19: directory for output

    line 22: permitted error rate for primer sequence identification
    line 23: quality score threshold for primer sequence
    line 24: permitted rate of "unknown" base in primer sequence
    line 25: minimum length of the peptide encoding region
    line 26: maximum length of the peptide encoding region

    line 27: minimum length of the peptide
    line 28: maximum length of the peptide 

    line 32: round on which sorting of the sequence abundance is based.
  
- Activate a virtual environment "preprocessing" (use preprocessing.yml)
- Run preproc_fw.sh

This will create the summary of the sequence and counts sorted based on the round specified on line 28.

 2) clustering
 
- Edit and Run Df_to_fasta.py
    line 5: fasta_folder
    line 6: file_path of tsv
- Activate a virtual environment "clustering" (use clustering.yml)
- Run CD-HIT clustering (eaxmple code is given) 
   >cd-hit -i [input].fasta -o [output].fasta -c 0.65 -sc 1 -d 100 
		*-c score: sequence identity score for clustering
		*-sc 1: sort clusters by decreasing cluster size (if library size is the same, this does not matter)
    *-d: description length limit


  
