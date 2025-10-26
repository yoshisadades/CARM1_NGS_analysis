#! /bin/bash
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

#obtain variables from config file (change location and filename if needed!)
source D/preproc.config
#sudo mkdir -p $PREPROC_DIR_LOC/fastqc_reports

conda activate preprocessing
sudo mkdir -m 777 -p $PREPROC_DIR_LOC/adapted_seq_data


#for trimming primers only (ignore adapters, as they will be cut anyway), paired, and with quality
for READ1 in ${SEQ_DATA_FOLDER}/$RAW_FILE_PATTERN; do 
    R1FILENAME=$(basename "$READ1" .fastq)

    cutadapt --max-n $MAX_N -q $QUALITY -e $ERROR_RATE --minimum-length $MIN_READ_LEN --maximum-length $MAX_READ_LEN\
    -a ^$PRIMER1...$PRIMER2  \
    -o $PREPROC_DIR_LOC/adapted_seq_data/Adapted${R1FILENAME}.fastq $READ1 $READ2 \
    --untrimmed-output $PREPROC_DIR_LOC/adapted_seq_data/NOTAdapted${R1FILENAME}.fastq \
    --too-short-output $PREPROC_DIR_LOC/adapted_seq_data/TooShort${R1FILENAME}.fastq \
    --too-long-output $PREPROC_DIR_LOC/adapted_seq_data/TooLong${R1FILENAME}.fastq \
    > $PREPROC_DIR_LOC/adapted_seq_data/${R1FILENAME}.txt
done 

echo "Trimmed forward and reverse primer of the forward read. Discarded reads not meeting quality treshold. "

sudo mkdir -m 777 -p $PREPROC_DIR_LOC/peptide_seq_data
sudo mkdir -m 777 -p $PREPROC_DIR_LOC/peptide_count_data

# translate the DNA sequences to protein sequence and save. Count all protein sequences. 
for file in $PREPROC_DIR_LOC/adapted_seq_data/Adapted*$RAW_FILE_PATTERN; do 
    filename=$(basename "$file" .fastq)
    python3 translate.py $file "$PREPROC_DIR_LOC/peptide_seq_data/${filename}.fasta"

    cat "$PREPROC_DIR_LOC/peptide_seq_data/${filename}.fasta" | grep -x '[A-Z\*]\+' |sort | uniq -c | sort -r | sed -E 's/^ *//; s/ /\t/' > $PREPROC_DIR_LOC/peptide_count_data/Translated${filename}_pepcount.tsv

done

echo "Reads translated to protein sequence. Counted."

# summarize counts of all samples in one file
python3 summarize_counts.py --inputfolder $PREPROC_DIR_LOC/peptide_count_data \
--out $WORKING_DIR/combined_counts_${EXPERIMENT_NAME}_e$ERROR_RATE.tsv --minlen $MIN_PEP_LEN --maxlen $MAX_PEP_LEN --filepattern "TranslatedAdapted*_R1_*" --sort_on $RD

conda deactivate
echo "Summary file created."
