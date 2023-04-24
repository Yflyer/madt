#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Yufei,2023
# github.com/Yflyer

from Bio import SeqIO
import os,subprocess,re,pickle,shutil


######################################################################################
### Input parameters of script ###
raw_dir = 'raw_data'
se_dir = 'se_data'
pe_dir = 'pe_data'
up_dir = 'up_data'
raw_se_path = 'raw_se_dict.pkl'
raw_pe_path = 'pe_dict.pkl'

pe_dict = {}
se_dict = {}
up_dict = {}
um_dict = {}

def mkdir(dir):
	if not os.path.exists(dir):
		os.makedirs(dir)
		print ("--- Create a new directory ---")
	else:
		print ("--- The directory has already existed! ---")
mkdir(up_dir)
mkdir(se_dir)
mkdir(pe_dir)

######################################################################################
### index pair_end dataset
for file in os.listdir(pe_dir):
    # check if the file ends with .fastq
    if file.endswith('.fastq'):
        # extract the data ID from the file name
        data_id = re.split(r'[^a-zA-Z0-9]+', file)[0]

        # check if the file name contains '_1' or '_2'
        if '_1.fastq' in file or '_2.fastq' in file:
            # classify the file as potential_pair
            pe_dict.setdefault(data_id, []).append(file)
            pe_dict[data_id].sort()
        elif 'unpaired_1.fastq' in file or 'unpaired_2.fastq' in file:
            # classify the file as up_dict
            up_dict.setdefault(data_id, []).append(file)
            up_dict[data_id].sort()

            # move unpaired files
            if len(up_dict[data_id]) ==2: 
                shutil.move(f"{pe_dir}/{up_dict[data_id][0]}", f"{up_dir}/{up_dict[data_id][0]}")
                shutil.move(f"{pe_dir}/{up_dict[data_id][1]}", f"{up_dir}/{up_dict[data_id][1]}")
                print(f'* Unpaired sequences have been moved to {up_dir} for {data_id}: {up_dict[data_id]}')

#######################################################################################
# PE check
for data_id in list(pe_dict.keys()):
    ####################################################################################### 
    # check duplication
    print('### Now Check pairs from ',data_id)
    fastq1 = SeqIO.parse(f"{pe_dir}/{pe_dict[data_id][0]}", "fastq")
    fastq2 = SeqIO.parse(f"{pe_dir}/{pe_dict[data_id][1]}", "fastq")
    
    # We check the first 20 sequences to determine whether a duplication.
    seq1 = seq2 = ''
    check20=0
    for record1, record2 in zip(fastq1, fastq2):
        # Extract the sequence from the two records
        seq1 += str(record1.seq)
        seq2 += str(record2.seq)
        check20 += 1 
        if check20 == 20: break

    # sort out duplications and merge other pair-end
    if seq1 == seq2:
        # move duplications to se and remove from pair dict
        cmd = f"bbduk.sh in={pe_dir}/{pe_dict[data_id][0]} out={se_dir}/{data_id}_Dup.fasta qtrim=rl trimq=20 minlen=150"
        result = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
        print(f'------Duplication sequences have been moved to {se_dir} of {pe_dict[data_id][0]} and cleaned as Single-End data: {data_id}_SE.fasta')
        os.remove(f"{pe_dir}/{pe_dict[data_id][0]}")
        os.remove(f"{pe_dir}/{pe_dict[data_id][1]}")
        pe_dict.pop(data_id)
    else:
        # Try to merge
        fwd_seq = pe_dict[data_id][0]
        rev_seq = pe_dict[data_id][1]
        cmd = f"bbmerge.sh in1={pe_dir}/{fwd_seq} in2={pe_dir}/{rev_seq} out={pe_dir}/{data_id}_merged.fastq"
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        
        # count sequence for the merged file
        cmd=f"cat {pe_dir}/{data_id}_merged.fastq |  grep \"^@\" | wc -l"
        merge_count = int(subprocess.getoutput(cmd))
        cmd=f"cat {pe_dir}/{fwd_seq} |  grep \"^@\" | wc -l"
        raw_count = int(subprocess.getoutput(cmd))

        # filter the merged files
        if merge_count/raw_count<0.1:
            os.remove(f"{pe_dir}/{fwd_seq}")
            os.remove(f"{pe_dir}/{rev_seq}")
            os.remove(f"{pe_dir}/{data_id}_merged.fastq")
            pe_dict.pop(data_id)

            um_dict[data_id] = [fwd_seq,rev_seq]
            print('* Unmerged ',data_id, ' (Merged rate: ',round(merge_count/raw_count,3),') was removed')
            
        else:
            os.remove(f"{pe_dir}/{fwd_seq}")
            os.remove(f"{pe_dir}/{rev_seq}")

            # PE QC
            cmd = f"bbduk.sh in={pe_dir}/{data_id}_merged.fastq out={pe_dir}/{data_id}_PE.fasta qtrim=rl trimq=20 minlen=150"
            result = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
            os.remove(f"{pe_dir}/{data_id}_merged.fastq")

            pe_dict[data_id]=f"{data_id}_PE.fasta"
            print(f'------Merged pair-end sequences of {pe_dict[data_id]} have been cleaned: {pe_dict[data_id]}')

#######################################################################################
# SE check
with open('raw_se_dict.pkl', 'rb') as f:
    raw_se_dict = pickle.load(f)

for data_id in list(raw_se_dict.keys()):
    cmd = f"bbduk.sh in={raw_dir}/{raw_se_dict[data_id]} out={se_dir}/{data_id}_SE.fasta qtrim=rl trimq=20 minlen=150"
    result = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
    print(f'------Single-end sequences of {raw_se_dict[data_id]} have been cleaned: {data_id}_SE.fasta')

f.close()




