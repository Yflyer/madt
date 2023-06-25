#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Yufei,2023
# github.com/Yflyer

from Bio import SeqIO
import argparse
import os,subprocess,re,shutil


######################################################################################
# madt_LibCheck.py -s single_data -p pair_data
### Input parameters of script ###
def get_parser():
    parser = argparse.ArgumentParser(description="""Meta-Amplicon Database Toolkit\n
github.com/Yflyer/ \n
----------
Module: Check library\n
This scripts is to index rawdata downloaded from sra\n
Sort out single-end and pair-end sequences; remove files without enough data\n
""",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s','--SingleDir',default='single_data',help='Directory of outputing single data')
    parser.add_argument('-p','--PairDir',default='pair_data',help='Directory of outputing pair-end data')
    parser.add_argument('-up','--UnpairDir',default='unpair_data',help='Directory of outputing pair-end data')
    parser.add_argument('-cl','--CountLimit',default=1000,help='The limitation of remaining data')
    parser.add_argument('-mr','--MergeRatio',default=0.2,help='The limitation of remaining data')
    return parser

###################
parser = get_parser()
args = parser.parse_args()
pe_dir = args.PairDir
se_dir = args.SingleDir
up_dir = args.UnpairDir
count_limit = args.CountLimit
merge_ratio = args.MergeRatio

def mkdir(dir):
	if not os.path.exists(dir):
		os.makedirs(dir)
		print ("--- Create a new directory ---")
	else:
		print ("--- The directory has already existed! ---")

pe_dict = {}
se_dict = {}
up_dict = {}
#um_dict = {}

mkdir(up_dir)
mkdir(se_dir)
mkdir(pe_dir)

######################################################################################
### index pair data
for file in os.listdir(pe_dir):
    # check if the file ends with .fastq
    if file.endswith('.fastq'):
        # extract the data ID from the file name
        data_id = re.split(r'[^a-zA-Z0-9]+', file)[0]

        # check if the file name contains '_1' or '_2'
        # UP check
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

######################################################################################
### index single data
for file in os.listdir(se_dir):
    if file.endswith('.fastq'):
        # extract the data ID from the file name
        data_id = re.split(r'[^a-zA-Z0-9]+', file)[0]
        se_dict[data_id] = file

#######################################################################################
# SE check
for data_id in list(se_dict.keys()):
    cmd = f"bbduk.sh in={se_dir}/{se_dict[data_id]} out={se_dir}/{data_id}_SE.fasta qtrim=rl trimq=20 minlen=150"
    result = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
    os.remove(f"{se_dir}/{se_dict[data_id]}")
    print(f'------Single-end sequences of {se_dict[data_id]} have been cleaned: {data_id}_SE.fasta')

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
        print(f'------Duplication sequences have been moved to {se_dir} of {pe_dict[data_id][0]} and cleaned as Single-End data: {data_id}_Dup.fasta')
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
        if merge_count/raw_count<merge_ratio or merge_count < count_limit:
            os.remove(f"{pe_dir}/{fwd_seq}")
            os.remove(f"{pe_dir}/{rev_seq}")
            os.remove(f"{pe_dir}/{data_id}_merged.fastq")
            pe_dict.pop(data_id)

            #um_dict[data_id] = [fwd_seq,rev_seq]
            print('* Unmerged ',data_id, ' (Merged rate: ',round(merge_count/raw_count,3),'; reads count:',merge_count,') was removed')
            
        else:
            os.remove(f"{pe_dir}/{fwd_seq}")
            os.remove(f"{pe_dir}/{rev_seq}")

            # PE QC
            cmd = f"bbduk.sh in={pe_dir}/{data_id}_merged.fastq out={pe_dir}/{data_id}_PE.fasta qtrim=rl trimq=20 minlen=150"
            result = subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
            os.remove(f"{pe_dir}/{data_id}_merged.fastq")

            pe_dict[data_id]=f"{data_id}_PE.fasta"
            print(f'------Merged pair-end sequences of {pe_dict[data_id]} have been cleaned: {pe_dict[data_id]}')







