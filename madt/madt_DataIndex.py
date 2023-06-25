#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Yufei,2023
# github.com/Yflyer
import os,subprocess,pickle,re
import argparse
from Bio import SeqIO
import shutil

######################################################################################
# madt_DataIndex.py -i raw_data
### Input parameters of script ###
def get_parser():
    parser = argparse.ArgumentParser(description="""Meta-Amplicon Database Toolkit\n
github.com/Yflyer/ \n
----------
Module: Data index\n
This scripts is to index rawdata downloaded from sra \n
Sort out single-end and pair-end sequences; remove files without enough data\n
""",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i','--Input',required=True,help='Directory of raw data downloaded from SRA')
    parser.add_argument('-s','--SingleDir',default='single_data',help='Directory of outputing single data')
    parser.add_argument('-p','--PairDir',default='pair_data',help='Directory of outputing pair data')
    parser.add_argument('-cl','--CountLimit',default=1000,help='The limitation of remaining data by sequence numbers')
    return parser

def mkdir(dir):
	if not os.path.exists(dir):
		os.makedirs(dir)
		print ("--- Create a new directory ---")
	else:
		print ("--- The directory has already existed! ---")
                
def IndexSe(raw_dir,count_limit):
    raw_se_dict = {}
    for file in os.listdir(raw_dir):
        # check if the file ends with .fastq
        if file.endswith('_pass.fastq'):
            # extract the data ID from the file name
            data_id = re.split(r'[^a-zA-Z0-9]+', file)[0]

            # count sequence for the file
            cmd=f"cat {raw_dir}/{file} |  grep \"^@\" | wc -l"
            read_count = int(subprocess.getoutput(cmd))

            # classify the file as single_end
            if read_count>count_limit:
                raw_se_dict[data_id] = file
                print("### single-end dataset constructed:",data_id,":",raw_se_dict[data_id])
            else:
                print(file,"have not-enough sequence to following work!")
    
    return(raw_se_dict)

def IndexPe(raw_dir,count_limit):
    raw_pe_dict = {}
    for file in os.listdir(raw_dir):
        # check if the file name contains '_1' or '_2'
        if '_1' in file or '_2' in file:
            # extract the data ID from the file name
            data_id = re.split(r'[^a-zA-Z0-9]+', file)[0]

            # count sequence for the file
            cmd=f"cat {raw_dir}/{file} |  grep \"^@\" | wc -l"
            read_count = int(subprocess.getoutput(cmd))

            # select the file for pair_end
            if read_count>count_limit:
                raw_pe_dict.setdefault(data_id, []).append(file)
                raw_pe_dict[data_id].sort()
            else:
                print(file,"have not-enough sequence to following work!")
    
    return(raw_pe_dict)

def SortFile(raw_se_dict,raw_pe_dict,se_dir,pe_dir):
    for data_id in list(raw_pe_dict.keys()):
        # sort out non-pair files as single-end file
        if len(raw_pe_dict[data_id]) == 1:
            raw_se_dict[data_id] = raw_pe_dict[data_id][0]
            raw_pe_dict.pop(data_id)
            print(f'{raw_se_dict[data_id]} was the only used file in {data_id} pair files and regarded as single-end dict')
        else:
            ### output the well paired sequences and unpaired sequences to pe_dir
            fwd_seq, rev_seq = raw_pe_dict[data_id]

            cmd = f"seqkit pair --id-regexp '^(\S+)\.\s?' -1 {raw_dir}/{fwd_seq} -2 {raw_dir}/{rev_seq} -O {pe_dir} -t dna -u"
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            #os.remove(f"{raw_dir}/{fwd_seq}")
            #os.remove(f"{raw_dir}/{rev_seq}")

            print('--- sorting out paired sequences from:',data_id, ' ---')
            print("### Pair-end dataset constructed:",data_id,":",raw_pe_dict[data_id])

    for data_id in list(raw_se_dict.keys()):
        shutil.move(f"{raw_dir}/{raw_se_dict[data_id]}", f"{se_dir}/{raw_se_dict[data_id]}")
        print("### Single-file dataset constructed:",data_id,":",raw_se_dict[data_id])


################################################################
# Run script

parser = get_parser()
args = parser.parse_args()
raw_dir = args.Input
se_dir = args.SingleDir
pe_dir = args.PairDir
count_limit = args.CountLimit

################################################################
# # create a dictionary to store the classification information
mkdir(se_dir)
mkdir(pe_dir)

######################
# index types of library
raw_se_dict = IndexSe(raw_dir,count_limit)
raw_pe_dict = IndexPe(raw_dir,count_limit)

######################################################################################
### sort out sequences
SortFile(raw_se_dict,raw_pe_dict,se_dir,pe_dir)
        
'''
# Import the dictionary from the file using pickle
with open('raw_se_dict.pkl', 'wb') as f:
    pickle.dump(raw_se_dict, f)

with open('raw_pe_dict.pkl', 'wb') as f:
    pickle.dump(raw_pe_dict, f)'
'''
