#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Yufei,2023
# github.com/Yflyer
import os,subprocess,pickle,re
from Bio import SeqIO

######################################################################################
### Input parameters of script ###
raw_dir = 'raw_data'
se_dir = 'se_data'
pe_dir = 'pe_data'

def mkdir(dir):
	if not os.path.exists(dir):
		os.makedirs(dir)
		print ("--- Create a new directory ---")
	else:
		print ("--- The directory has already existed! ---")
mkdir(pe_dir)

# create a dictionary to store the classification information
raw_se_dict = {}
raw_pe_dict = {}

######################################################################################
# index types of library
for file in os.listdir(raw_dir):

    # check if the file ends with .fastq
    if file.endswith('_pass.fastq'):
        # extract the data ID from the file name
        data_id = re.split(r'[^a-zA-Z0-9]+', file)[0]

        # count sequence for the file
        cmd=f"cat {raw_dir}/{file} |  grep \"^@\" | wc -l"
        read_count = int(subprocess.getoutput(cmd))

        # classify the file as single_end
        if read_count>1000:
            raw_se_dict[data_id] = file
            print("### single-end dataset constructed:",data_id,":",raw_se_dict[data_id])
        else:
            print(file,"have not-enough sequence to following work!")

    # check if the file name contains '_1' or '_2'
    if '_1' in file or '_2' in file:
        # extract the data ID from the file name
        data_id = re.split(r'[^a-zA-Z0-9]+', file)[0]

        # count sequence for the file
        cmd=f"cat {raw_dir}/{file} |  grep \"^@\" | wc -l"
        read_count = int(subprocess.getoutput(cmd))

        # select the file for pair_end
        if read_count>1000:
            raw_pe_dict.setdefault(data_id, []).append(file)
            raw_pe_dict[data_id].sort()
        else:
            print(file,"have not-enough sequence to following work!")

######################################################################################
### sort out pair_end sequences
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
        print("### pair-end dataset constructed:",data_id,":",raw_pe_dict[data_id])
        
# Import the dictionary from the file using pickle
with open('raw_se_dict.pkl', 'wb') as f:
    pickle.dump(raw_se_dict, f)

with open('raw_pe_dict.pkl', 'wb') as f:
    pickle.dump(raw_pe_dict, f)