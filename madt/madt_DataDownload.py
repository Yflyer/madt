#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Yufei,2023
# github.com/Yflyer
import os,subprocess,pickle,re
import argparse
from Bio import SeqIO
import shutil

######################################################################################
# madt_DataIndex.py -i SRA.list.txt -o raw_data
### Input parameters of script ###
def get_parser():
    parser = argparse.ArgumentParser(description="""Meta-Amplicon Database Toolkit\n
github.com/Yflyer/ \n
----------
Module: Data download\n
This scripts is to download rawdata and unzip them into fastq from sra \n
Script will continue to download until all SRA accession ID downloaded\n
""",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i','--Input',required=True,help='The txt file of SRA accession ID by line')
    parser.add_argument('-o','--SingleDir',default='single_data',help='Directory of outputing single data')

parser = get_parser()
args = parser.parse_args()
input_list = args.input
output_dir = args.output

### workflow:
# 1. prefetch download
# 2. check prefetch
# 3. (parallel) fastq-dump to extract fastq data 
print('----------------------------------------------------')

### 1. prefetch
# prefetch download:
print('sra files are downloading.')
subprocess.run(['prefetch', '--option-file', input_list, '--output-directory', output_dir])

### 2. check
# two kind of errors in SRA download:
# a. net conection failed
# b. missed by prefecch --file
# create completed and uncompleted list
completed = set()
for file in os.listdir(output_dir):
    if file.endswith('.sra*'):
        completed.add(file.split('.')[0])

with open(input_list) as f:
    all_ids = set(line.strip() for line in f)

uncompleted = all_ids - completed
with open('uncompleted.run.txt', 'w') as f:
    for item in uncompleted:
        f.write(item + '\n')

print(len(uncompleted), 'sra files are undownloaded. Now start to re-download')

# use while to repeatly check whether uncompleted sra number left. And then refresh it every time after a successful prefetch.
while len(uncompleted) > 0:
    subprocess.run(['prefetch', '--option-file', 'uncompleted.run.txt', '--output-directory', output_dir])
    # refresh the uncompleted list
    completed = set()
    for file in os.listdir('.'):
        if file.endswith('.sra'):
            completed.add(file.split('.')[0])
    uncompleted = all_ids - completed
    with open('uncompleted.run.txt', 'w') as f:
        for item in uncompleted:
            f.write(item + '\n')

print('All SRA list has been checked and downloaded')

### 3. fastq-dump
print('Now start fastq dump')
for prefix in completed:
    prefix = prefix.strip()
    subprocess.run(['fastq-dump', '-I', '--split-e', prefix + '*.sra*', '--outdir', output_dir])
#print('Used threads: ', os.sys.argv[2])
print('-----------------All work done!---------------------')
print('----------------------------------------------------')
