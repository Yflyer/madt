#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Yufei,2023
# github.com/Yflyer
import os,time,subprocess,re
import pandas as pd
import shutil
from Bio import SeqIO
import argparse
######################################################################################
# madt_LibMap.py -i pair_data -b blast -r ../Ref/Ecoil/Ecoil
### Input parameters of script ###
def get_parser():
    parser = argparse.ArgumentParser(description="""Meta-Amplicon Database Toolkit\n
github.com/Yflyer/ \n
----------
Module: Check library\n
This scripts is to index rawdata downloaded from sra\n
Sort out single-end and pair-end sequences; remove files without enough data\n
""",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i','--input_dir',required=True,help='Directory of input data')
    parser.add_argument('-r','--Ref_path',required=True,help='Path of blast db for reference')
    parser.add_argument('-b','--Blast_dir',default='blast',help='Directory of blast output data')
    parser.add_argument('-ke','--Keep_empty',default=False,help='Keep the record of sequence file which can not be mapped')
    parser.add_argument('-cl','--Coverage_length',default=200,help='Requirement of overlap length for both query sequences and reference sequences')
    return parser

###################
parser = get_parser()
args = parser.parse_args()
input_dir = args.input_dir
blast_dir = args.Blast_dir
Ref_path = args.Ref_path
Keep_empty = args.Keep_empty
Coverage_len  = int(args.Coverage_length)
#input_dir = 'pe_data'
#blast_dir = 'blast'
#Ref_path = '../Ref/Ecoil/Ecoil'
#Ref_path = '../Ref/mmseqs_80c/mmseqs_80c'

def mkdir(dir):
	if not os.path.exists(dir):
		os.makedirs(dir)
		print ("--- Create a new directory ---")
	else:
		print ("--- The directory has already existed! ---")
def bac_region_estimate(sstart,send):
    V3_f = range(325, 360)
    V4_f = range(500, 545)
    V4_r = range(775, 820)
    V5_r = range(890, 935)
    
    # detect forward primer region
    if sstart in V3_f:
        Forward_primer_region = 'V3f'
    elif sstart in V4_f:
        Forward_primer_region = 'V4f'
    elif sstart > V4_f.stop:
        Forward_primer_region = 'V4nf'
    elif sstart > V3_f.stop:
        Forward_primer_region = 'V3nf'
    else:
        Forward_primer_region = 'Unf'
    # detect reverse primer region
    if send in V4_r:
        Reverse_primer_region = 'V4r'
    elif send in V5_r:
        Reverse_primer_region = 'V5r'
    elif send < V4_r.start:
        Reverse_primer_region = 'V4nr'
    elif send > V5_r.start:
        Reverse_primer_region = 'V5nr'
    else:
        Reverse_primer_region = 'Unf'

    Region = Forward_primer_region +'-'+Reverse_primer_region
    return(Region)


mkdir(blast_dir)

# Create an empty DataFrame to store the mean values
mean_df = pd.DataFrame(columns=['Run_ID','Lib_type', 'Mean_sstart', 'Mean_send','Mean_length','Mean_qlen','Match_ratio','Region','Fwd_Pri_len','Rev_Pri_len'
])

# Import the dictionary from the file using pickle

seq_dict = {}
### index single data
for file in os.listdir(input_dir):
    if file.endswith('.fasta'):
        # extract the data ID from the file name
        data_id = re.split(r'[^a-zA-Z0-9]+', file)[0]
        seq_dict[data_id] = file

mapped_num = rev_num = unmapped_num = 0
for data_id in seq_dict:
    
    ######################################################################################
    # Run BLAST and count the number of bacterial hits
    filename = seq_dict[data_id]
    blast_output = data_id + "_blast.tsv"

    print('------------------------------------------------------')
    print('------{}: now start blast file {}'.format(time.strftime("%Y-%m-%d %X"),filename))
    
    cmd = f"seqkit sample -n 100 {input_dir}/{filename} > {blast_dir}/{data_id}_100.fasta"
    subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
    cmd = f"blastn -query {blast_dir}/{data_id}_100.fasta -db {Ref_path} -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length evalue bitscore'  -num_alignments 1 -out {blast_dir}/{blast_output} -num_threads 12"
    subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )

    ######################################################################################
    # Check the mapped status of each data
    if os.path.getsize(f"{blast_dir}/{blast_output}") != 0:
        df = pd.read_csv(f"{blast_dir}/{blast_output}",sep='\t', header=None) 
        Match_ratio = round(df.shape[0]/100,3)
        Mean_qlen,Mean_qstart,Mean_qend,Mean_slen,Mean_sstart,Mean_send,Mean_length=round(df.iloc[:,[1,2,3,5,6,7,8]].mean(),0)

        ###########################################
        # check requirement of map length  
        if abs(Mean_send-Mean_sstart) >= Coverage_len and abs(Mean_qend-Mean_qstart) >= Coverage_len:
            # check reverse between R1 and R2
            if  Mean_sstart>Mean_send:
                print(f'* R1 and R2 files are exchange in {data_id} according to the ref position from {Mean_sstart} to {Mean_send}. Now reverse the merged file.')
                rev_num += 1
                
                Input_fasta = SeqIO.parse(f"{input_dir}/{filename}", "fasta")

                # Create a list of records to write to the output file
                new_records = []

                # Reverse complement the sequence and update the record
                for record in Input_fasta:
                    record.seq = record.seq[::-1]
                    new_records.append(record)

                with open(f"{input_dir}/Rev_{filename}", "w") as Rev_fasta: 
                    SeqIO.write(new_records, Rev_fasta, "fasta")
                # replace reverse file
                shutil.move(f"{input_dir}/Rev_{filename}", f"{input_dir}/{filename}")
                data_stat = pd.DataFrame({'Run_ID': data_id, 
                                          'Lib_type': input_dir,
                                          'Mean_sstart': Mean_send, 
                                          'Mean_send': Mean_sstart,
                                          'Mean_length':Mean_length,
                                          'Mean_qlen':Mean_qlen,
                                          'Match_ratio':Match_ratio,
                                          'Region':Region,
                                          'Fwd_Pri_len':[Mean_qstart],
                                          'Rev_Pri_len':[Mean_qlen-Mean_qend]}, index=[data_id])
                mapped_num += 1
                # concat new row to output
                Region = bac_region_estimate(Mean_send,Mean_sstart)
                mean_df = pd.concat([mean_df, data_stat], ignore_index=True)
            else:
                mapped_num += 1
                # print basic message
                print(f"{data_id} (query position:{Mean_qstart}-{Mean_qend}) matched at the ref position from {Mean_sstart} to {Mean_send} (length: {Mean_length})")
                Region = bac_region_estimate(Mean_sstart,Mean_send)
                data_stat = pd.DataFrame({'Run_ID': data_id, 
                                        'Lib_type': input_dir,
                                        'Mean_sstart': Mean_sstart, 
                                        'Mean_send': Mean_send,
                                        'Mean_length':Mean_length,
                                        'Mean_qlen':Mean_qlen,
                                        'Match_ratio':Match_ratio,
                                        'Region':Region,
                                        'Fwd_Pri_len':[Mean_qstart],
                                        'Rev_Pri_len':[Mean_qlen-Mean_qend]}, index=[data_id])
                
                # concat new row to output
                mean_df = pd.concat([mean_df, data_stat], ignore_index=True)
                 
            

        else:
            print(f"{data_id} mapped by ref positions from {Mean_sstart} to {Mean_send} can not meet the requirement of overlap coverage")
            unmapped_num += 1 
            if Keep_empty :
                Mean_qstart = Mean_qend = 0
                data_stat = pd.DataFrame({'Run_ID': data_id, 
                                          'Lib_type': input_dir,
                                          'Mean_sstart': Mean_sstart, 
                                          'Mean_send': Mean_send,
                                          'Mean_length':Mean_length,
                                          'Mean_qlen':Mean_qlen,
                                          'Match_ratio':Match_ratio,                                        
                                          'Region':'Unmapped',
                                          'Fwd_Pri_len':0,
                                          'Rev_Pri_len':0}, index=[data_id])
                
                # concat new row to output
                mean_df = pd.concat([mean_df, data_stat], ignore_index=True)    

    else:
        print(f"{data_id} blast result is empty!")
        unmapped_num += 1 
        if Keep_empty :
                Mean_qstart = Mean_qend = 0
                data_stat = pd.DataFrame({'Run_ID': data_id, 
                                        'Lib_type': input_dir,
                                        'Mean_sstart': 0, 
                                        'Mean_send': 0,
                                        'Mean_length':0,
                                        'Mean_qlen':0,
                                        'Match_ratio':0,
                                        'Region':'Unmapped',
                                        'Fwd_Pri_len':0,
                                        'Rev_Pri_len':0}, index=[data_id])
                
                # concat new row to output
                mean_df = pd.concat([mean_df, data_stat], ignore_index=True)   

    ######################################################################################
            
    
print('------------------------------------------------------')
mean_df.to_csv(f'{input_dir}_Ecoilmap_summary.csv', index=False)
print('----------------------------------------------------')
print('----------------------Summary-----------------------')
print(f'In {len(seq_dict)} data of {input_dir}:')
print(f'{mapped_num} data are mapped to reference sequences.')
print(f'{rev_num} data are mapped but reverse.')
print(f'{unmapped_num} data can not be mapped.')
print('----------------------------------------------------')
            


