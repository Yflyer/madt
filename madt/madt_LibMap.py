#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Yufei,2023
# github.com/Yflyer
import pickle,os,time,subprocess
import pandas as pd
import shutil
from Bio import SeqIO

######################################################################################
### Input parameters of script ###
input_dir = 'pe_data'
dict_path = 'pe_dict.pkl'
blast_dir = 'blast'

Ref_path = '../Ref/Ecoil/Ecoil'
#Ref_path = '../Ref/mmseqs_80c/mmseqs_80c'

def mkdir(dir):
	if not os.path.exists(dir):
		os.makedirs(dir)
		print ("--- Create a new directory ---")
	else:
		print ("--- The directory has already existed! ---")
mkdir(blast_dir)

# Create an empty DataFrame to store the mean values
mean_df = pd.DataFrame(columns=['Run_ID','Lib_type', 'Mean_sstart', 'Mean_send','Match_ratio'])

# Import the dictionary from the file using pickle
with open(dict_path, 'rb') as f:
    seq_dict = pickle.load(f)

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
    cmd = f"blastn -query {blast_dir}/{data_id}_100.fasta -db {Ref_path} -outfmt '6 qseqid length qstart qend sseqid sstart send evalue bitscore'  -num_alignments 1 -out {blast_dir}/{blast_output} -num_threads 12"
    subprocess.run(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )

    ######################################################################################
    # Check the mapped status of each data
    if os.path.getsize(f"{blast_dir}/{blast_output}") != 0:
        df = pd.read_csv(f"{blast_dir}/{blast_output}",sep='\t', header=None) 
        Match_ratio = round(df.shape[0]/100,3)
        Mean_length,Mean_qstart,Mean_qend,Mean_sstart,Mean_send=round(df.iloc[:,[1,2,3,5,6]].mean(),0)
                                                
        # check requirement of map length  
        if Mean_length > 150:
            print(f"{data_id} (query position:{Mean_qstart}-{Mean_qend}) matched at the ref position from {Mean_sstart} to {Mean_send} (length: {Mean_length})")
            data_stat = pd.DataFrame({'Run_ID': data_id, 'Lib_type': input_dir,'Mean_sstart': Mean_sstart, 'Mean_send': Mean_send,'Match_ratio':Match_ratio}, index=[data_id])
            mapped_num += 1
        else:
            print(f"{data_id} can not match to reference sequences")
            Mean_qstart = Mean_qend = 0
            data_stat = pd.DataFrame({'Run_ID': data_id, 'Lib_type': input_dir,'Mean_sstart': 0, 'Mean_send': 0,'Match_ratio':0}, index=[data_id])
            unmapped_num += 1

    else:
        print(f"{data_id} can not match to reference sequences")
        Mean_qstart = Mean_qend = 0
        data_stat = pd.DataFrame({'Run_ID': data_id, 'Lib_type': input_dir,'Mean_sstart': 0, 'Mean_send': 0,'Match_ratio':0}, index=[data_id])
        unmapped_num += 1

    mean_df = pd.concat([mean_df, data_stat], ignore_index=True)

    ######################################################################################
    # check reverse between R1 and R2
    if  Mean_qstart>Mean_qend:
        print(f'* R1 and R2 files are exchange in {data_id}. Now reverse merged files to forward orientation.')
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
        
        #shutil.move(f"{input_dir}/Rev_{filename}", f"{input_dir}/{filename}")

    print('------------------------------------------------------')

mean_df.to_csv(f'{input_dir}_Ecoilmap_summary.csv', index=False)


print('----------------------------------------------------')
print('----------------------Summary-----------------------')
print(f'In {len(seq_dict)} data of {input_dir}:')
print(f'{mapped_num} data are mapped to reference sequences.')
print(f'{rev_num} data are mapped but reverse.')
print(f'{unmapped_num} data can not be mapped.')
print('----------------------------------------------------')
            


