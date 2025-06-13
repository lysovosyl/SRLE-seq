#%%
from collections import defaultdict
from Bio import SeqIO
import argparse
import os
from itertools import product
import numpy as np
import csv
from scipy.stats import ttest_ind


parse = argparse.ArgumentParser()
parse.add_argument('-nuc_fq1_file')
parse.add_argument('-nuc_fq2_file')
parse.add_argument('-cyto_fq1_file')
parse.add_argument('-cyto_fq2_file')
parse.add_argument('-save_path')
parse.add_argument('-up_primer_z')
parse.add_argument('-down_primer_z')
args = parse.parse_args()

up_primer_z = args.up_primer_z
down_primer_z = args.down_primer_z
nuc_fq1_file  = args.nuc_fq1_file
nuc_fq2_file  = args.nuc_fq2_file
cyto_fq1_file  = args.cyto_fq1_file
cyto_fq2_file  = args.cyto_fq2_file
save_path = args.save_path

# up_primer_z = 'ATTATGAT'
# down_primer_z='GCTTAGTG'
#
# nuc_fq1_file  = '/mnt/dfc_data1/home/linyusen/tmp/FP180000985TL_L01_249_1.fq'
# nuc_fq2_file  = '/mnt/dfc_data1/home/linyusen/tmp/FP180000985TL_L01_249_2.fq'
# cyto_fq1_file  = '/mnt/dfc_data1/home/linyusen/tmp/FP180000948BL_L01_248_1.fq'
# cyto_fq2_file  = '/mnt/dfc_data1/home/linyusen/tmp/FP180000948BL_L01_248_2.fq'
# save_path = '/mnt/dfc_data1/home/linyusen/tmp/6mer_diversity_nuc_1.pkl'


from  tqdm import tqdm

def data_kmer_count(fq1,fq2):
    kmer_count = defaultdict(int)
    for record in tqdm(SeqIO.parse(fq1, 'fastq')):
        seq = str(record.seq)
        if up_primer_z in seq and down_primer_z in seq:
            start = seq.find(up_primer_z)+len(up_primer_z)
            end = seq.find(down_primer_z)
            kmer = seq[start:end]
            kmer_count[kmer]+=1

    for record in tqdm(SeqIO.parse(fq2, 'fastq')):
        seq = str(record.seq)
        if up_primer_z in seq and down_primer_z in seq:
            start = seq.find(up_primer_z)+len(up_primer_z)
            end = seq.find(down_primer_z)
            kmer = seq[start:end]
            kmer_count[kmer]+=1

    for kmer in list(kmer_count.keys()):
        if 'N' in kmer:
            kmer_count.pop(kmer)
        elif len(kmer) != 6:
            kmer_count.pop(kmer)
    return kmer_count


data = {}
data['nuc'] = data_kmer_count(nuc_fq1_file,nuc_fq2_file)
data['cyto'] = data_kmer_count(cyto_fq1_file,cyto_fq2_file)

f = open(os.path.join(save_path,'cyto.kmer.count.csv'),'w')
f.write('kmer,log2fc\n')
for kmer in data['cyto']:
    f.write(f"{kmer},{data['cyto'][kmer]}\n")
f.close()

f = open(os.path.join(save_path,'nuc.kmer.count.csv'),'w')
f.write('kmer,log2fc\n')
for kmer in data['nuc']:
    f.write(f"{kmer},{data['nuc'][kmer]}\n")
f.close()
#%%


for types in data:
    total_count  = 0
    for kmer in data[types]:
        total_count+=data[types][kmer]
    for kmer in data[types]:
        data[types][kmer] = data[types][kmer] / total_count * 1000000


def generate_6mers():
    bases = ['A', 'T', 'C', 'G']
    six_mers = [''.join(p) for p in product(bases, repeat=6)]
    return six_mers
kmer_list = generate_6mers()


kmer_info = {}
for kmer in kmer_list:
    temp = [kmer]
    cyto = []
    nuc = []

    for types in data:
        temp.append(data[types][kmer])
        if types == 'nuc':
            nuc.append(data[types][kmer])
        else:
            cyto.append(data[types][kmer])
    temp.append(np.log2(np.mean(cyto)/np.mean(nuc)))
    kmer_info[kmer] = [np.log2(np.mean(cyto)/np.mean(nuc))]

f = open(os.path.join(save_path,'cyto.kmer.csv'),'w')
f.write('kmer,log2fc\n')
for key in kmer_info:
    f.write(f'{key},{kmer_info[key][0]}\n')
f.close()

f = open(os.path.join(save_path,'nuc.kmer.csv'),'w')
f.write('kmer,log2fc\n')
for key in kmer_info:
    f.write(f'{key},{kmer_info[key][0]}\n')
f.close()


