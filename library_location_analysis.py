#%%
from collections import defaultdict
from Bio import SeqIO
import argparse
import os
from itertools import product
import numpy as np
import csv
from scipy.stats import ttest_ind
from tqdm import tqdm

parse = argparse.ArgumentParser(description="K-mer counting and plotting script for nuc/cyto comparison")
parse.add_argument('-nuc_fq1_file', type=str, required=True, help="nuclear sample R1 fastq file")
parse.add_argument('-nuc_fq2_file', type=str, required=True, help="nuclear sample R2 fastq file")
parse.add_argument('-cyto_fq1_file', type=str, required=True, help="cytoplasmic sample R1 fastq file")
parse.add_argument('-cyto_fq2_file', type=str, required=True, help="cytoplasmic sample R2 fastq file")
parse.add_argument('-s', type=str, required=True, help="Path to save result figures or tables")
parse.add_argument('-up_flanking', type=str, required=True, help="Upstream primer sequence (z)")
parse.add_argument('-down_flanking', type=str, required=True, help="Downstream primer sequence (z)")
parse.add_argument('-mode', choices=['kmer_complexity', 'fragment_complexity'], required=True,help="Processing mode: library_complexity or cell_complexity")
parse.add_argument('-kmer_length', type=int, required=True, help="Length of k-mers to count")
args = parse.parse_args()

up_primer_z = args.up_flanking
down_primer_z = args.down_flanking
nuc_fq1_file  = args.nuc_fq1_file
nuc_fq2_file  = args.nuc_fq2_file
cyto_fq1_file  = args.cyto_fq1_file
cyto_fq2_file  = args.cyto_fq2_file
save_path = args.s
mode = args.mode
kmer_length = args.kmer_length

# up_primer_z = 'ATTATGAT'
# down_primer_z='GCTTAGTG'
# nuc_fq1_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/kmer/6mer_Nuc_Sample3/FP180000985TL_L01_249_1.fq'
# nuc_fq2_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/kmer/6mer_Nuc_Sample3/FP180000985TL_L01_249_2.fq'
# cyto_fq1_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/kmer/6mer_Cyto_Sample3/FP180000948BL_L01_248_1.fq'
# cyto_fq2_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/kmer/6mer_Cyto_Sample3/FP180000948BL_L01_248_2.fq'
# save_path = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/'
# mode = 'library_complexity'
# kmer_length = 6

if mode == 'kmer_complexity':
    if not isinstance(kmer_length, int):
        raise ValueError("when in library_complexity mode, the kmer_length should be a number")
#%%
def data_kmer_count(fq1,fq2,kmer_length,mode):
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

    if mode == 'library_complexity':
        for kmer in list(kmer_count.keys()):
            if 'N' in kmer:
                kmer_count.pop(kmer)
            elif len(kmer) != kmer_length:
                kmer_count.pop(kmer)
    else:
        pass
    return kmer_count


data = {}
data['nuc'] = data_kmer_count(nuc_fq1_file,nuc_fq2_file,kmer_length,mode)
data['cyto'] = data_kmer_count(cyto_fq1_file,cyto_fq2_file,kmer_length,mode)
# data['total'] = data_kmer_count(total_fq1_file,total_fq2_file)
#%%
total_count = 0
for i in data['nuc']:
    total_count += data['nuc'][i]
for i in data['cyto']:
    total_count += data['cyto'][i]

#%%
f = open(os.path.join(save_path,'cyto.kmer.count.csv'),'w')
f.write('kmer,count,CPM\n')
for kmer in data['cyto']:
    cpm = data['cyto'][kmer] / total_count * 1000000
    f.write(f"{kmer},{data['cyto'][kmer]},{cpm}\n")
f.close()

f = open(os.path.join(save_path,'nuc.kmer.count.csv'),'w')
f.write('kmer,count,CPM\n')
for kmer in data['nuc']:
    cpm = data['nuc'][kmer] / total_count * 1000000
    f.write(f"{kmer},{data['nuc'][kmer]},{cpm}\n")
f.close()
#%%

import matplotlib.pyplot as plt
import numpy as np

nuc_dict = data['nuc']
cyto_dict = data['cyto']

# 确保 kmer 键一致
common_kmers = set(nuc_dict) & set(cyto_dict)

x = []
y = []
colors = []

for kmer in common_kmers:
    cyto_val = cyto_dict[kmer]
    nuc_val = nuc_dict[kmer]

    # 避免除0或log0，加一个极小值
    cyto_val = cyto_val if cyto_val > 0 else 1e-6
    nuc_val = nuc_val if nuc_val > 0 else 1e-6

    log2fc = np.log2(nuc_val / cyto_val)

    if log2fc < -0.58:
        colors.append("blue")
    elif log2fc > 0.56:
        colors.append("red")
    else:
        colors.append("gray")

    x.append(cyto_dict[kmer])
    y.append(nuc_dict[kmer])

# 绘图（log scale 更直观）
plt.figure(figsize=(6, 6))
plt.scatter(x, y, c=colors, alpha=0.6, s=10)
plt.plot([min(x + y), max(x + y)], [min(x + y), max(x + y)], color='black', linestyle='--')  # 对角线

plt.xlabel("cyto")
plt.ylabel("nuc")
plt.title("K-mer abundance scatter plot with log2FC coloring")
plt.tight_layout()
plt.savefig(os.path.join(save_path,'nuc.cyto.kmer.count.png'))
#%%
for types in data:
    total_count  = 0
    for kmer in data[types]:
        total_count+=data[types][kmer]
    for kmer in data[types]:
        data[types][kmer] = data[types][kmer] / total_count * 1000000


def generate_kmers(kmer_length):
    bases = ['A', 'T', 'C', 'G']
    six_mers = [''.join(p) for p in product(bases, repeat=kmer_length)]
    return six_mers
kmer_list = generate_kmers(kmer_length)


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
    nes = np.log2((np.mean(nuc) + 0.5) / (np.mean(cyto) + 0.5))
    kmer_info[kmer] = [np.log2(np.mean(cyto)/np.mean(nuc)),nes]

f = open(os.path.join(save_path,'kmer.info.csv'),'w')
f.write('kmer,log2fc,nes\n')
for key in kmer_info:
    f.write(f'{key},{kmer_info[key][0]},{kmer_info[key][1]}\n')
f.close()

