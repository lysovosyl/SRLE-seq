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
import matplotlib.pyplot as plt
import numpy as np
import csv

parse = argparse.ArgumentParser(description="K-mer counting and plotting script for nuc/cyto comparison")
parse.add_argument('-nuc_fq1', type=str, required=True, help="nuclear sample R1 fastq file")
parse.add_argument('-nuc_fq2', type=str, required=True, help="nuclear sample R2 fastq file")
parse.add_argument('-cyto_fq1', type=str, required=True, help="cytoplasmic sample R1 fastq file")
parse.add_argument('-cyto_fq2', type=str, required=True, help="cytoplasmic sample R2 fastq file")
parse.add_argument('-s', type=str, required=True, help="Path to save result figures or tables")
parse.add_argument('-up_flanking', type=str, required=True, help="Upstream primer sequence (z)")
parse.add_argument('-down_flanking', type=str, required=True, help="Downstream primer sequence (z)")
parse.add_argument('-mode', choices=['kmer_complexity', 'fragment_complexity'], default='kmer_complexity',help="Processing mode: library_complexity or cell_complexity")
parse.add_argument('-kmer_length', type=int, default=None, help="Length of k-mers to count")
args = parse.parse_args()

up_primer_z = args.up_flanking
down_primer_z = args.down_flanking
nuc_fq1_file  = args.nuc_fq1
nuc_fq2_file  = args.nuc_fq2
cyto_fq1_file  = args.cyto_fq1
cyto_fq2_file  = args.cyto_fq2
save_path = args.s
mode = args.mode
kmer_length = args.kmer_length

# up_primer_z = 'ATTATGAT'
# down_primer_z='GCTTAGTG'
# nuc_fq1_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/kmer/6mer_Nuc_Sample3/test_1.fq'
# nuc_fq2_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/kmer/6mer_Nuc_Sample3/test_2.fq'
# cyto_fq1_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/kmer/6mer_Cyto_Sample3/test_1.fq'
# cyto_fq2_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/kmer/6mer_Cyto_Sample3/test_2.fq'
# save_path = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/'
# mode = 'kmer_complexity'
# kmer_length = 6

if mode == 'kmer_complexity':
    if not isinstance(kmer_length, int):
        raise ValueError("when in library_complexity mode, the kmer_length should be a number")

#%%
# Function to count k-mers extracted between two primer sequences from paired-end FASTQ files
def data_kmer_count(fq1, fq2, kmer_length, mode):
    kmer_count = defaultdict(int)  # Dictionary to count occurrences of each k-mer

    # Process forward reads (fq1)
    for record in tqdm(SeqIO.parse(fq1, 'fastq'), desc="Processing FASTQ", unit="reads"):
        seq = str(record.seq)
        # Check if both primers are found in the read
        if up_primer_z in seq and down_primer_z in seq:
            start = seq.find(up_primer_z) + len(up_primer_z)  # Start position after the upstream primer
            end = seq.find(down_primer_z)                     # End position before the downstream primer
            kmer = seq[start:end]                             # Extract the k-mer between primers
            kmer_count[kmer] += 1                             # Count the k-mer

    # Process reverse reads (fq2) in the same way
    for record in tqdm(SeqIO.parse(fq2, 'fastq'), desc="Processing FASTQ", unit="reads"):
        seq = str(record.seq)
        if up_primer_z in seq and down_primer_z in seq:
            start = seq.find(up_primer_z) + len(up_primer_z)
            end = seq.find(down_primer_z)
            kmer = seq[start:end]
            kmer_count[kmer] += 1

    # If the mode is set to 'library_complexity', filter out invalid kmers
    if mode == 'library_complexity':
        for kmer in list(kmer_count.keys()):
            if 'N' in kmer:                      # Remove kmers containing ambiguous base 'N'
                kmer_count.pop(kmer)
            elif len(kmer) != kmer_length:       # Remove kmers that are not of the expected length
                kmer_count.pop(kmer)
    else:
        pass  # No filtering for other modes

    return kmer_count  # Return the dictionary of k-mer counts

#%%
def generate_kmers(kmer_length):
    bases = ['A', 'T', 'C', 'G']
    six_mers = [''.join(p) for p in product(bases, repeat=kmer_length)]
    return six_mers
kmer_complexity_out_dict = {}
kmer_list = []
if mode == 'kmer_complexity':
    if isinstance(kmer_length, int):
        kmer_list = generate_kmers(kmer_length)
        for kmer in kmer_list:
            kmer_complexity_out_dict[kmer] = {'cyto':{},'nuc':{}}

#%%
data = {}
data['nuc'] = data_kmer_count(nuc_fq1_file, nuc_fq2_file, kmer_length, mode)   # k-mer counting for nuclear sample
data['cyto'] = data_kmer_count(cyto_fq1_file, cyto_fq2_file, kmer_length, mode) # k-mer counting for cytoplasmic sample

#%%

total_count = 0
for i in data['nuc']:
    total_count += data['nuc'][i]
for i in data['cyto']:
    total_count += data['cyto'][i]

for kmer in kmer_complexity_out_dict:
    if kmer in data['cyto']:
        cpm = data['cyto'][kmer] / total_count * 1000000
        kmer_complexity_out_dict[kmer]['cyto']['cpm'] = cpm
        kmer_complexity_out_dict[kmer]['cyto']['count'] = data['cyto'][kmer]
    else:
        kmer_complexity_out_dict[kmer]['cyto']['cpm'] = 0
        kmer_complexity_out_dict[kmer]['cyto']['count'] = 0


for kmer in kmer_complexity_out_dict:
    if kmer in data['nuc']:
        cpm = data['nuc'][kmer] / total_count * 1000000
        kmer_complexity_out_dict[kmer]['nuc']['cpm'] = cpm
        kmer_complexity_out_dict[kmer]['nuc']['count'] = data['nuc'][kmer]
    else:
        kmer_complexity_out_dict[kmer]['nuc']['cpm'] = 0
        kmer_complexity_out_dict[kmer]['nuc']['count'] = 0


#%%
if mode == 'kmer_complexity':
    for kmer in kmer_complexity_out_dict:
        count_cyto = kmer_complexity_out_dict[kmer]['cyto']['count']
        cpm_cyto = kmer_complexity_out_dict[kmer]['cyto']['cpm']
        count_nuc = kmer_complexity_out_dict[kmer]['nuc']['count']
        cpm_nuc = kmer_complexity_out_dict[kmer]['cyto']['cpm']
        nes = np.log2((cpm_nuc + 0.5) / (cpm_cyto + 0.5))
        log2fc = np.log2((count_nuc+0.000001) / (count_cyto+0.000001))
        kmer_complexity_out_dict[kmer]['nes'] = nes
        kmer_complexity_out_dict[kmer]['log2fc'] = log2fc

    f = open(os.path.join(save_path,'kmer_complexity.info.csv'),'w')
    w = csv.writer(f)
    w.writerow(['kmer','cyto_count','nuc_count','log2fc','nes'])
    for kmer in kmer_complexity_out_dict:
        w.writerow([kmer,kmer_complexity_out_dict[kmer]['cyto']['count'],kmer_complexity_out_dict[kmer]['nuc']['count'],
                    kmer_complexity_out_dict[kmer]['log2fc'],kmer_complexity_out_dict[kmer]['nes']])
    f.close()

    nuc_dict = data['nuc']
    cyto_dict = data['cyto']
    common_kmers = set(nuc_dict) & set(cyto_dict)
    x = []
    y = []
    colors = []
    for kmer in common_kmers:
        cyto_val = cyto_dict[kmer]
        nuc_val = nuc_dict[kmer]
        cyto_val = cyto_val if cyto_val > 0 else 1e-6
        nuc_val = nuc_val if nuc_val > 0 else 1e-6
        log2fc = np.log2(nuc_val / cyto_val)

        if log2fc < -0.58:
            colors.append("blue")
        elif log2fc > 0.56:
            colors.append("red")
        else:
            colors.append("gray")

        if log2fc > 10 or log2fc < -10:
            continue
        x.append(cyto_dict[kmer])
        y.append(nuc_dict[kmer])
    plt.figure(figsize=(6, 6))
    plt.scatter(x, y, c=colors, alpha=0.6, s=10)
    plt.plot([min(x + y), max(x + y)], [min(x + y), max(x + y)], color='black', linestyle='--')
    plt.xlabel("cyto")
    plt.ylabel("nuc")
    plt.title("K-mer abundance scatter plot with log2FC coloring")
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, 'kmer_complexity.kmer.count.png'))
