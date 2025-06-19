#%%
import os
import numpy as np
from collections import defaultdict
from Bio import SeqIO
import argparse
import os
from itertools import product
import csv
from scipy.stats import ttest_ind
from tqdm import tqdm
from multiprocessing import Pool
import re

parse = argparse.ArgumentParser(description="K-mer counting and plotting script for nuc/cyto comparison")
parse.add_argument('-fq1', type=str, required=True, help="R1 fastq file")
parse.add_argument('-fq2', type=str, required=True, help="R2 fastq file")
parse.add_argument('-s', type=str, required=True, help="Output directory for saving result figures and tables")
parse.add_argument('-up_flanking', type=str, required=True, help="Sequence of the upstream (5') anchor primer")
parse.add_argument('-down_flanking', type=str, required=True, help="Sequence of the downstream (3') anchor primer")
parse.add_argument('-kmer_length', type=int, required=True, help="Length of k-mers to count")
args = parse.parse_args()

fq1_file = args.fq1
fq2_file = args.fq2
save_path = args.s
up_flanking_region = args.up_flanking
down_flanking_region = args.down_flanking
kmer_length = args.kmer_length

# fq1_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/matched_fq1.fastq'
# fq2_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/matched_fq2.fastq'
# save_path = '/mnt/dfc_data1/home/linyusen/tmp/xinquan'
# up_flanking_region = 'ATTATGAT'
# down_flanking_region='GCTTAGTG'
# mode = 'kmer_diversity'
# kmer_length = 6

def data_kmer_count_fast(fq1, fq2, pattern):
    kmer_count = defaultdict(int)
    for file in [fq1, fq2]:
        for record in tqdm(SeqIO.parse(file, 'fastq')):
            seq = str(record.seq)
            match = pattern.search(seq)
            if match:
                kmer = match.group(1)
                kmer_count[kmer] += 1
    return kmer_count


pattern = re.compile(f"{up_flanking_region}([ACGT]+?){down_flanking_region}")
data = data_kmer_count_fast(fq1_file, fq2_file,pattern)
for i in list(data.keys()):
    if len(i) != kmer_length:
        data.pop(i)
filtered_data = {k: v for k, v in data.items() if len(k) == kmer_length}
sorted_data = dict(sorted(filtered_data.items(), key=lambda item: item[1], reverse=True))

total_count = sum(sorted_data.values())

f = open(os.path.join(save_path,f'kmer_diversity.validation.tsv'),'w')
w = csv.writer(f,delimiter='\t')
w.writerow(['kmer','count','CPM'])
for i in sorted_data:
    w.writerow([i,sorted_data[i],sorted_data[i]/total_count*1000000])
f.close()
print(f'{kmer_length}mer total_count:{total_count}')
