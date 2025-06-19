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
import matplotlib.pyplot as plt
import pysam
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import islice
import yaml

#%%
parse = argparse.ArgumentParser(description="K-mer counting and plotting script for nuc/cyto comparison")
parse.add_argument('-fq1', type=str, required=True, help="R1 fastq file")
parse.add_argument('-fq2', type=str, required=True, help="R2 fastq file")
parse.add_argument('-s', type=str, required=True, help="Path to save result tables")
parse.add_argument('-up_flanking', type=str, required=True, help="Upstream primer sequence (z)")
parse.add_argument('-down_flanking', type=str, required=True, help="Downstream primer sequence (z)")
parse.add_argument('-gene_chr', type=str, required=True, help="Downstream primer sequence (z)")
parse.add_argument('-gene_region_start', type=int, required=True, help="Downstream primer sequence (z)")
parse.add_argument('-gene_region_end', type=int, required=True, help="Downstream primer sequence (z)")
parse.add_argument('-bowtie', type=str, required=True, help="Downstream primer sequence (z)")
parse.add_argument('-genome_index', type=str, required=True, help="Downstream primer sequence (z)")
parse.add_argument('-samtools', type=str, required=True, help="Downstream primer sequence (z)")
parse.add_argument('-thread', type=int, required=True, help="Downstream primer sequence (z)")
parse.add_argument('-config',required=True,help='')

args = parse.parse_args()

fq1_file = args.fq1
fq2_file = args.fq2

save_path = args.s
up_flanking_region = args.up_flanking
down_flanking_region = args.down_flanking
gene_chr = args.gene_chr
gene_region_start = args.gene_region_start
gene_region_end = args.gene_region_end
thread = args.thread
config_file = args.config
with open(config_file, 'r') as file:
    config_data = yaml.safe_load(file)

bowtie = config_data['software']['bowtie']
genome_index = config_data['database']['genome_index']
samtools = config_data['software']['samtools']
pear = config_data['software']['pear']
#%%
# fq1_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/MALAT1.library.R1.fastq'
# fq2_file = '/mnt/dfc_data1/home/linyusen/tmp/xinquan/MALAT1.library.R2.fastq'
# save_path = '/mnt/dfc_data1/home/linyusen/tmp/xinquan'
# up_flanking_region = 'ATTATGAT'
# down_flanking_region='GCTTAGTG'
# gene_chr = 'chr11'
# gene_region_start = 65497688
# gene_region_end = 65506516
# thread = 64
# bowtie = '/mnt/dfc_data1/home/linyusen/miniconda3/envs/lys/bin/bowtie'
# genome_index = '/mnt/dfc_data2/project/linyusen/database/02_hg/hg38/bowtie_index/hg38'
# samtools = '/mnt/dfc_data1/software/anaconda2/envs/samtools-1.9/bin/samtools'
# pear = '/mnt/dfc_data1/home/linyusen/miniconda3/envs/lys/bin/pear'
#%%
def read_fastq_pairs(fq1_path, fq2_path, batch_size=10000):
    with open(fq1_path) as fq1, open(fq2_path) as fq2:
        fq1_iter = FastqGeneralIterator(fq1)
        fq2_iter = FastqGeneralIterator(fq2)

        while True:
            batch1 = list(islice(fq1_iter, batch_size))
            batch2 = list(islice(fq2_iter, batch_size))

            if not batch1 or not batch2:
                break

            # 补全少的一侧（可选，根据你对残缺配对的容忍度）
            min_len = min(len(batch1), len(batch2))
            for i in range(min_len):
                title1, seq1, qual1 = batch1[i]
                title2, seq2, qual2 = batch2[i]
                yield (title1, seq1, qual1, title2, seq2, qual2)
def read_fastq_single(fq_path, batch_size=10000):
    with open(fq_path) as fq:
        fq_iter = FastqGeneralIterator(fq)

        while True:
            batch1 = list(islice(fq_iter, batch_size))

            if not batch1:
                break

            # 补全少的一侧（可选，根据你对残缺配对的容忍度）
            min_len = len(batch1)
            for i in range(min_len):
                title1, seq1, qual1 = batch1[i]
                yield (title1, seq1, qual1)
def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]
up_flanking_region_re = reverse_complement(up_flanking_region)
down_flanking_region_re = reverse_complement(down_flanking_region)

#%%
pear_out = os.path.join(save_path,'pear')
cmd = f'{pear} -f {fq1_file} -r {fq2_file} -j {thread} -o {pear_out}'
os.system(cmd)
single_fastq = os.path.join(save_path,'pear.assembled.fastq')

pattern = re.compile(f"{up_flanking_region}([ACGT]+?){down_flanking_region}")
pattern_re = re.compile(f"{down_flanking_region_re}([ACGT]+?){up_flanking_region_re}")
save_fq = os.path.join(save_path,"output.fastq")

with open(save_fq, "w") as fq_out:
    for i, (title1, seq1, qual1) in enumerate(read_fastq_single(single_fastq)):
        seq1 = str(seq1)
        match1 = pattern.search(seq1)
        match1_re = pattern_re.search(seq1)
        fq_seq = False

        if match1:
            start = match1.start(1)
            end = match1.end(1)
            fq_seq = match1.group(1)
            fq_qual = qual1[start:end]


        if match1_re:
            start = match1_re.start(1)
            end = match1_re.end(1)
            fq_seq = match1_re.group(1)
            fq_qual = qual1[start:end]

        if fq_seq:
            if len(fq_seq) > 10:
                fq_out.write(f"@{title1}\n{fq_seq}\n+\n{fq_qual}\n")



bam_file = os.path.join(save_path,'dna_fragment.alignment.bam')
cmd = f'{bowtie} -v 2 -a -k 10 -p {thread} -S {genome_index} {save_fq} | {samtools} view -bS > {bam_file}'
os.system(cmd)
sorted_bam_file = os.path.join(save_path,'dna_fragment.alignment.sorted.bam')
cmd = f'{samtools} sort -@ {thread} -o {sorted_bam_file} {bam_file} '
os.system(cmd)
cmd = f'{samtools} index {sorted_bam_file}'
os.system(cmd)

bam = pysam.AlignmentFile(sorted_bam_file, "rb")
coverage = [0] * (gene_region_end - gene_region_start)
coverage_x = []
for read in bam.fetch(gene_chr, gene_region_start, gene_region_end):
    print(read)
    if read.is_unmapped:
        continue
    for pos in read.get_reference_positions():
        if gene_region_start <= pos < gene_region_end:
            coverage[pos - gene_region_start] += 1

output_bed = os.path.join(save_path,'dna_fragment.coverage.bed')
with open(output_bed, "w") as out:
    for i, cov in enumerate(coverage):
        chrom_start = gene_region_start + i
        chrom_end = chrom_start + 1
        coverage_x.append(gene_region_start + i)
        out.write(f"{gene_chr}\t{chrom_start}\t{chrom_end}\t{cov}\n")
#%%
bwith = 3
fig, (ax1) = plt.subplots(1, 1, figsize=(16, 3),dpi=300)
labelsize = 14
ax1.bar(coverage_x,coverage,width=1)
ax1.set_title('{}:{}-{}'.format(gene_chr,gene_region_start,gene_region_end))
ax1.set_xlim([gene_region_start,gene_region_end])
ax1.spines['left'].set_linewidth(bwith)
ax1.spines['top'].set_linewidth(bwith)
ax1.spines['right'].set_linewidth(bwith)
ax1.spines['bottom'].set_linewidth(bwith)
ax1.tick_params(axis='both', which='major', width=bwith, length=bwith*3,labelsize=labelsize)
figure_save_path = os.path.join(save_path,'dna_fragment.coverage.png')
plt.tight_layout()
plt.savefig(figure_save_path)
plt.show()
