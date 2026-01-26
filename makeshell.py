

data = {'0h-1':['TTAAGTAGAGGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'24h-1':['TGATGCACATCTGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'48h-1':['TTCGATAGCAATTCGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'72h-1':['TGATCGATCCAGTTAGGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'96h-1':['TACGATCGATACACGATCGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'24h-2':['TATCATGCTTAGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'48h-2':['TCGATTGCTCGACGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'72h-2':['TATCGATAGTTGCTTGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'96h-2':['TCGATCGATTTGAGCCTGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'24h-3':['TTAAGTAGAGGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'48h-3':['TGATGCACATCTGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'72h-3':['TTCGATAGCAATTCGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT'],
'96h-3':['TGATCGATCCAGTTAGGCTTTATATATCTTGTGGAAAGGACGA','TTAAAGCAGCGTATCCACATAGCGT']}


import os
pathon = '/mnt/dfc_data1/home/linyusen/miniconda/envs/lys/bin/python'
for key in data:
    if '-1' in key:
        fq1 = '/mnt/dfc_data3/project/linyusen/00.CleanData/BJ-MORF-A/BJ-MORF-A_1.clean.fq'
        fq2 = '/mnt/dfc_data3/project/linyusen/00.CleanData/BJ-MORF-A/BJ-MORF-A_2.clean.fq'
        save_path = '/mnt/dfc_data3/project/linyusen/00.CleanData/BJ-MORF-A'
        up = data[key][0][-6:]
        cmd = f'{pathon} /mnt/dfc_data2/project/linyusen/project/72_xinquan_rnaloc/github/kmer_diversity_evaluation.py \
        -fq1 {fq1} \
        -fq2 {fq2} \
        -s {os.path.join(save_path,key)} \
        -up_flanking {up} \
        -down_flanking ACGCGT \
        -primer1 {data[key][0]} \
        -primer2 {data[key][1]}'
    elif '-2' in key:
        fq1 = '/mnt/dfc_data3/project/linyusen/00.CleanData/BJ-MORF-B/BJ-MORF-B_1.clean.fq'
        fq2 = '/mnt/dfc_data3/project/linyusen/00.CleanData/BJ-MORF-B/BJ-MORF-B_2.clean.fq'
        save_path = '/mnt/dfc_data3/project/linyusen/00.CleanData/BJ-MORF-B'
        up = data[key][0][-6:]
        cmd = f'{pathon} /mnt/dfc_data2/project/linyusen/project/72_xinquan_rnaloc/github/kmer_diversity_evaluation.py \
                -fq1 {fq1} \
                -fq2 {fq2} \
                -s {os.path.join(save_path, key)} \
                -up_flanking {up} \
                -down_flanking ACGCGT \
                -primer1 {data[key][0]} \
                -primer2 {data[key][1]}'
    elif '-3' in key:
        fq1 = '/mnt/dfc_data3/project/linyusen/00.CleanData/BJ-MORF-C/BJ-MORF-C_1.clean.fq'
        fq2 = '/mnt/dfc_data3/project/linyusen/00.CleanData/BJ-MORF-C/BJ-MORF-C_2.clean.fq'
        save_path = '/mnt/dfc_data3/project/linyusen/00.CleanData/BJ-MORF-C'
        up = data[key][0][-6:]
        cmd = f'{pathon} /mnt/dfc_data2/project/linyusen/project/72_xinquan_rnaloc/github/kmer_diversity_evaluation.py \
                -fq1 {fq1} \
                -fq2 {fq2} \
                -s {os.path.join(save_path, key)} \
                -up_flanking {up} \
                -down_flanking ACGCGT \
                -primer1 {data[key][0]} \
                -primer2 {data[key][1]}'
    print(cmd)
