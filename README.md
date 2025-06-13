# 🔬 6-mer Diversity Analysis from RNA-seq Fastq Files

This repository provides a Python script for extracting 6-mer sequences between defined primer regions from paired-end FASTQ files and comparing their relative abundance between nuclear and cytoplasmic RNA fractions. The output includes normalized counts and log2 fold-change values for each 6-mer.

## 📌 Features

- Extract 6-mers between user-defined upstream and downstream primer sequences
- Support for paired-end FASTQ input (nucleus and cytoplasm)
- Normalization to counts per million (CPM)
- Outputs:
  - Raw 6-mer counts (`.count.csv`)
  - Log2 fold change between cytoplasm and nucleus (`.csv`)

---

## 📁 Input Requirements

- Four paired-end FASTQ files:
  - Nuclear sample: `nuc_fq1`, `nuc_fq2`
  - Cytoplasmic sample: `cyto_fq1`, `cyto_fq2`
- Primer sequences:
  - `up_primer_z`: upstream anchor primer
  - `down_primer_z`: downstream anchor primer
- Output directory: `save_path`

---

## 🧪 Installation

Requires Python ≥ 3.6 and the following libraries:

```bash
pip install biopython numpy scipy tqdm
