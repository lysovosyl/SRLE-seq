# 🧬 Abstract
RNA molecules localize to specific subcellular compartments to perform their functions. While most mRNAs are translated at the endoplasmic reticulum, many lncRNAs act within organelles. To better understand RNA localization mechanisms, we developed SRLE-seq (Screening RNA Localization Elements by Sequencing), a high-throughput method to identify functional localization elements. 

---
# 💻 Installation

Requires Python ≥ 3.8 and the following libraries:
```bash
pip install biopython numpy scipy tqdm
```
---
# ⚙️ Functional

## 📚 Kmer-based library diversity validation
### 1.Introduction
This section is designed to validate whether the constructed Kmer-based plasmid library is suitable for downstream experiments. The distribution of plasmid types corresponding to each k-mer should be relatively uniform. No single k-mer type should dominate the library in abundance, as this would indicate amplification bias or cloning artifacts.


### 2.Input Requirements
The script accepts the following required command-line arguments:

| Argument         | Type     | Description                                              | Default |
|------------------|----------|----------------------------------------------------------|---------|
| `-fq1`           | `str`    | Path to the nuclear sample R1 (forward reads) FASTQ file |         |
| `-fq2`           | `str`    | Path to the nuclear sample R2 (reverse reads) FASTQ file | None    |
| `-s`             | `str`    | Output directory for saving result figures and tables    |         |
| `-up_flanking`   | `str`    | Sequence of the upstream (5') anchor primer              |         |
| `-down_flanking` | `str`    | Sequence of the downstream (3') anchor primer            |         |
| `-kmer_length`   | `str`    | Length of k-mers to count                                |         |
> ⚠️ For paired-end reads with Paired-End sequences, at least one of the primer-binding regions should be longer than 150 bp. Otherwise, the reads should be merged using PEAR before further processing. 

### 3.Quick Start
Follow these steps to quickly analyse library diversity:
```angular2html
python kmer_diversity_validation.py \
  -fq1 R1.fastq \
  -fq2 R2.fastq \
  -s ./ \
  -up_flanking ATTATGAT \
  -down_flanking GCTTAGTG \
  -kmer_length 6 \
```
### 4.Check the Output

After completion, the output directory will contain:

kmer_diversity.validation.tsv – records the raw count and CPM (counts per million) for each k-mer.

---

## 📚 Fragment-based library diversity validation
### 1.Introduction
This section is designed to validate whether the constructed fragment-based plasmid library is suitable for downstream experiments. The inserted fragments should collectively achieve broad and uniform coverage across the full-length source sequence. Ideally, all regions of the original sequence should be represented within the library, ensuring comprehensive functional screening.

### 2.Input Requirements
The script accepts the following required command-line arguments:

| Argument             | Type  | Description                                                | Default |
|----------------------|-------|------------------------------------------------------------|---------|
| `-fq1`               | `str` | Path to the nuclear sample R1 (forward reads) FASTQ file   |         |
| `-fq2`               | `str` | Path to the nuclear sample R2 (reverse reads) FASTQ file   |         |
| `-s`                 | `str` | Output directory for saving result figures and tables      |         |
| `-up_flanking`       | `str` | Sequence of the upstream (5') anchor primer                |         |
| `-down_flanking`     | `str` | Sequence of the downstream (3') anchor primer              |         |
| `-gene_chr`          | `str` | Chromosome name of the gene locus (e.g., `chr11`)          |         |
| `-gene_region_start` | `int` | Start coordinate of the target gene region (1-based)       |         |
| `-gene_region_end`   | `int` | End coordinate of the target gene region (1-based)         |         |
| `-config`            | `str` | Path to a configuration file containing parameter presets  |         |
| `-thread`            | `int` | Number of threads to use for parallel processing           |         |
> ⚠️ For paired-end reads with Paired-End sequences, at least one of the primer-binding regions should be longer than 150 bp. Otherwise, the reads should be merged using PEAR before further processing. 

### 3.Quick Start
Follow these steps to quickly analyse library diversity:
```angular2html
python sequence_diversity_validation.py \
  -fq1 R1.fastq \
  -fq2 R2.fastq \
  -s ./ \
  -up_flanking ATTATGAT \
  -down_flanking GCTTAGTG \
  -gene_chr chr11 \
  -gene_region_start 65497688
  -gene_region_end 65506516
  -config ./config.yaml
  -thread 64
```
### 4.Check the Output

After completion, the output directory will contain:

dna_fragment.coverage.bed - BED-format file showing per-base coverage of extracted fragments across the specified gene region.

dna_fragment.coverage.png - Coverage plot visualizing the distribution of extracted fragments.

---
## 📍 Library location analysis
### 1.Introduction
This section aims to investigate the subcellular localization preferences of individual fragments following plasmid transfection into cells. As a prerequisite, quality control (QC) is performed on the post-transfection library, which involves two key components:

**For the k-mer library:** The presence and abundance of each k-mer in the transfected cells are assessed. To pass QC, a k-mer must exhibit a count-per-million (CPM) value greater than 1 in at least one of the Total, Nuclear (Nuc), or Cytoplasmic (Cyto) compartments, indicating successful transfection and sufficient expression.

**For the sequence fragment library:** Each fragment is required to have a minimum sequencing coverage of 100× to ensure reliable quantification.

Following quality control, the abundance, CPM, and nuclear enrichment score (NES) of each k-mer or sequence fragment are calculated to facilitate subsequent analysis of their subcellular distribution patterns.

### 2.Input Requirements
The script accepts the following required command-line arguments:

| Argument         | Type  | Description                                                 | Default |
|------------------|-------|-------------------------------------------------------------|---------|
| `-nuc_fq1_file`  | `str` | Path to the nuclear sample R1 (forward reads) FASTQ file    |         |
| `-nuc_fq2_file`  | `str` | Path to the nuclear sample R2 (reverse reads) FASTQ file    |         |
| `-cyto_fq1_file` | `str` | Path to the cytoplasm sample R1 (forward reads) FASTQ file  |         |
| `-cyto_fq2_file` | `str` | Path to the cytoplasm sample R2 (forward reads) FASTQ file  |         |
| `-s`             | `str` | Output directory for saving result figures and tables       |         |
| `-up_flanking`   | `str` | Sequence of the upstream (5') anchor primer                 |         |
| `-down_flanking` | `str` | Sequence of the downstream (3') anchor primer               |         |
| `-mode`          | `str` | Processing mode: `kmer_complexity` or `fragment_complexity` |         |
| `-kmer_length`   | `int` | Length of k-mers to count                                   |         |
> ⚠️ For paired-end reads with Paired-End sequences, at least one of the primer-binding regions should be longer than 150 bp. Otherwise, the reads should be merged using PEAR before further processing. 

### 3.Quick Start
Follow these steps to quickly analyse fragment location analysis:
```angular2html
python library_location_analysis.py \
  -nuc_fq1_file nuc.R1.fastq \
  -nuc_fq2_file nuc.R2.fastq \
  -cyto_fq1_file cyto.R1.fastq \
  -cyto_fq2_file cyto.R2.fastq \
  -s ./ \
  -up_flanking ATTATGAT \
  -down_flanking GCTTAGTG \
  -mode kmer_complexity \
  -kmer_length 6 \
```

### 4.Check the Output

After completion, the output directory will contain:

cyto.kmer.count.csv - the counts of each kmer which split from cytoplasm 

nuc.kmer.count.csv  - the counts of each kmer which split from nuclear 

nuc.cyto.kmer.count.png  - Scatter plot comparing the abundance of each k-mer in the cytoplasm (x-axis) versus the nucleus (y-axis).

kmer.info.csv - Summary table reporting the log₂ fold change (cytoplasm vs. nucleus) and Nuclear Retention Score (NES) for each k-mer.


---

# 📝 Please Cite

If you use this script or parts of it in your research or project, please cite the repository or acknowledge the author appropriately. A suggested citation format:

> [Zeng. x. (2025). *SRLE-seq k-mer comparison and visualization pipeline*](https://doi.org/10.5281/zenodo.1234567)

---

# 👨‍💻 Maintainer

**[Lin Yusen](https://github.com/lysovosyl)**  
Email: `linyusen1688@gmail.com`  

Feel free to open an issue or pull request for improvements or bug fixes.

---

# 🤝 Contributors

> Current contributors:
> - Yusen Lin
> - Xingqaun Zeng
> - Jiajian Zhou
---

# 📄 License

This project is licensed under the [MIT License](LICENSE.txt).  
You are free to use, modify, and distribute it with attribution.
