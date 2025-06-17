# 🔬 Abstract

This script performs **k-mer extraction, counting, normalization, and comparison** between nuclear and cytoplasmic RNA-seq paired-end reads, specifically designed for the analysis of **SRLE-seq (Spatially Resolved Ligation-based Expression sequencing)** data.

It focuses on analyzing short sequence fragments located between user-defined **upstream and downstream primers**, allowing precise quantification of ligation-based tags captured during SRLE-seq. This facilitates in-depth profiling of subcellular transcript distribution.

## 🧪 Installation

Requires Python ≥ 3.6 and the following libraries:
```bash
pip install biopython numpy scipy tqdm
```

## 📁 Input Requirements
The script accepts the following required command-line arguments:

| Argument             | Type     | Description                                                                 |
|----------------------|----------|-----------------------------------------------------------------------------|
| `-nuc_fq1_file`      | `str`    | Path to the nuclear sample R1 (forward reads) FASTQ file                    |
| `-nuc_fq2_file`      | `str`    | Path to the nuclear sample R2 (reverse reads) FASTQ file                    |
| `-cyto_fq1_file`     | `str`    | Path to the cytoplasmic sample R1 (forward reads) FASTQ file                |
| `-cyto_fq2_file`     | `str`    | Path to the cytoplasmic sample R2 (reverse reads) FASTQ file                |
| `-save_path`         | `str`    | Output directory for saving result figures and tables                       |
| `-up_primer_z`       | `str`    | Sequence of the upstream (5') anchor primer                                 |
| `-down_primer_z`     | `str`    | Sequence of the downstream (3') anchor primer                               |
| `-mode`              | `str`    | Processing mode: `library_complexity` or `cell_complexity`                  |
| `-kmer_length`       | `int`    | Length of k-mers to extract and count                                       |

> ⚠️ All arguments are **required**.

## ⚡ Quick Start
Follow these steps to quickly run the SRLE-seq k-mer comparison pipeline:
```angular2html
python analyse_kmer_library.py \
  -nuc_fq1_file nuc_R1.fastq \
  -nuc_fq2_file nuc_R2.fastq \
  -cyto_fq1_file cyto_R1.fastq \
  -cyto_fq2_file cyto_R2.fastq \
  -up_primer_z ATTATGAT \
  -down_primer_z GCTTAGTG \
  -kmer_length 6 \
  -mode library_complexity \
  -save_path ./results/
```
### Check the Output
After completion, the output directory will contain:

nuc.kmer.count.csv: Nuclear k-mer counts and CPM values

cyto.kmer.count.csv: Cytoplasmic k-mer counts and CPM values

kmer.log2fc.csv: Log2 fold-change between nuclear and cytoplasmic

nuc.cyto.kmer.count.png: Scatter plot of k-mer abundance comparison

---

## 📝 Please Cite

If you use this script or parts of it in your research or project, please cite the repository or acknowledge the author appropriately. A suggested citation format:

> Zeng. x. (2025). *SRLE-seq k-mer comparison and visualization pipeline*. GitHub repository: [https://github.com/your_github_name/KmerCountCompare](https://github.com/your_github_name/KmerCountCompare)

---

## 👨‍💻 Maintainer

**Lin Yusen**  
Email: `linyusen1688@gmail.com`  
GitHub: [lysovosyl](https://github.com/lysovosyl)

Feel free to open an issue or pull request for improvements or bug fixes.

---

## 🤝 Contributors

> Current contributors:
- Yusen Lin
- Xingqaun Zeng
- Jiajian Zhou
---

## 📄 License

This project is licensed under the [MIT License](LICENSE.txt).  
You are free to use, modify, and distribute it with attribution.
