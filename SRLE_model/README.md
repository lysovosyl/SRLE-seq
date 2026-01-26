# 🧬 Abstract
This project implements a Transformer-based deep learning model for binary classification of numeric sequences.

---
# ⚙️ Functional

## 📚 Deep learning model classify sequence
### 1.Introduction
This framework provides a PyTorch-based sequence classification pipeline using a Transformer architecture.

### 2.Input Requirements

| Argument       | Type     | Description                         | Default |
|----------------|----------|-------------------------------------|---------|
| `--file0`      | `str` | TSV file containing class 0 samples |         |
| `--file1`      | `str` | TSV file containing class 1 samples |         |
| `--out`        | `str` | Output directory for results        |         |
| `--epochs`     | `str` | Number of training epochs           |         |
| `--batch`      | `str` | Batch size                          |         |
| `--lr`         | `str` | Learning rate                       |         |
| `--test_size`  | `int` | Proportion of data used for testing |         |


### 3.Quick Start
Follow these steps to train the model:
```angular2html
python train.py \
  --file0 class0.tsv \
  --file1 class1.tsv \
  --out results \
  --epochs 20 \
  --batch 8 \
  --lr 1e-4
```
### 4.Check the Output
After training completes, the output directory will contain:
```angular2html
model.pt
accuracy.csv
```