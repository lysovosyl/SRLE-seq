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

---
## 📍 Sequence Prediction
### 1.Introduction
This script performs inference using a pre-trained sequence classifier.

### 2.Input Requirements
The script accepts the following required command-line arguments:

| Argument       | Type    | Description                                                            | Default |
|----------------|---------|------------------------------------------------------------------------|---------|
| `--input`      | `str`   | Input CSV file containing sequences with seq_id in the first column    |         |
| `--model`      | `str`   | Pre-trained model file path (.pt / .pth)                               |         |
| `--out`        | `str`   | Output CSV file for prediction results                                 |         |
| `--device`     | `str`   | Device to run prediction (cuda or cpu)                                 |         |


### 3.Quick Start
Run inference using the following command:
```angular2html
python infer_sequence.py \
  --input sequence.tsv \
  --model model.pt \
  --out results.csv \
  --device cuda \
```

### 4.Check the Output
After the script finishes, the output file will contain:
results.csv - A table of predict location for each reads.
---
```
