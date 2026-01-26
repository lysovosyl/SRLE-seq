## 📍 Sequence location classify model
### 1.Introduction
This section aims to train a deep learning–based classification model to distinguish the subcellular localization of DNA sequences within cells.

### 2.Input Requirements
The script accepts the following required command-line arguments:

| Argument       | Type   | Description                         | Default |
|----------------|--------|-------------------------------------|---------|
| `--file0`      | `str`  | TSV file containing class 0 samples |         |
| `--file1`      | `str`  | TSV file containing class 1 samples |         |
| `--out`        | `str`  | Output directory for results        |         |
| `--epochs`     | `str`  | Number of training epochs           |         |
| `--batch`      | `str`  | Batch size                          |         |
| `--lr`         | `str`  | Learning rate                       |         |
| `--test_size`  | `int`  | Proportion of data used for testing |         |

### 3.Quick Start
Follow these steps to quickly 
```angular2html
python train_model.py \
  --file0 cyto.csv (class0) \
  --file1 nuc.csv (class1) \
  --out results \
  --epochs 20 \
  --batch 8 \
  --lr 1e-4
```

### 4.Check the Output
After training completes, the output directory will contain:
model.pt – The saved model weights (PyTorch checkpoint).
accuracy.csv – A log table of training progress for each epoch.
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
