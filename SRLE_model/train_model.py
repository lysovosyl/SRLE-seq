import argparse, os, csv
import numpy as np
import torch
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from model import TransformerClassifier


# ---------------- 数据读取 ----------------
def load_tsv(file):
    x = []
    with open(file) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            temp = row[1:]
            temp = [float(i) for i in temp]
            x.append(temp)

    return x


# ---------------- Dataset ----------------
class SeqDataset(Dataset):
    def __init__(self, X, Y):
        self.X, self.Y = X, Y

    def __len__(self):
        return len(self.X)

    def __getitem__(self, i):
        return torch.tensor(self.X[i], dtype=torch.float32), torch.tensor(self.Y[i], dtype=torch.long)


def collate_fn(batch):
    seqs, labels = zip(*batch)
    max_len = max(len(s) for s in seqs)

    padded = torch.zeros(len(seqs), max_len, dtype=torch.float32)
    for i, s in enumerate(seqs):
        padded[i, :len(s)] = s

    labels = torch.stack(labels)
    return padded, labels


# ---------------- 主程序 ----------------
parser = argparse.ArgumentParser()
parser.add_argument("--file0", required=True, help="tsv for class 0")
parser.add_argument("--file1", required=True, help="tsv for class 1")
parser.add_argument("--out", default="results", help="output folder")
parser.add_argument("--epochs", type=int, default=20)
parser.add_argument("--batch", type=int, default=8)
parser.add_argument("--lr", type=float, default=1e-4)
parser.add_argument("--test_size", type=float, default=0.2)
args = parser.parse_args()

os.makedirs(args.out, exist_ok=True)

# 读取数据
x, y = [], []
for i in load_tsv(args.file0):
    x.append(i)
    y.append(0)
for i in load_tsv(args.file1):
    x.append(i)
    y.append(1)

# x.append(load_tsv(args.file0)); y.append(0)
# x.append(load_tsv(args.file1)); y.append(1)

# 训练/测试切分
train_x, test_x, train_y, test_y = train_test_split(
    x, y, test_size=args.test_size, random_state=42, shuffle=True
)

# DataLoader
train_loader = DataLoader(
    SeqDataset(train_x, train_y),
    batch_size=args.batch,
    shuffle=True,
    collate_fn=collate_fn
)
test_loader = DataLoader(
    SeqDataset(test_x, test_y),
    batch_size=args.batch,
    shuffle=False,
    collate_fn=collate_fn
)

# 模型/优化器/损失
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = TransformerClassifier().to(device)
opt = torch.optim.Adam(model.parameters(), lr=args.lr)
loss_fn = torch.nn.CrossEntropyLoss()

records = []

# 训练循环
for epoch in range(args.epochs):
    model.train()
    train_loss = 0
    for x_batch, y_batch in train_loader:
        x_batch, y_batch = x_batch.to(device), y_batch.to(device)

        opt.zero_grad()
        logits = model(x_batch)
        loss = loss_fn(logits, y_batch)
        loss.backward()
        opt.step()

        train_loss += loss.item()

    # 测试
    model.eval()
    correct = 0
    total = 0
    with torch.no_grad():
        for x_batch, y_batch in test_loader:
            x_batch, y_batch = x_batch.to(device), y_batch.to(device)
            logits = model(x_batch)
            pred = logits.argmax(dim=1)
            correct += (pred == y_batch).sum().item()
            total += y_batch.size(0)

    acc = correct / total
    avg_loss = train_loss / len(train_loader)

    records.append({
        "epoch": epoch + 1,
        "train_loss": avg_loss,
        "test_acc": acc
    })

    print(f"Epoch {epoch+1} | Loss {avg_loss:.4f} | Test Acc {acc:.4f}")

# 保存模型 & 结果
torch.save(model.state_dict(), os.path.join(args.out, "model.pt"))
pd.DataFrame(records).to_csv(os.path.join(args.out, "accuracy.csv"), index=False)

print("Training finished.")
print(f"Saved model to {os.path.join(args.out, 'model.pt')}")
print(f"Saved accuracy to {os.path.join(args.out, 'accuracy.csv')}")
