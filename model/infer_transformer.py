import argparse
import csv
import torch
import numpy as np
import pandas as pd
from model import TransformerClassifier


# ------------------ 读取 CSV（带 seq_id） ------------------
def load_sequences_with_id(csv_file):
    seq_ids = []
    sequences = []

    with open(csv_file) as f:
        r = csv.reader(f)
        header = next(r)

        for row in r:
            seq_ids.append(row[0])
            sequences.append(np.array(row[1:], dtype=float))

    return seq_ids, sequences


# ------------------ 主程序 ------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="input CSV with seq_id")
    parser.add_argument("--model", required=True, help="trained model .pth")
    parser.add_argument("--out", default="prediction.csv")
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")

    # 加载模型
    model = TransformerClassifier().to(device)
    model.load_state_dict(torch.load(args.model, map_location=device))
    model.eval()

    # 加载数据
    seq_ids, sequences = load_sequences_with_id(args.input)

    results = []

    with torch.no_grad():
        for seq_id, seq in zip(seq_ids, sequences):
            x = torch.tensor(seq, dtype=torch.float32).unsqueeze(0).to(device)

            logits = model(x)
            probs = torch.softmax(logits, dim=1)[0].cpu().numpy()
            pred = int(probs.argmax())

            results.append({
                "seq_id": seq_id,
                "length": len(seq),
                "pred_class": pred,
                "prob_0": float(probs[0]),
                "prob_1": float(probs[1]),
            })

    df = pd.DataFrame(results)
    df.to_csv(args.out, index=False)
    print(f"Saved predictions to {args.out}")
    print(df.head())


if __name__ == "__main__":
    main()
