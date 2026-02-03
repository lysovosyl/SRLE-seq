import math
import torch
import torch.nn as nn
import torch.nn.functional as F

class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len).unsqueeze(1).float()
        div = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0)/d_model))
        pe[:, 0::2] = torch.sin(position * div)
        pe[:, 1::2] = torch.cos(position * div)
        self.register_buffer("pe", pe.unsqueeze(0))

    def forward(self, x):
        return x + self.pe[:, :x.size(1)]


class MultiHeadAttention(nn.Module):
    def __init__(self, d_model, num_heads):
        super().__init__()
        assert d_model % num_heads == 0
        self.d_k = d_model // num_heads
        self.num_heads = num_heads

        self.q = nn.Linear(d_model, d_model)
        self.k = nn.Linear(d_model, d_model)
        self.v = nn.Linear(d_model, d_model)
        self.o = nn.Linear(d_model, d_model)

    def forward(self, q, k, v, mask=None):
        B = q.size(0)

        def split(x):
            return x.view(B, -1, self.num_heads, self.d_k).transpose(1, 2)

        q, k, v = split(self.q(q)), split(self.k(k)), split(self.v(v))
        scores = (q @ k.transpose(-2, -1)) / math.sqrt(self.d_k)

        if mask is not None:
            scores = scores.masked_fill(mask[:, None, None, :] == 0, -1e9)

        attn = F.softmax(scores, dim=-1)
        out = attn @ v
        out = out.transpose(1, 2).contiguous().view(B, -1, self.num_heads * self.d_k)
        return self.o(out)


class EncoderLayer(nn.Module):
    def __init__(self, d_model, heads, d_ff):
        super().__init__()
        self.attn = MultiHeadAttention(d_model, heads)
        self.ff = nn.Sequential(
            nn.Linear(d_model, d_ff),
            nn.ReLU(),
            nn.Linear(d_ff, d_model)
        )
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)

    def forward(self, x, mask):
        x = self.norm1(x + self.attn(x, x, x, mask))
        x = self.norm2(x + self.ff(x))
        return x


class TransformerClassifier(nn.Module):
    def __init__(self, num_classes=2, d_model=256, layers=4, heads=4, d_ff=1024):
        super().__init__()
        self.input_proj = nn.Linear(1, d_model)
        self.pos = PositionalEncoding(d_model)
        self.encoders = nn.ModuleList([
            EncoderLayer(d_model, heads, d_ff) for _ in range(layers)
        ])
        self.cls = nn.Linear(d_model, num_classes)

    def forward(self, x):
        B, L = x.shape
        mask = (x != 0)
        x = self.input_proj(x.unsqueeze(-1))

        # === CLS token ===
        cls_token = self.cls_token.expand(B, -1, -1)  # (B, 1, d_model)
        x = torch.cat([cls_token, x], dim=1)  # (B, L+1, d_model)

        # === mask 同步扩展 ===
        cls_mask = torch.ones(B, 1, device=x.device, dtype=mask.dtype)
        mask = torch.cat([cls_mask, mask], dim=1)  # (B, L+1)


        x = self.pos(x)
        for layer in self.encoders:
            x = layer(x, mask)
        return self.cls(x[:, 0])
