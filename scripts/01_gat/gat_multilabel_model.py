#!/usr/bin/env python
"""Three-syndrome multi-label GAT used for herb prioritization.

This public version preserves the recovered model architecture and split logic while
making paths explicit and adding deterministic seed controls. The original prescription
workbooks were not present in the recovered archive, so exact rerunning of the historical
GAT fit requires the original source workbooks. Final reported performance is preserved
in results/reference/gat_performance_table_s4.tsv.
"""
from __future__ import annotations
import argparse
from collections import Counter
from pathlib import Path
import random
import re

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    hamming_loss,
    label_ranking_average_precision_score,
    label_ranking_loss,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split

REPO_ROOT = Path(__file__).resolve().parents[2]
SYNDROMES = ["风热犯肺证", "痰热壅肺证", "痰湿阻肺证"]
COOC_MIN_COUNT = 6  # recovered final setting: co-occurrence >5


def norm_cn(value: object) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    value = str(value).strip()
    value = re.sub(r"[\s\t\r\n]+", "", value)
    return re.sub(r"[，,、;；:：/\\|]+", "；", value)


def split_herb_list(cell: object) -> list[str]:
    text = norm_cn(cell)
    return [x for x in re.split(r"[；,，、\s]+", text) if x] if text else []


def extract_prescriptions(df: pd.DataFrame, known_herbs: set[str]) -> list[list[str]]:
    """Recover prescription-level herb lists from several common worksheet layouts."""
    prescriptions: list[list[str]] = []
    candidate_text_cols = []
    for col in df.columns:
        if df[col].dtype != object:
            continue
        sample = df[col].dropna().astype(str).head(30).tolist()
        score = 0
        for value in sample:
            herbs = [h for h in split_herb_list(value) if h in known_herbs]
            if len(herbs) >= 2:
                score += 1
        if score >= 3:
            candidate_text_cols.append(col)
    if candidate_text_cols:
        col = candidate_text_cols[0]
        for value in df[col].dropna().astype(str):
            herbs = list(dict.fromkeys(h for h in split_herb_list(value) if h in known_herbs))
            if len(herbs) >= 2:
                prescriptions.append(herbs)
        if prescriptions:
            return prescriptions

    cols_as_herbs = [c for c in df.columns if str(c).strip() in known_herbs]
    if len(cols_as_herbs) >= 5:
        for _, row in df[cols_as_herbs].iterrows():
            herbs = [str(c).strip() for c in cols_as_herbs if row.get(c, 0) in [1, "1", True]]
            herbs = list(dict.fromkeys(herbs))
            if len(herbs) >= 2:
                prescriptions.append(herbs)
        if prescriptions:
            return prescriptions

    obj_cols = [c for c in df.columns if df[c].dtype == object]
    for _, row in df[obj_cols].iterrows():
        herbs: list[str] = []
        for col in obj_cols:
            herbs.extend(h for h in split_herb_list(row.get(col, "")) if h in known_herbs)
        herbs = list(dict.fromkeys(herbs))
        if len(herbs) >= 2:
            prescriptions.append(herbs)
    return prescriptions


def read_feature_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path, sheet_name=0)
    else:
        df = pd.read_csv(path, sep="\t")
    herb_col = next((c for c in ["Herb_ChineseName", "中药名"] if c in df.columns), None)
    if herb_col is None:
        raise ValueError("Feature table must contain 'Herb_ChineseName' or '中药名'.")
    return df.rename(columns={herb_col: "Herb_ChineseName"})


def set_all_seeds(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    import torch
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def main() -> None:
    ap = argparse.ArgumentParser(description="Train the recovered three-syndrome multi-label GAT.")
    ap.add_argument("--features", type=Path, default=REPO_ROOT / "data/derived/herb_features_table_s3.tsv")
    ap.add_argument("--wind-heat", type=Path, default=REPO_ROOT / "data/restricted_inputs/syndrome_wind_heat.xlsx")
    ap.add_argument("--tanre-yongfei", type=Path, default=REPO_ROOT / "data/restricted_inputs/syndrome_tanre_yongfei.xlsx")
    ap.add_argument("--phlegm-dampness", type=Path, default=REPO_ROOT / "data/restricted_inputs/syndrome_phlegm_dampness.xlsx")
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "results/generated/gat")
    ap.add_argument("--seed", type=int, default=42, help="Deterministic public-repository seed. Historical weight initialization was not saved.")
    ap.add_argument("--cooc-min-count", type=int, default=COOC_MIN_COUNT)
    ap.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    args = ap.parse_args()

    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
        from torch_geometric.data import Data
        from torch_geometric.nn import GATConv
    except ImportError as exc:
        raise SystemExit("GAT stage requires PyTorch and PyTorch Geometric. See environment/requirements-gat.txt") from exc

    for path in [args.features, args.wind_heat, args.tanre_yongfei, args.phlegm_dampness]:
        if not path.exists():
            raise SystemExit(f"Missing required GAT input: {path}")

    set_all_seeds(args.seed)
    args.out.mkdir(parents=True, exist_ok=True)

    feat_df = read_feature_table(args.features)
    feat_df["Herb_ChineseName"] = feat_df["Herb_ChineseName"].astype(str).str.strip()
    feature_cols = [c for c in feat_df.columns if c != "Herb_ChineseName"]
    numeric = feat_df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    herb_names = feat_df["Herb_ChineseName"].tolist()
    known_herbs = set(herb_names)
    herb2id = {h: i for i, h in enumerate(herb_names)}
    X_np = numeric.to_numpy(dtype=np.float32)

    syndrome_files = {
        "风热犯肺证": args.wind_heat,
        "痰热壅肺证": args.tanre_yongfei,
        "痰湿阻肺证": args.phlegm_dampness,
    }
    Y_np = np.zeros((len(herb_names), len(SYNDROMES)), dtype=np.float32)
    cooc_counter: Counter[tuple[int, int]] = Counter()
    prescription_counts = {}
    for si, syndrome in enumerate(SYNDROMES):
        df = pd.read_excel(syndrome_files[syndrome])
        prescs = extract_prescriptions(df, known_herbs)
        prescription_counts[syndrome] = len(prescs)
        if not prescs:
            raise RuntimeError(f"No prescriptions could be parsed for {syndrome}: {syndrome_files[syndrome]}")
        herbs_in_syndrome: set[str] = set()
        for herbs in prescs:
            herbs_in_syndrome.update(herbs)
            ordered = sorted(set(herbs), key=herb2id.get)
            for i in range(len(ordered)):
                for j in range(i + 1, len(ordered)):
                    a, b = herb2id[ordered[i]], herb2id[ordered[j]]
                    cooc_counter[(a, b)] += 1
        for herb in herbs_in_syndrome:
            Y_np[herb2id[herb], si] = 1.0

    edges = []
    for (a, b), count in cooc_counter.items():
        if count >= args.cooc_min_count:
            edges.extend([(a, b), (b, a)])
    if not edges:
        raise RuntimeError("No co-occurrence edges passed the configured threshold; verify the prescription files.")

    labeled_idx = np.where(Y_np.sum(axis=1) > 0)[0]
    if len(labeled_idx) < 10:
        raise RuntimeError("Fewer than 10 labeled herbs were recovered; verify herb-name matching.")
    train_idx, temp_idx = train_test_split(labeled_idx, test_size=0.4, random_state=42)
    val_idx, test_idx = train_test_split(temp_idx, test_size=0.5, random_state=42)

    X = torch.tensor(X_np, dtype=torch.float32)
    Y = torch.tensor(Y_np, dtype=torch.float32)
    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
    masks = {}
    for name, idxs in [("train", train_idx), ("val", val_idx), ("test", test_idx)]:
        mask = torch.zeros(len(herb_names), dtype=torch.bool)
        mask[idxs] = True
        masks[name] = mask
    data = Data(x=X, edge_index=edge_index, y=Y)
    data.train_mask, data.val_mask, data.test_mask = masks["train"], masks["val"], masks["test"]

    class GATNet(nn.Module):
        def __init__(self, in_dim: int, hidden_dim: int, out_dim: int, heads: int = 8, dropout: float = 0.3):
            super().__init__()
            self.dropout = dropout
            self.gat1 = GATConv(in_dim, hidden_dim, heads=heads, dropout=dropout, concat=True)
            self.gat2 = GATConv(hidden_dim * heads, out_dim, heads=1, dropout=dropout, concat=False)
        def forward(self, x, edges_):
            x = F.dropout(x, p=self.dropout, training=self.training)
            x = F.elu(self.gat1(x, edges_))
            x = F.dropout(x, p=self.dropout, training=self.training)
            return self.gat2(x, edges_)

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    elif args.device == "cuda" and not torch.cuda.is_available():
        raise SystemExit("--device cuda requested but CUDA is unavailable.")
    else:
        device = torch.device(args.device)
    model = GATNet(X.shape[1], 256, len(SYNDROMES), heads=8, dropout=0.3).to(device)
    data = data.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-4)
    criterion = nn.BCEWithLogitsLoss()

    def macro_ap(logits, y_true):
        truth = y_true.detach().cpu().numpy()
        prob = torch.sigmoid(logits).detach().cpu().numpy()
        return float(average_precision_score(truth, prob, average="macro"))

    best_val_ap, best_state, patience_count = -1.0, None, 0
    history = []
    for epoch in range(1, 101):
        model.train(); optimizer.zero_grad()
        logits = model(data.x, data.edge_index)
        loss = criterion(logits[data.train_mask], data.y[data.train_mask])
        loss.backward(); optimizer.step()
        model.eval()
        with torch.no_grad():
            val_logits = model(data.x, data.edge_index)[data.val_mask]
            val_ap = macro_ap(val_logits, data.y[data.val_mask])
        history.append({"epoch": epoch, "loss": float(loss.item()), "val_macro_auprc": val_ap})
        if val_ap > best_val_ap + 1e-4:
            best_val_ap = val_ap
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            patience_count = 0
        else:
            patience_count += 1
        if patience_count >= 10:
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    model = model.to(device); model.eval()
    with torch.no_grad():
        all_logits = model(data.x, data.edge_index)
        probs = torch.sigmoid(all_logits).cpu().numpy()
        test_prob = probs[test_idx]
        test_true = Y_np[test_idx].astype(int)

    pred = pd.DataFrame({"Herb_ChineseName": herb_names})
    for i, syndrome in enumerate(SYNDROMES):
        pred[f"P({syndrome})"] = probs[:, i]
        pred[f"label_{syndrome}"] = Y_np[:, i].astype(int)
    pred.to_csv(args.out / "GAT_syndrome_predictions.tsv", sep="\t", index=False)
    pred.to_excel(args.out / "GAT_syndrome_predictions.xlsx", index=False, sheet_name="predictions")

    metric_rows = []
    for i, syndrome in enumerate(SYNDROMES):
        y, p = test_true[:, i], test_prob[:, i]
        yhat = (p >= 0.5).astype(int)
        auc = roc_auc_score(y, p) if len(np.unique(y)) == 2 else np.nan
        metric_rows.append({
            "Syndrome": syndrome, "Test n": len(y), "Positive labels": int(y.sum()), "Negative labels": int((1-y).sum()),
            "AUROC": auc, "AUPRC": average_precision_score(y, p),
            "Precision (0.50)": precision_score(y, yhat, zero_division=0), "Recall (0.50)": recall_score(y, yhat, zero_division=0),
            "F1 (0.50)": f1_score(y, yhat, zero_division=0), "Accuracy (0.50)": accuracy_score(y, yhat),
        })
    metric_df = pd.DataFrame(metric_rows)
    macro = {"Syndrome": "Macro average", "Test n": len(test_idx), "Positive labels": int(test_true.sum()), "Negative labels": int(test_true.size-test_true.sum())}
    for col in ["AUROC", "AUPRC", "Precision (0.50)", "Recall (0.50)", "F1 (0.50)", "Accuracy (0.50)"]:
        macro[col] = float(metric_df[col].mean())
    metric_df = pd.concat([metric_df, pd.DataFrame([macro])], ignore_index=True)
    metric_df.to_csv(args.out / "heldout_per_syndrome_metrics.tsv", sep="\t", index=False)

    overall = pd.DataFrame([
        {"metric": "LRAP", "value": label_ranking_average_precision_score(test_true, test_prob)},
        {"metric": "Label ranking loss", "value": label_ranking_loss(test_true, test_prob)},
        {"metric": "Micro F1 (0.50)", "value": f1_score(test_true.ravel(), (test_prob >= 0.5).astype(int).ravel(), average="micro")},
        {"metric": "Hamming loss (0.50)", "value": hamming_loss(test_true, (test_prob >= 0.5).astype(int))},
        {"metric": "Exact-match subset accuracy (0.50)", "value": accuracy_score(test_true, (test_prob >= 0.5).astype(int))},
    ])
    overall.to_csv(args.out / "heldout_overall_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(history).to_csv(args.out / "training_history.tsv", sep="\t", index=False)
    split = pd.DataFrame({"Herb_ChineseName": herb_names, "split": "unlabeled"})
    split.loc[train_idx, "split"] = "train"; split.loc[val_idx, "split"] = "validation"; split.loc[test_idx, "split"] = "test"
    split.to_csv(args.out / "node_splits.tsv", sep="\t", index=False)
    graph_summary = pd.DataFrame({
        "item": ["num_herbs", "num_features", "num_edges_directed", "cooc_min_count", "num_labeled_herbs", "seed", "device"] + [f"prescriptions_{s}" for s in SYNDROMES],
        "value": [len(herb_names), X.shape[1], edge_index.shape[1], args.cooc_min_count, len(labeled_idx), args.seed, str(device)] + [prescription_counts[s] for s in SYNDROMES],
    })
    graph_summary.to_csv(args.out / "graph_summary.tsv", sep="\t", index=False)
    torch.save(best_state, args.out / "best_model_state.pt")
    print(f"Saved GAT outputs to {args.out}")
    print(metric_df.to_string(index=False))
    print(overall.to_string(index=False))

if __name__ == "__main__":
    main()
