#!/usr/bin/env python
"""Apply the final two-category Tanre Yongfei herb-retention rule.

The recovered development script contained a third 'cross-syndrome control' branch.
The submitted final supplementary Table S5 contains only A (label-positive core) and
B (high-probability annotation-supported) categories; cross-syndrome status is retained
only as descriptive metadata. This public script follows that final reporting logic.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TARGET = "痰热壅肺证"
PROB_COL = f"P({TARGET})"
LABEL_COL = f"label_{TARGET}"
TOP_PERCENT = 0.15
KEYWORDS = [
    "清热", "化痰", "祛痰", "燥湿化痰", "温化寒痰", "止咳", "润肺", "清肺", "宣肺", "肃肺",
    "泻肺", "平喘", "降气", "下气", "豁痰", "痰", "喘", "咳", "肺热", "痰热", "痰壅"
]

def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=0)
    return pd.read_csv(path, sep="\t")

def normalize_text(x) -> str:
    return "" if pd.isna(x) else str(x).strip()

def hits(text: str) -> list[str]:
    return [k for k in KEYWORDS if k in normalize_text(text)]

def main() -> None:
    ap = argparse.ArgumentParser(description="Select final Tanre Yongfei-associated herbs from GAT outputs.")
    ap.add_argument("--predictions", type=Path, default=ROOT / "results/generated/gat/GAT_syndrome_predictions.tsv")
    ap.add_argument("--annotations", type=Path, default=ROOT / "data/derived/herb_annotations_table_s2.tsv")
    ap.add_argument("--out", type=Path, default=ROOT / "results/generated/herb_selection/selected_tanre_yongfei_herbs.tsv")
    ap.add_argument("--top-percent", type=float, default=TOP_PERCENT)
    args = ap.parse_args()
    pred, ann = read_table(args.predictions), read_table(args.annotations)
    for df in [pred, ann]:
        if "中药名" in df.columns and "Herb_ChineseName" not in df.columns:
            df.rename(columns={"中药名": "Herb_ChineseName"}, inplace=True)
    if "Herb_ChineseName" not in pred.columns or "Herb_ChineseName" not in ann.columns:
        raise ValueError("Both tables must contain Herb_ChineseName (or 中药名).")
    if PROB_COL not in pred.columns or LABEL_COL not in pred.columns:
        raise ValueError(f"Prediction table must contain {PROB_COL!r} and {LABEL_COL!r}.")

    # Map standardized supplementary annotation headers onto the names used in the recovered rule.
    aliases = {
        "Nature, Flavor and Toxicity": "性味", "Meridian Tropism": "归经",
        "Traditional Actions": "功能", "Traditional Indications": "主治",
    }
    ann = ann.rename(columns={k: v for k, v in aliases.items() if k in ann.columns})
    df = pred.merge(ann, on="Herb_ChineseName", how="left", suffixes=("", "_annotation"))
    for col in ["性味", "归经", "功能", "主治"]:
        if col not in df.columns:
            df[col] = ""
        df[col] = df[col].map(normalize_text)

    other_label_cols = [c for c in df.columns if c.startswith("label_") and c != LABEL_COL]
    df["is_shared_other_syndrome"] = df[other_label_cols].fillna(0).sum(axis=1).ge(1) if other_label_cols else False
    df["文本合并"] = df["功能"] + "；" + df["主治"]
    df["命中关键词"] = df["文本合并"].map(lambda x: "、".join(hits(x)))
    df["功效匹配"] = df["命中关键词"].str.len().gt(0)
    df["功能主治完整"] = df["功能"].str.len().gt(0) & df["主治"].str.len().gt(0)
    threshold = float(df[PROB_COL].quantile(1 - args.top_percent))
    df["高分候选"] = df[PROB_COL].ge(threshold)

    core = df[df[LABEL_COL].eq(1)].copy(); core["入选层级"] = "A_原始标签核心药"
    candidate = df[df[LABEL_COL].eq(0) & df["高分候选"] & df["功效匹配"] & df["功能主治完整"]].copy()
    candidate["入选层级"] = "B_GAT高分未标注候选药"
    result = pd.concat([core, candidate], ignore_index=True).drop_duplicates("Herb_ChineseName")
    result["_order"] = result["入选层级"].map({"A_原始标签核心药": 1, "B_GAT高分未标注候选药": 2})
    result = result.sort_values(["_order", PROB_COL], ascending=[True, False]).drop(columns="_order")
    result["Retention rationale"] = result["入选层级"].map({
        "A_原始标签核心药": "Original label-positive core herb",
        "B_GAT高分未标注候选药": "High-probability unlabeled candidate retained after annotation-guided screening",
    })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)
    summary = pd.DataFrame([
        {"metric": "total herbs", "value": len(df)},
        {"metric": "Tanre Yongfei label-positive herbs", "value": int(df[LABEL_COL].eq(1).sum())},
        {"metric": "top-probability threshold", "value": threshold},
        {"metric": "A core herbs", "value": len(core)},
        {"metric": "B retained candidates", "value": len(candidate)},
        {"metric": "final retained herbs", "value": len(result)},
    ])
    summary.to_csv(args.out.with_name("selection_summary.tsv"), sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"Saved {args.out}")

if __name__ == "__main__":
    main()
