# Recovered STRING v12 conversion script; public CLI entry point restored.
import argparse
import gzip
from pathlib import Path
import pandas as pd


def read_alias_map(alias_gz: str, preferred_sources=None):
    """
    从 STRING aliases 文件构建：protein_id -> gene_symbol 的映射
    文件格式（tsv）通常是：protein_id, alias, source
    一个 protein 会有多个 alias，这里按 source 优先级挑一个“最像基因符号”的。
    """
    if preferred_sources is None:
        # 不同版本 STRING source 命名略有差异；这里给一个常用优先级列表
        preferred_sources = [
            "Ensembl_HGNC_symbol",
            "HGNC",
            "Ensembl_gene",
            "Ensembl",
            "Gene_Name",
            "gene name",
            "BioMart_HUGO",
            "BLAST_UniProt_GN",
        ]

    # 读取 aliases
    alias_gz = str(alias_gz)
    df = pd.read_csv(
        alias_gz,
        sep="\t",
        compression="gzip",
        header=0,
        dtype=str
    )

    # 兼容列名（STRING 文件列名一般为：#string_protein_id / string_protein_id）
    cols = [c.strip() for c in df.columns]
    df.columns = cols

    # 猜列名
    prot_col = None
    for c in df.columns:
        if "protein" in c and "id" in c:
            prot_col = c
            break
    if prot_col is None:
        prot_col = df.columns[0]

    alias_col = None
    for c in df.columns:
        if c.lower() in ["alias", "preferred_name", "preferredname"]:
            alias_col = c
            break
    if alias_col is None:
        alias_col = df.columns[1]

    source_col = None
    for c in df.columns:
        if "source" in c.lower():
            source_col = c
            break
    if source_col is None:
        source_col = df.columns[2]

    df = df[[prot_col, alias_col, source_col]].rename(
        columns={prot_col: "protein_id", alias_col: "alias", source_col: "source"}
    )

    # 清洗
    df["alias"] = df["alias"].astype(str).str.strip()
    df["source"] = df["source"].astype(str).str.strip()

    # 先挑 preferred_sources 中出现的
    df["source_rank"] = df["source"].apply(
        lambda s: preferred_sources.index(s) if s in preferred_sources else 10**9
    )
    preferred = df.sort_values(["protein_id", "source_rank"]).drop_duplicates("protein_id")

    # 对还没映射到的 protein，fallback：用“看起来像基因符号”的 alias（全大写字母数字/短横线）
    mapped = dict(zip(preferred["protein_id"], preferred["alias"]))

    # 补全：对没映射的 protein 用启发式挑一个 alias
    remain = df[~df["protein_id"].isin(mapped.keys())].copy()
    if len(remain) > 0:
        # 简单启发式：优先全大写+数字+短横线，且长度 2~15
        def score_alias(a: str) -> int:
            if not isinstance(a, str):
                return 10**9
            a = a.strip()
            if len(a) < 2 or len(a) > 20:
                return 10**8
            # 排除明显不是基因名的（包含空格、点号太多等）
            if " " in a:
                return 10**8
            # 基因符号常见特征
            ok = all(ch.isalnum() or ch in ["-", "_"] for ch in a)
            if not ok:
                return 10**7
            # 更偏好大写（但不强制）
            upper_bonus = 0 if a.upper() == a else 100
            return upper_bonus + len(a)

        remain["alias_score"] = remain["alias"].apply(score_alias)
        remain2 = remain.sort_values(["protein_id", "alias_score"]).drop_duplicates("protein_id")
        for pid, alias in zip(remain2["protein_id"], remain2["alias"]):
            mapped[pid] = alias

    return mapped


def build_edges(links_gz: str, protein2gene: dict, min_score: int, out_path: str):
    """
    从 protein links detailed 构建 gene-gene edges（用 combined_score）
    输出去重（同一 gene 对取最大 score），并规范 geneA < geneB（字典序）方便后续处理。
    """
    links_gz = str(links_gz)
    out_path = str(out_path)

    # 读取 links detailed
    df = pd.read_csv(
        links_gz,
        sep=r"\s+",
        compression="gzip",
        header=0,
        dtype={"protein1": str, "protein2": str}
    )

    # 兼容列名：通常 protein1 protein2 combined_score ...
    cols = [c.strip() for c in df.columns]
    df.columns = cols

    if "combined_score" not in df.columns:
        raise ValueError("未在 links.detailed 文件中找到 combined_score 列。请检查文件版本/列名。")
    if "protein1" not in df.columns or "protein2" not in df.columns:
        # 有些版本列名会带 # 号
        p1 = next((c for c in df.columns if "protein1" in c), None)
        p2 = next((c for c in df.columns if "protein2" in c), None)
        if p1 and p2:
            df = df.rename(columns={p1: "protein1", p2: "protein2"})
        else:
            raise ValueError("未找到 protein1/protein2 列，请检查 links 文件列名。")

    df["combined_score"] = pd.to_numeric(df["combined_score"], errors="coerce")
    df = df.dropna(subset=["combined_score"])
    df = df[df["combined_score"] >= min_score].copy()

    # protein -> gene
    df["gene1"] = df["protein1"].map(protein2gene)
    df["gene2"] = df["protein2"].map(protein2gene)
    df = df.dropna(subset=["gene1", "gene2"])

    # 去掉自环
    df = df[df["gene1"] != df["gene2"]].copy()

    # 规范顺序 geneA < geneB
    geneA = df[["gene1", "gene2"]].min(axis=1)
    geneB = df[["gene1", "gene2"]].max(axis=1)
    df["geneA"] = geneA
    df["geneB"] = geneB

    # 同一 gene 对多条边：取最大 score
    df = df.groupby(["geneA", "geneB"], as_index=False)["combined_score"].max()
    df = df.rename(columns={"combined_score": "score"})

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)

    return df


def main():
    parser = argparse.ArgumentParser(description="Build gene-gene PPI edges from STRING v12 protein links + aliases.")
    parser.add_argument("--links", required=True, help="Path to 9606.protein.links.detailed.v12.0.txt.gz")
    parser.add_argument("--aliases", required=True, help="Path to 9606.protein.aliases.v12.0.txt.gz")
    parser.add_argument("--min_score", type=int, default=700, help="Min combined_score threshold (e.g., 400/700/900)")
    parser.add_argument("--out", default="string_ppi_edges.tsv", help="Output TSV path")
    args = parser.parse_args()

    print("Reading aliases:", args.aliases)
    protein2gene = read_alias_map(args.aliases)
    print("protein2gene mapped:", len(protein2gene))

    print("Building edges from links:", args.links)
    edges = build_edges(args.links, protein2gene, args.min_score, args.out)
    print("Saved:", args.out)
    print("Edges:", edges.shape[0])
    print(edges.head(5).to_string(index=False))


if __name__ == "__main__":
    main()
