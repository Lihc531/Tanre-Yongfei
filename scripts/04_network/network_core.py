"""Core network analysis functions for the final manuscript.

The equations and statistical definitions match the recovered final scripts.
The implementation uses SciPy sparse matrices for speed. RWR uses the same
row-normalized weighted adjacency convention as the recovered code:
    p_new = alpha * p0 + (1 - alpha) * W @ p
Final alpha = 0.5 (author-confirmed).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, diags
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr


def load_cap_vector(path):
    df = pd.read_csv(path, sep="\t")
    df = df[["gene_symbol", "D_meta"]].dropna().copy()
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.strip().str.upper()
    df["D_meta"] = pd.to_numeric(df["D_meta"], errors="coerce")
    return df.dropna().drop_duplicates("gene_symbol")


def load_targets(path):
    with open(path, "r", encoding="utf-8") as f:
        return sorted({x.strip().upper() for x in f if x.strip()})


def load_ppi_edges(path):
    e = pd.read_csv(path, sep="\t", compression="infer")
    e["geneA"] = e["geneA"].astype(str).str.strip().str.upper()
    e["geneB"] = e["geneB"].astype(str).str.strip().str.upper()
    e["score"] = pd.to_numeric(e["score"], errors="coerce").fillna(1.0)
    return e.dropna(subset=["geneA", "geneB"]).query("geneA != geneB").copy()


def build_graph(ppi_edges, node_list):
    genes = sorted(set(node_list))
    idx = {g: i for i, g in enumerate(genes)}
    sub = ppi_edges[ppi_edges["geneA"].isin(idx) & ppi_edges["geneB"].isin(idx)]
    i = sub["geneA"].map(idx).to_numpy(dtype=int)
    j = sub["geneB"].map(idx).to_numpy(dtype=int)
    w = sub["score"].to_numpy(dtype=float)
    W_raw = coo_matrix((np.r_[w, w], (np.r_[i, j], np.r_[j, i])), shape=(len(genes), len(genes))).tocsr()
    rs = np.asarray(W_raw.sum(axis=1)).ravel()
    inv = np.divide(1.0, rs, out=np.zeros_like(rs, dtype=float), where=rs != 0)
    W = (diags(inv) @ W_raw).tocsr()
    return genes, idx, W


def make_seed_vector(idx, targets):
    p0 = np.zeros(len(idx), dtype=float)
    mapped = [g for g in targets if g in idx]
    if not mapped:
        raise ValueError("No targets mapped into the CAP-PPI graph.")
    p0[[idx[g] for g in mapped]] = 1.0
    return p0 / p0.sum(), mapped


def rwr(W, p0, alpha=0.5, max_iter=100, tol=1e-8):
    p = p0.copy()
    for _ in range(max_iter):
        p_new = alpha * p0 + (1.0 - alpha) * W.dot(p)
        if np.linalg.norm(p_new - p, ord=1) < tol:
            return p_new
        p = p_new
    return p


def fisher_enrichment(universe, top_set, sig_set):
    universe = set(universe); top_set = set(top_set) & universe; sig_set = set(sig_set) & universe
    a = len(top_set & sig_set); b = len(top_set - sig_set); c = len(sig_set - top_set); d = len(universe - top_set - sig_set)
    oddsratio, p = fisher_exact([[a,b],[c,d]], alternative="greater")
    return oddsratio, p, a


def safe_spearman(x, y):
    x=np.asarray(x,float); y=np.asarray(y,float); m=np.isfinite(x)&np.isfinite(y)
    if m.sum()<3: return np.nan, np.nan
    return spearmanr(x[m],y[m])


def permutation_proximity(W, target_size, test_idx, bg_idx, alpha=0.5, max_iter=100, tol=1e-8,
                          n_perm=1000, seed=42, block_size=128):
    """Permutation null for proximity using blockwise sparse RWR.

    Seed sets are drawn in exactly the same order as the recovered sequential implementation
    (`numpy.random.default_rng(seed).choice(..., replace=False)`). Multiple seed vectors are
    propagated together to reduce runtime on larger STRING graphs. Each column is stopped
    independently once its L1 change falls below `tol`.
    """
    rng = np.random.default_rng(seed)
    n = W.shape[0]
    out = np.empty(n_perm, dtype=float)
    for start in range(0, n_perm, block_size):
        b = min(block_size, n_perm - start)
        P0 = np.zeros((n, b), dtype=float)
        for j in range(b):
            P0[rng.choice(n, size=target_size, replace=False), j] = 1.0 / target_size
        P = P0.copy()
        active = np.ones(b, dtype=bool)
        for _ in range(max_iter):
            cols = np.flatnonzero(active)
            if len(cols) == 0:
                break
            P_new = alpha * P0[:, cols] + (1.0 - alpha) * W.dot(P[:, cols])
            delta = np.sum(np.abs(P_new - P[:, cols]), axis=0)
            P[:, cols] = P_new
            active[cols[delta < tol]] = False
        out[start:start+b] = P[test_idx, :].mean(axis=0) - P[bg_idx, :].mean(axis=0)
    return out


def analyze_target_space(cap_d_meta_path, ppi_edge_path, target_genes_path, mode="diffusion", alpha=0.5,
                         max_iter=100, tol=1e-8, n_perm=1000, random_seed=42, up_thr=1.0, down_thr=-1.0,
                         top_n=200, return_details=False, null_up_override=None, null_down_override=None):
    cap=load_cap_vector(cap_d_meta_path); targets=load_targets(target_genes_path); ppi=load_ppi_edges(ppi_edge_path)
    ppi_nodes=set(ppi.geneA)|set(ppi.geneB); genes0=sorted(set(cap.gene_symbol)&ppi_nodes)
    genes,idx,W=build_graph(ppi,genes0)
    cap_sub=cap[cap.gene_symbol.isin(genes)].drop_duplicates('gene_symbol').set_index('gene_symbol').loc[genes].reset_index()
    cap_up=set(cap_sub.loc[cap_sub.D_meta>=up_thr,'gene_symbol']); cap_down=set(cap_sub.loc[cap_sub.D_meta<=down_thr,'gene_symbol'])
    p0,mapped=make_seed_vector(idx,targets)
    score=rwr(W,p0,alpha,max_iter,tol) if mode=='diffusion' else p0.copy()
    if mode not in {'diffusion','raw'}: raise ValueError('mode must be diffusion or raw')
    df=pd.DataFrame({'gene_symbol':genes,'D_meta':cap_sub.D_meta.to_numpy(),'score':score})
    up=df.gene_symbol.isin(cap_up); down=df.gene_symbol.isin(cap_down); bg=~df.gene_symbol.isin(cap_up|cap_down)
    u_up,p_up=mannwhitneyu(df.loc[up,'score'],df.loc[bg,'score'],alternative='greater')
    u_dn,p_dn=mannwhitneyu(df.loc[down,'score'],df.loc[bg,'score'],alternative='greater')
    up_idx=np.array([idx[g] for g in cap_up],int); dn_idx=np.array([idx[g] for g in cap_down],int); bg_idx=np.array([idx[g] for g in df.loc[bg,'gene_symbol']],int)
    obs_up=float(score[up_idx].mean()-score[bg_idx].mean()); obs_dn=float(score[dn_idx].mean()-score[bg_idx].mean())
    # Preserve recovered ablation definition: permutation null is diffusion-based even when mode='raw'.
    # Optional overrides allow identical null distributions to be reused for configurations that
    # share the same graph and target-space size, reducing runtime without changing results.
    if null_up_override is None:
        null_up=permutation_proximity(W,len(mapped),up_idx,bg_idx,alpha,max_iter,tol,n_perm,random_seed)
    else:
        null_up=np.asarray(null_up_override,dtype=float)
    if null_down_override is None:
        null_dn=permutation_proximity(W,len(mapped),dn_idx,bg_idx,alpha,max_iter,tol,n_perm,random_seed+1)
    else:
        null_dn=np.asarray(null_down_override,dtype=float)
    if len(null_up)==0 or len(null_dn)==0:
        raise ValueError("Permutation null arrays must contain at least one value.")
    emp_up=(np.sum(null_up>=obs_up)+1)/(len(null_up)+1); emp_dn=(np.sum(null_dn>=obs_dn)+1)/(len(null_dn)+1)
    top=set(df.sort_values('score',ascending=False).head(top_n).gene_symbol)
    or_up,fp_up,ov_up=fisher_enrichment(genes,top,cap_up); or_dn,fp_dn,ov_dn=fisher_enrichment(genes,top,cap_down)
    absD=df.D_meta.abs(); posD=df.D_meta.clip(lower=0)
    rho_abs,pr_abs=safe_spearman(df.score,absD); rho_pos,pr_pos=safe_spearman(df.score,posD)
    metrics=dict(mode=mode,ppi_edges=len(ppi),ppi_nodes=len(ppi_nodes),cap_genes_in_ppi=len(genes),targets_total=len(targets),
                 targets_in_ppi=len(set(targets)&ppi_nodes),targets_mapped_into_cap_ppi=len(mapped),cap_up_genes=len(cap_up),cap_down_genes=len(cap_down),
                 CAP_up_vs_bg_U=float(u_up),CAP_up_vs_bg_p=float(p_up),CAP_down_vs_bg_U=float(u_dn),CAP_down_vs_bg_p=float(p_dn),
                 observed_proximity_up=obs_up,observed_proximity_down=obs_dn,empirical_p_up=float(emp_up),empirical_p_down=float(emp_dn),
                 fisher_top200_CAP_up_overlap=ov_up,fisher_top200_CAP_up_or=float(or_up),fisher_top200_CAP_up_p=float(fp_up),
                 fisher_top200_CAP_down_overlap=ov_dn,fisher_top200_CAP_down_or=float(or_dn),fisher_top200_CAP_down_p=float(fp_dn),
                 corr_absD_spearman=float(rho_abs),corr_absD_p=float(pr_abs),corr_posD_spearman=float(rho_pos),corr_posD_p=float(pr_pos))
    if return_details:
        return metrics, {'df':df,'null_up':null_up,'null_down':null_dn,'cap_up':cap_up,'cap_down':cap_down,'mapped_targets':mapped,'W':W,'idx':idx}
    return metrics
