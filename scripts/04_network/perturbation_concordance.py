#!/usr/bin/env python
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]

def zscore_sample(x):
    x=np.asarray(x,float); return (x-x.mean())/x.std(ddof=1)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',default=str(ROOT/'results/generated/main_network/all_genes_with_diffusion.tsv')); ap.add_argument('--out',default=str(ROOT/'results/generated/concordance')); ap.add_argument('--n-perm',type=int,default=1000); ap.add_argument('--seed',type=int,default=42)
    a=ap.parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True); df=pd.read_csv(a.input,sep='\t')
    zD=zscore_sample(df.D_meta); zP=zscore_sample(df.P_diffusion); comp=zD*zP; obs=float(comp.sum())
    rng=np.random.RandomState(a.seed); null=np.empty(a.n_perm)
    for i in range(a.n_perm): null[i]=float((zD*zP[rng.permutation(len(zP))]).sum())
    p=float((np.sum(null>=obs)+1)/(a.n_perm+1))
    pd.DataFrame([{'pci_mode':'positive_alignment','observed_pci':obs,'empirical_p':p,'n_perm':a.n_perm,'seed':a.seed}]).to_csv(out/'pci_summary.tsv',sep='\t',index=False)
    pd.DataFrame({'null_pci':null}).to_csv(out/'pci_null.tsv',sep='\t',index=False)
    c=df[['gene_symbol','D_meta','P_diffusion']].copy(); c['comp_score']=comp; c=c.sort_values('comp_score',ascending=False)
    c.to_csv(out/'all_genes_composite_scores.tsv',sep='\t',index=False); c.head(18).to_csv(out/'top18_composite_genes.tsv',sep='\t',index=False)
    print(f'PCI={obs:.6f}; empirical P={p:.6g}')
if __name__=='__main__': main()
