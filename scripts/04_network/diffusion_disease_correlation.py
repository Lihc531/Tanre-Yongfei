#!/usr/bin/env python
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
ROOT=Path(__file__).resolve().parents[2]

def safe_corr(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); m=np.isfinite(x)&np.isfinite(y); x=x[m];y=y[m]
    if len(x)<3: return len(x),np.nan,np.nan,np.nan,np.nan
    sr,sp=spearmanr(x,y); pr,pp=pearsonr(x,y); return len(x),sr,sp,pr,pp

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',default=str(ROOT/'results/generated/main_network/all_genes_with_diffusion.tsv')); ap.add_argument('--out',default=str(ROOT/'results/generated/correlation'))
    a=ap.parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True); df=pd.read_csv(a.input,sep='\t')
    df['abs_D_meta']=df.D_meta.abs(); df['positive_D']=df.D_meta.clip(lower=0); df['negative_D_abs']=(-df.D_meta).clip(lower=0)
    rows=[]
    specs=[('all_genes',df,[('D_meta','D_meta'),('abs_D_meta','abs_D_meta'),('positive_D','positive_D'),('negative_D_abs','negative_D_abs')]),
           ('CAP_up',df[df.group=='CAP_up'],[('D_meta','D_meta'),('abs_D_meta','abs_D_meta'),('positive_D','positive_D')]),
           ('CAP_down',df[df.group=='CAP_down'],[('D_meta','D_meta'),('abs_D_meta','abs_D_meta'),('negative_D_abs','negative_D_abs')])]
    for subset,d,ys in specs:
        for yname,col in ys:
            n,sr,sp,pr,pp=safe_corr(d.P_diffusion,d[col]); rows.append(dict(subset=subset,x='P_diffusion',y=yname,n=n,spearman_rho=sr,spearman_p=sp,pearson_r=pr,pearson_p=pp))
    res=pd.DataFrame(rows); res.to_csv(out/'correlation_summary.tsv',sep='\t',index=False); print(res.to_string(index=False))
if __name__=='__main__': main()
