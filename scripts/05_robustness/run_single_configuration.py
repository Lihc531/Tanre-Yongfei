#!/usr/bin/env python
from pathlib import Path
import argparse, json, sys
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/04_network'))
from network_core import analyze_target_space

def load_null(path):
    if not path: return None
    return pd.read_csv(path,sep='\t').iloc[:,0].to_numpy(float)

def main():
    ap=argparse.ArgumentParser(description='Run one network-analysis configuration in an isolated process.')
    ap.add_argument('--ppi',required=True); ap.add_argument('--targets',required=True)
    ap.add_argument('--mode',choices=['raw','diffusion'],default='diffusion')
    ap.add_argument('--cap',default=str(ROOT/'data/derived/CAP_D_meta.tsv'))
    ap.add_argument('--alpha',type=float,default=0.5); ap.add_argument('--n-perm',type=int,default=1000); ap.add_argument('--seed',type=int,default=42)
    ap.add_argument('--null-up'); ap.add_argument('--null-down'); ap.add_argument('--out',required=True)
    a=ap.parse_args()
    nu=load_null(a.null_up); nd=load_null(a.null_down)
    r=analyze_target_space(a.cap,a.ppi,a.targets,mode=a.mode,alpha=a.alpha,n_perm=a.n_perm,random_seed=a.seed,
                           null_up_override=nu,null_down_override=nd)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    Path(a.out).write_text(json.dumps(r,indent=2,sort_keys=True),encoding='utf-8')
if __name__=='__main__': main()
