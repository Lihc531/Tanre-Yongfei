#!/usr/bin/env python
from pathlib import Path
import argparse, json, subprocess, sys, tempfile
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]

def run_one(label, ppi, targets, mode, n_perm, seed, tmpdir, null_up=None, null_down=None):
    out=tmpdir/f'{label}.json'
    cmd=[sys.executable,str(ROOT/'scripts/05_robustness/run_single_configuration.py'),
         '--ppi',str(ppi),'--targets',str(targets),'--mode',mode,
         '--alpha','0.5','--n-perm',str(n_perm),'--seed',str(seed),'--out',str(out)]
    if null_up is not None and null_down is not None:
        cmd += ['--null-up',str(null_up),'--null-down',str(null_down)]
    print('Running',label,'...',flush=True)
    subprocess.run(cmd,cwd=ROOT,check=True)
    return json.loads(out.read_text(encoding='utf-8'))

def main():
    ap=argparse.ArgumentParser(description='Reproduce the final robustness/ablation configurations (alpha=0.5).')
    ap.add_argument('--n-perm',type=int,default=1000); ap.add_argument('--seed',type=int,default=42)
    ap.add_argument('--out',default=str(ROOT/'results/generated/robustness'))
    ap.add_argument('--main-null-dir',default=str(ROOT/'results/generated/main_network'),help='Reuse the already generated 700/core null distributions when available.')
    a=ap.parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    cap=ROOT/'data/derived/CAP_D_meta.tsv'; core=ROOT/'data/derived/target_genes_core.txt'; full=ROOT/'data/derived/target_genes_full.txt'
    ppi={t:ROOT/f'data/derived/string_ppi_edges_{t}.tsv.gz' for t in ('400','700','900')}
    null_dir=Path(a.main_null_dir); up_file=null_dir/'null_up.tsv'; dn_file=null_dir/'null_down.tsv'
    reuse = up_file.exists() and dn_file.exists() and len(pd.read_csv(up_file,sep='\t'))==a.n_perm and len(pd.read_csv(dn_file,sep='\t'))==a.n_perm
    nu=up_file if reuse else None; nd=dn_file if reuse else None
    with tempfile.TemporaryDirectory(prefix='robustness_',dir=out) as td:
        td=Path(td)
        # The 400-threshold network is the largest configuration, so run it first.
        # This avoids late-run memory/CPU pressure on constrained review environments.
        r400=run_one('ppi400',ppi['400'],core,'diffusion',a.n_perm,a.seed,td)
        r900=run_one('ppi900',ppi['900'],core,'diffusion',a.n_perm,a.seed,td)
        full700=run_one('full700',ppi['700'],full,'diffusion',a.n_perm,a.seed,td)
        diff700=run_one('diff700',ppi['700'],core,'diffusion',a.n_perm,a.seed,td,nu,nd)
        raw700=run_one('raw700',ppi['700'],core,'raw',a.n_perm,a.seed,td,nu,nd)
    rows=[
        {'section':'raw_vs_diffusion','configuration':'raw_targets',**raw700},
        {'section':'raw_vs_diffusion','configuration':'diffusion_targets',**diff700},
        {'section':'ppi_threshold','configuration':'400',**r400},
        {'section':'ppi_threshold','configuration':'700',**diff700},
        {'section':'ppi_threshold','configuration':'900',**r900},
        {'section':'target_space','configuration':'full_syndrome',**full700},
        {'section':'target_space','configuration':'core_herbs',**diff700},
    ]
    df=pd.DataFrame(rows); df.to_csv(out/'robustness_all_configurations.tsv',sep='\t',index=False)
    print(df.to_string(index=False))
if __name__=='__main__': main()
