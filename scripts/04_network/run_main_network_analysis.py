#!/usr/bin/env python
from pathlib import Path
import argparse
import pandas as pd
from scipy.stats import mannwhitneyu
from network_core import analyze_target_space, fisher_enrichment

ROOT=Path(__file__).resolve().parents[2]

def main():
    ap=argparse.ArgumentParser(description='Run final alpha=0.5 RWR/module-localization analysis.')
    ap.add_argument('--cap', default=str(ROOT/'data/derived/CAP_D_meta.tsv'))
    ap.add_argument('--ppi', default=str(ROOT/'data/derived/string_ppi_edges_700.tsv.gz'))
    ap.add_argument('--targets', default=str(ROOT/'data/derived/target_genes_core.txt'))
    ap.add_argument('--out', default=str(ROOT/'results/generated/main_network'))
    ap.add_argument('--alpha', type=float, default=0.5)
    ap.add_argument('--n-perm', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=42)
    args=ap.parse_args()
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    metrics,det=analyze_target_space(args.cap,args.ppi,args.targets,mode='diffusion',alpha=args.alpha,n_perm=args.n_perm,random_seed=args.seed,return_details=True)
    df=det['df'].rename(columns={'score':'P_diffusion'}).copy()
    df['group']='background'; df.loc[df.gene_symbol.isin(det['cap_up']),'group']='CAP_up'; df.loc[df.gene_symbol.isin(det['cap_down']),'group']='CAP_down'
    df=df.sort_values('P_diffusion',ascending=False).reset_index(drop=True)
    df.to_csv(out/'all_genes_with_diffusion.tsv',sep='\t',index=False)
    pd.DataFrame({'null_up':det['null_up']}).to_csv(out/'null_up.tsv',sep='\t',index=False)
    pd.DataFrame({'null_down':det['null_down']}).to_csv(out/'null_down.tsv',sep='\t',index=False)
    pd.DataFrame([{'metric':k,'value':v} for k,v in metrics.items() if k!='mode']).to_csv(out/'final_summary.tsv',sep='\t',index=False)
    for n in (200,500):
        top=df.head(n); top_set=set(top.gene_symbol); up=det['cap_up']; dn=det['cap_down']; genes=set(df.gene_symbol)
        rows=[]
        for label,sig in [('CAP_up',up),('CAP_down',dn)]:
            orr,p,ov=fisher_enrichment(genes,top_set,sig)
            a=ov; b=len(top_set-sig); c=len(sig-top_set); d=len(genes-top_set-sig)
            rows.append(dict(top_n=n,signature=label,a_top_and_sig=a,b_top_not_sig=b,c_sig_not_top=c,d_other=d,oddsratio=orr,pvalue=p))
        pd.DataFrame(rows).to_csv(out/f'fisher_top{n}.tsv',sep='\t',index=False)
        top.to_csv(out/f'top{n}_diffusion_genes.tsv',sep='\t',index=False)
    print(pd.DataFrame([metrics]).to_string(index=False))
    print('Saved:',out)
if __name__=='__main__': main()
