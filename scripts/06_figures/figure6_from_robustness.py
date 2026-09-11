#!/usr/bin/env python
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--metrics',default=str(ROOT/'results/reference/robustness_ablation_table_s7.tsv')); ap.add_argument('--out',default=str(ROOT/'results/generated/figures/Figure6_reconstructed.png'))
    a=ap.parse_args(); df=pd.read_csv(a.metrics,sep='\t'); out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
    fig,axs=plt.subplots(2,3,figsize=(16,10)); red='#C64E4E'; blue='#4F7DB8'; green='#74B59B'; peach='#E7A39B';
    # A
    ax=axs[0,0]; rows=pd.concat([df[(df.section=='ppi_threshold')],df[(df.section=='target_space') & (df.configuration=='full_syndrome')],df[(df.section=='target_space') & (df.configuration=='core_herbs')]])
    labels=['PPI 400','PPI 700','PPI 900','Broader target space','Core target space']; x=np.arange(5); w=.25
    for off,col,key,lab in [(-w,green,'targets_total','Total targets'),(0,blue,'targets_in_ppi','In PPI'),(w,peach,'targets_mapped_into_cap_ppi','In CAP-PPI')]:
        vals=rows[key].to_numpy(float); ax.bar(x+off,vals,w,label=lab)
        for xx,v in zip(x+off,vals): ax.text(xx,v+5,f'{int(v)}',ha='center',va='bottom',fontsize=8)
    ax.set_xticks(x,labels,rotation=25,ha='right'); ax.set_ylabel('Number of targets'); ax.set_title('A  Target-space coverage',loc='left',fontweight='bold'); ax.legend(frameon=False,fontsize=8)
    # B
    ax=axs[0,1]; r=df[df.section=='raw_vs_diffusion'].copy(); y=np.arange(2); labs=['Raw targets','Network diffusion']
    for i,row in enumerate(r.itertuples(index=False)):
        vals=[row.observed_proximity_down*1e5,row.observed_proximity_up*1e5]; ax.plot(vals,[i,i],color='0.7',lw=1); ax.scatter(vals,[i,i],s=[30,130],c=[blue,red],edgecolor='0.3');
        ax.text(vals[0]+.1,i,f'P={row.empirical_p_down:.3g}',va='center',fontsize=8); ax.text(vals[1]+.1,i,'P<0.001' if row.empirical_p_up<.0015 else f'P={row.empirical_p_up:.3g}',va='center',fontsize=8)
    ax.set_yticks(y,labs); ax.set_xlabel('Observed proximity (×10⁻⁵)'); ax.set_title('B  Raw versus diffusion localization',loc='left',fontweight='bold')
    # C
    ax=axs[0,2]; r=df[df.section=='ppi_threshold'].copy(); xx=r.configuration.astype(int); ax.plot(xx,r.observed_proximity_up*1e5,'o-',c=red,label='CAP-up'); ax.plot(xx,r.observed_proximity_down*1e5,'o-',c=blue,label='CAP-down');
    for x0,y0,p in zip(xx,r.observed_proximity_up*1e5,r.empirical_p_up):
        if p<.05: ax.text(x0,y0+.2,'*',ha='center')
    ax.set_xticks(xx); ax.set_xlabel('STRING confidence threshold'); ax.set_ylabel('Observed proximity (×10⁻⁵)'); ax.set_title('C  PPI-threshold module proximity',loc='left',fontweight='bold'); ax.legend(frameon=False)
    # D
    ax=axs[1,0]; ax.plot(xx,r.fisher_top200_CAP_up_or,'o-',c=red,label='CAP-up'); ax.plot(xx,r.fisher_top200_CAP_down_or,'o-',c=blue,label='CAP-down'); ax.axhline(1,color='0.6',ls='--',lw=1)
    for x0,y0,n,p in zip(xx,r.fisher_top200_CAP_up_or,r.fisher_top200_CAP_up_overlap,r.fisher_top200_CAP_up_p): ax.text(x0,y0+.02,f'n={int(n)}'+('*' if p<.05 else ''),ha='center',fontsize=8)
    for x0,y0,n,p in zip(xx,r.fisher_top200_CAP_down_or,r.fisher_top200_CAP_down_overlap,r.fisher_top200_CAP_down_p): ax.text(x0,y0-.05,f'n={int(n)}'+('*' if p<.05 else ''),ha='center',fontsize=8)
    ax.set_xticks(xx); ax.set_xlabel('STRING confidence threshold'); ax.set_ylabel('Odds ratio'); ax.set_title('D  PPI-threshold top-200 enrichment',loc='left',fontweight='bold'); ax.legend(frameon=False)
    # E
    ax=axs[1,1]; s=df[(df.section=='target_space')].set_index('configuration'); cats=['core_herbs','full_syndrome']; y=np.array([1,0]); labs=['Core target space','Broader target space']; metrics=[('observed_proximity_up',1e5,'CAP-up proximity'),('observed_proximity_down',1e5,'CAP-down proximity'),('fisher_top200_CAP_up_or',1,'Top-200 CAP-up OR'),('corr_posD_spearman',1,'Positive-D Spearman rho')]
    # compact 2x2 inset-style using normalized x positions within one axis
    ax.axis('off'); positions=[(0.05,.58,.38,.34),(0.55,.58,.38,.34),(0.05,.08,.38,.34),(0.55,.08,.38,.34)]
    for (key,scale,title),pos in zip(metrics,positions):
        ia=ax.inset_axes(pos); vals=np.array([s.loc[c,key]*scale for c in cats],float); ia.hlines(y,0,vals,color='0.8'); ia.scatter(vals,y,c=['#E4A35D',blue],s=55,edgecolor='0.3'); ia.set_yticks(y,labs if pos[0]<.1 else []); ia.set_title(title,fontsize=9,fontweight='bold'); ia.tick_params(labelsize=7)
    ax.set_title('E  Target-space definition sensitivity',loc='left',fontweight='bold')
    # F
    ax=axs[1,2]; order=[('target_space','core_herbs','Core target space'),('target_space','full_syndrome','Broader target space'),('ppi_threshold','900','PPI 900'),('ppi_threshold','700','PPI 700'),('ppi_threshold','400','PPI 400'),('raw_vs_diffusion','diffusion_targets','Network diffusion'),('raw_vs_diffusion','raw_targets','Raw targets')]
    cols=[('CAP_up_vs_bg_p','CAP-up\nMWU P'),('CAP_down_vs_bg_p','CAP-down\nMWU P'),('observed_proximity_up','CAP-up\nproximity'),('observed_proximity_down','CAP-down\nproximity'),('fisher_top200_CAP_up_or','Top-200\nCAP-up OR'),('corr_posD_spearman','Positive-D\nSpearman rho')]
    mat=[]; labels=[]
    for sec,conf,lab in order:
        rr=df[(df.section==sec)&(df.configuration.astype(str)==str(conf))].iloc[0]; mat.append([float(rr[k]) for k,_ in cols]); labels.append(lab)
    mat=np.array(mat); scaled=np.zeros_like(mat)
    for j in range(mat.shape[1]):
        lo,hi=np.nanmin(mat[:,j]),np.nanmax(mat[:,j]); scaled[:,j]=.5 if hi==lo else (mat[:,j]-lo)/(hi-lo)
    im=ax.imshow(scaled,aspect='auto',cmap='OrRd',vmin=0,vmax=1); ax.set_yticks(range(len(labels)),labels); ax.set_xticks(range(len(cols)),[x[1] for x in cols],rotation=55,ha='right',fontsize=7); ax.set_title('F  Multi-metric stability map',loc='left',fontweight='bold')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v=mat[i,j]; txt='<0.001' if j<2 and v<.001 else (f'{v:.2e}' if j in (2,3) else f'{v:.3g}'); ax.text(j,i,txt,ha='center',va='center',fontsize=6)
    for i,(sec,conf,_) in enumerate(order):
        if conf in {'core_herbs','700','diffusion_targets'}:
            ax.add_patch(plt.Rectangle((-0.5,i-0.5),len(cols),1,fill=False,edgecolor=red,lw=1.5))
    fig.colorbar(im,ax=ax,fraction=.046,pad=.04,label='Within-metric relative signal')
    fig.tight_layout(); fig.savefig(out,dpi=300,bbox_inches='tight'); print('Saved:',out)
if __name__=='__main__': main()
