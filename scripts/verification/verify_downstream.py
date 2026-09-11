#!/usr/bin/env python
from pathlib import Path
import argparse, math
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]

def close(a,b,rtol=2e-6,atol=1e-12): return math.isclose(float(a),float(b),rel_tol=rtol,abs_tol=atol)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--generated',default=str(ROOT/'results/generated')); a=ap.parse_args(); g=Path(a.generated)
    errs=[]
    # main summary
    fs=pd.read_csv(g/'main_network/final_summary.tsv',sep='\t').set_index('metric').value
    expected={'ppi_edges':236837,'ppi_nodes':16194,'cap_genes_in_ppi':10245,'targets_total':337,'targets_in_ppi':332,'targets_mapped_into_cap_ppi':203,'cap_up_genes':2796,'cap_down_genes':2398,'observed_proximity_up':5.27224330582954e-05,'observed_proximity_down':7.37877154350432e-06}
    aliases={k:k for k in expected}
    for k,v in expected.items():
        got=fs[aliases[k]]
        if not close(got,v): errs.append(f'{k}: got {got}, expected {v}')
    fish=pd.read_csv(g/'main_network/fisher_top200.tsv',sep='\t'); up=fish[fish.signature=='CAP_up'].iloc[0]
    for k,v in [('a_top_and_sig',73),('oddsratio',1.54561463878711),('pvalue',0.00254183083187808)]:
        if not close(up[k],v): errs.append(f'top200 {k}: got {up[k]}, expected {v}')
    corr=pd.read_csv(g/'correlation/correlation_summary.tsv',sep='\t'); row=corr[(corr.subset=='all_genes')&(corr.y=='positive_D')].iloc[0]
    if not close(row.spearman_rho,0.149238974678315): errs.append(f'positive-D rho: {row.spearman_rho}')
    pci=pd.read_csv(g/'concordance/pci_summary.tsv',sep='\t').iloc[0]
    if not close(pci.observed_pci,571.5230339296656,rtol=2e-6): errs.append(f'PCI: {pci.observed_pci}')
    # robustness compare to final Table S7
    if (g/'robustness/robustness_all_configurations.tsv').exists():
        gen=pd.read_csv(g/'robustness/robustness_all_configurations.tsv',sep='\t'); ref=pd.read_csv(ROOT/'results/reference/robustness_ablation_table_s7.tsv',sep='\t')
        def norm_config(x):
            s=str(x).strip()
            try:
                f=float(s)
                if f.is_integer(): return str(int(f))
            except ValueError:
                pass
            return s
        gen['_configuration_norm']=gen.configuration.map(norm_config)
        ref['_configuration_norm']=ref.configuration.map(norm_config)
        for _,rr in ref.iterrows():
            m=gen[(gen.section==rr.section)&(gen._configuration_norm==rr._configuration_norm)]
            if len(m)!=1: errs.append(f'missing robustness row {rr.section}/{rr.configuration}'); continue
            gr=m.iloc[0]
            for col in ['ppi_edges','ppi_nodes','cap_genes_in_ppi','targets_total','targets_in_ppi','targets_mapped_into_cap_ppi','observed_proximity_up','observed_proximity_down','fisher_top200_CAP_up_or','corr_posD_spearman']:
                # Rank correlations can differ in the final decimal places across SciPy/R versions
                # because of tie handling; the tolerance remains far below manuscript precision.
                rtol = 1e-5 if col == 'corr_posD_spearman' else 3e-6
                if not close(gr[col],rr[col],rtol=rtol): errs.append(f'robustness {rr.section}/{rr.configuration} {col}: {gr[col]} vs {rr[col]}')
    if errs:
        print('VERIFICATION FAILED'); print('\n'.join('- '+e for e in errs)); raise SystemExit(1)
    print('VERIFICATION PASSED: generated downstream metrics agree with final reference values within tolerance.')
if __name__=='__main__': main()
