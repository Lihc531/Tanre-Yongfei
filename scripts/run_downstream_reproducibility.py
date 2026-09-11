#!/usr/bin/env python
from pathlib import Path
import argparse, subprocess, sys
ROOT=Path(__file__).resolve().parents[1]

def run(*args):
    print('\n$', ' '.join(map(str,args)))
    subprocess.run(list(map(str,args)),cwd=ROOT,check=True)

def main():
    ap=argparse.ArgumentParser(description='Reproduce downstream manuscript analyses from included derived inputs.'); ap.add_argument('--smoke',action='store_true',help='Use 5 permutations and skip full robustness verification.'); a=ap.parse_args()
    n=5 if a.smoke else 1000
    py=sys.executable
    run(py,ROOT/'scripts/04_network/run_main_network_analysis.py','--n-perm',n)
    run(py,ROOT/'scripts/04_network/diffusion_disease_correlation.py')
    run(py,ROOT/'scripts/04_network/perturbation_concordance.py','--n-perm',n)
    if not a.smoke:
        run(py,ROOT/'scripts/05_robustness/run_robustness.py','--n-perm',n)
    run(py,ROOT/'scripts/verification/verify_downstream.py')
    run(py,ROOT/'scripts/06_figures/figure6_from_robustness.py','--metrics',ROOT/'results/reference/robustness_ablation_table_s7.tsv')
if __name__=='__main__': main()
