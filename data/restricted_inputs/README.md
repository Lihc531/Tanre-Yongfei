# Restricted / external upstream inputs

These files are **not redistributed** in the public repository. They are needed only to rerun upstream stages from the original source material.

## Multi-label GAT

Place the following files in this directory (or edit the paths at the top of `scripts/01_gat/gat_multilabel_model.py`):

- `herb_features_encoded.xlsx` — encoded herb feature matrix used for GAT training
- `syndrome_wind_heat.xlsx` — prescription-derived data for wind-heat invading the lung
- `syndrome_tanre_yongfei.xlsx` — prescription-derived data for Tanre Yongfei
- `syndrome_phlegm_dampness.xlsx` — prescription-derived data for phlegm-dampness obstructing the lung

The final supplementary material contains the feature matrix and held-out performance summaries, but the original CNKI/WanFang prescription workbooks were not present in the recovered `code.zip`. Do not reconstruct them by guessing.

## CAP transcriptomics

The transcriptomic workflow uses public GEO source files for GSE103119 and GSE196399. Download the exact source files from GEO and place them here as:

- `GSE103119_series_matrix.txt.gz`
- `GSE196399_count_matrix.csv.gz`

The included `data/derived/CAP_D_meta.tsv` allows complete downstream network reproduction without downloading GEO data.

## STRING source files

To regenerate the derived PPI edge lists, obtain STRING v12 human files:

- `9606.protein.links.detailed.v12.0.txt.gz`
- `9606.protein.aliases.v12.0.txt.gz`

and run `scripts/03_ppi/build_string_gene_gene_edges.py` at thresholds 400, 700, and 900.

## Herb/target source acquisition

Archived acquisition helpers are in `scripts/optional_acquisition/`. Source websites can change and may impose access/licensing conditions. Verify current terms before rerunning acquisition. The raw CNKI/WanFang prescription corpus and the raw TCMID extraction workbook are not redistributed here.
