# Data availability for code reproduction

## Included derived inputs

The repository includes the final CAP meta-disease vector, final target-gene lists, and derived STRING gene-gene edge lists used for downstream analysis. These are sufficient to rerun the main network diffusion, correlation, PCI, and robustness analyses.

## Public source data not duplicated as raw files

The CAP transcriptomic source data are available from GEO under accession numbers **GSE103119** and **GSE196399**. The STRING PPI source files can be downloaded from STRING; the repository includes the derived thresholded gene-gene edge lists to make the downstream analysis reproducible without a large raw-data download.

## Inputs not available in the recovered code archive

The original three syndrome-specific prescription workbooks used for GAT graph construction were not present in the recovered code archive. The final supplementary material reports the herb feature matrix and held-out model performance, but those missing prescription-level source workbooks cannot be recreated faithfully from the available files alone.

## Third-party source material

The raw CNKI/WanFang prescription corpus and raw web-extracted TCMID records are not redistributed in this repository. The repository instead includes final derived target-gene sets for the downstream analysis and documents the acquisition/processing code for traceability.

Users should comply with the licensing and terms of use of all third-party data sources.
