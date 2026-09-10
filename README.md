# Mixed-strain *E. coli* benchmarking dataset and analysis

This repository contains the input data, tool outputs, and analysis code used to evaluate sequence type (MLST), serotype, antimicrobial resistance gene (ARG), and virulence gene (VG) prediction methods on simulated mixed-strain *Escherichia coli* metagenomes.

The analysis script reproduces:

- Summary metadata tables describing simulated datasets
- MLST prediction tables
- MLST performance figures
- Serotype prediction tables
- Serotype performance figures
- Antibiotic resistance gene (ARG) recovery tables
- Virulence gene (VG) recovery tables
- Rare (log-normal) dataset summary tables
- Supplementary LaTeX tables

All figures and tables are generated directly from the provided inputs and tool outputs.

The `figures/` and `tables/` directories contain the manuscript-ready outputs generated from `analysis.R`. 
Running the analysis script will recreate these files from the supplied input data and tool outputs.


## Repository structure

```text
├── analysis.R
├── data/
│   ├── st_serotype_mixed_range_fixed_manifest.tsv
│   └── st_serotype_mixed_range_rare_manifest.tsv
├── results/
│   ├── ecoli_torsten_mlst.tsv
│   ├── compiled_mlst-cge_mlst.tsv
│   ├── compiled_srst2_mlst.tsv
│   ├── assembly_ectyper_serotype.tsv
│   ├── assembly_serotypefinder_serotype.tsv
│   ├── ectyper_combined_serotypes.tsv
│   ├── compiled_srst2_serotype.tsv
│   ├── assembly-amrfinderplus.tsv
│   ├── assembly-rgi.tsv
│   ├── arg-srst2.tsv
│   ├── abricate.tsv
│   ├── virulencefinder.tsv
│   └── virulence-srst2.tsv
├── figures/
├── tables/
└── README.md
```


## Requirements
- R (≥ 4.2 recommended)

Required packages:
- dplyr
- readr
- stringr
- tidyverse
- ggplot2
- knitr
- kableExtra
- rstudioapi
These packages are installed by the R script directly.

## Running the analysis

Open `analysis.R` in RStudio or an R session and run:

```r
source("analysis.R")
```

The script automatically:
1. Imports simulation manifest files from `data/`
2. Imports MLST, serotype, ARG and virulence prediction outputs from `results/`
3. Calculates ARG gene recovery rates (GRR)
4. Calculates virulence gene recovery rates (GRR)
5. Computes average gene recovery rates (AGRR)
6. Generates MLST summary tables
7. Generates serotype summary tables
8. Produces figures
9. Generates manuscript and supplementary LaTeX tables


## Outputs

### Antibiotic resistance gene recovery
ARG recovery tables compare assembly-based and read-based approaches: AMRFinderPlus, RGI, SRST2.
Gene recovery rate (GRR) is calculated as the proportion of expected resistance genes (blaTEM, sul1-3, tetA, tetB, tetD) recovered in a dataset.
Average gene recovery rate (AGRR) is calculated as the mean GRR across all ten simulated datasets.

### Virulence gene recovery
Virulence recovery tables compare: ABRicate (VFDB), VirulenceFinder, SRST2.
For virulence benchmarking, Shiga toxin subtypes are collapsed into toxin-level presence/absence:
- stx1A and stx1B → stx1
- stx2A and stx2B → stx2
Gene recovery rate (GRR) is calculated from the presence of: eae, stx1, stx2.

### MLST
Performance of MLST prediction tools across dominant-strain abundance levels.
```text
figures/mixed_range_mlst_performance.png
```

###  Serotyping
Performance of serotype prediction tools across dominant-strain abundance levels.
```text
figures/mixed_range_serotype_performance.png
```

Performance of serotyping tools evaluated as: O-antigen recovery, H-antigen recovery, Overall antigen recovery, Complete serotype recovery.

The analysis produces LaTeX tables:
- ARG gene benchmarking from assembly tools
- Virulence gene recovery benchmarking from assembly tools
- ARG/VG recovery benchmarking from read-based SRST2
- MLST predictions
- MLST performance summary
- Serotype performance summary
- Overall antigen recovery
- O-antigen recovery
- H-antigen recovery
- Complete serotype recovery


## Reproducibility
All analyses are performed directly from the supplied manifest files and processed tool outputs included in this repository. 
Bash scripts were written to produce the utilised results files from the raw tool outputs.
No external databases, software or internet access are required to reproduce the tables and figures once the required R packages have been installed.


## Citation

If this repository contributes to your work, please cite the associated manuscript.


## Contact

Vaishnavi Divya Shridar

divya.shridar@it.uu.se

Department of Information Technology

Uppsala University
