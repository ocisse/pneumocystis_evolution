Pneumocystis evolution
==============================

Genomic analysis of host specific fungal pathogen Pneumocystis to mammals

[![DOI](https://zenodo.org/badge/263383521.svg)](https://zenodo.org/badge/latestdoi/263383521)

Project Organization
------------

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- Mostly figures
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── src                <- Source code for use in this project.
    │   ├── data           <- Scripts to download or generate data 
        |                  <- Pwk.sm (Pneumocystis wakefieldiae genome assembly notes - snakemake format)
    |   |                  <- Pcan.sh (Pneumocystis canis genome assembly notes)
    |   |                  <- Pmac.sh (Pneumocystis macacae genome assembly notes)
    │   |── envs
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations (R )

Data description (selected):

```
Pipelines:
    src/data/intra_inter_div.smk:   whole genome pairwise divergence computation (snakemake)
    src/data/Pmac.sh:   Pneumocystis macacae wgs assembly (bash)
    src/data/Pwk.smk:   Pneumocystis wakefieldiae assembly (snakemake)
    src/data/Pcan.sh:   Pneumocystis canis assembly (bash)

Datasets:
    docs/inter_intra_div/scores.csv:    whole genome pairwise divergence scores
```
