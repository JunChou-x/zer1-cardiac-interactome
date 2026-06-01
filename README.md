# ZER1-IP MS Candidate Prioritization Pipeline

## Overview

This repository contains an R-based workflow for prioritizing candidate ZER1-interacting proteins identified by immunoprecipitation mass spectrometry (IP-MS).

The pipeline integrates:

* UniProt protein sequence retrieval
* N-terminal residue annotation
* Initiator methionine cleavage prediction
* Identification of exposed N-terminal G/A/S/T/C residues
* GO-based functional prioritization
* Candidate ranking
* Publication-ready visualization

Originally developed for mouse cardiac IP-MS datasets, the workflow can be adapted to other protein interaction studies.

---

## Workflow

```text
IP-MS significant proteins
            │
            ▼
Retrieve UniProt sequences
            │
            ▼
Annotate N-terminal residues
            │
            ▼
Predict Met cleavage
            │
            ▼
Identify exposed G/A/S/T/C residues
            │
            ▼
GO functional annotation
            │
            ▼
Tier classification
            │
            ▼
Candidate ranking
            │
            ▼
Visualization
```

---

## Repository Structure

```text
ZER1_pipeline/
│
├── data/
│   ├── ZER1_IPMS_significant_hits.xlsx
│   └── Volcano_input.xlsx
│
├── results/
│
├── 01_sequence_annotation.R
├── 02_GO_prioritization.R
├── 03_visualization.R
│
└── README.md
```

---

## Required Input Files

### 1. IP-MS Candidate List

File:

```text
data/ZER1_IPMS_significant_hits.xlsx
```

Required columns:

| Column      | Description                    |
| ----------- | ------------------------------ |
| Gene Symbol | Official gene symbol           |
| Accession   | UniProt accession ID           |
| log2FC      | Log2 enrichment versus control |
| p-value     | Statistical significance       |

---

### 2. Volcano Plot Input

File:

```text
data/Volcano_input.xlsx
```

Required columns:

| Column      | Description              |
| ----------- | ------------------------ |
| Gene Symbol | Official gene symbol     |
| log2FC      | Log2 fold change         |
| p-value     | Statistical significance |

---

## Pipeline Components

### Step 1: Sequence Annotation

Script:

```r
01_sequence_annotation.R
```

Functions:

* Retrieve protein sequences from UniProt
* Extract first two amino acids
* Predict initiator methionine cleavage
* Identify proteins exposing G, A, S, T, or C after cleavage
* Automatically retry failed UniProt queries

Output:

```text
results/01_sequence_annotation.xlsx
```

---

### Step 2: GO-Based Prioritization

Script:

```r
02_GO_prioritization.R
```

Functions:

* Retrieve Gene Ontology annotations
* Filter proteins using predefined cardiac remodeling-related GO terms
* Assign Tier classifications

GO categories:

* ECM remodeling
* Muscle structure
* Mechanical stress response
* Cellular stress and survival

Outputs:

```text
results/ZER1_Supplementary_Tables.xlsx

results/Tier1_candidates.xlsx
```

---

### Step 3: Visualization

Script:

```r
03_visualization.R
```

Functions:

* Generate Tier 1 ranking lollipop plot
* Generate volcano plot

Outputs:

```text
results/Lollipop_plot.pdf

results/Volcano_plot.pdf
```

---

## Installation

Required R packages:

```r
install.packages(
  c(
    "readxl",
    "writexl",
    "dplyr",
    "ggplot2",
    "ggrepel"
  )
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(
  c(
    "org.Mm.eg.db",
    "GO.db",
    "AnnotationDbi"
  )
)
```

---

## Usage

Run the scripts sequentially:

```r
source("01_sequence_annotation.R")

source("02_GO_prioritization.R")

source("03_visualization.R")
```

---

## Tier Classification Strategy

### Tier 1

Proteins satisfying both criteria:

1. Predicted initiator methionine cleavage

```text
Second residue ∈ {A, C, G, P, S, T, V}
```

2. Functional association with cardiac remodeling-related GO terms

Examples:

* Extracellular matrix organization
* Collagen regulation
* Muscle structure
* Mechanotransduction
* Oxidative stress response

### Tier 2

All remaining proteins.

---

## Outputs

### Supplementary Table S1

All candidate proteins ranked by enrichment.

Columns include:

* Gene symbol
* UniProt accession
* log2FC
* p-value
* N-terminal residues
* Met cleavage prediction
* Tier assignment

---

### Supplementary Table S2

Prioritized Tier 1 candidates.

Columns include:

* Rank
* Gene symbol
* log2FC
* p-value
* Mature N-terminal residue
* GO category
* Matched GO terms

---

## Citation

If you use this workflow in your research, please cite the corresponding publication or repository.

---

## License

This project is distributed under the MIT License.
