# CellSparse: A Statistical Algorithm for Detecting Upregulated Genes in scRNA-seq Data

**Author:** Nina Gilshteyn, M.S.  
**Affiliation:** Deeds Lab, UCLA — Department of Integrative Biology & Physiology  
**Status:** Manuscript in preparation — *"The Biology is in The Tails: Skewness and Upregulation in scRNA-seq Data"* which is based off of my Master's thesis: Gilshteyn, N. (2024). The Biology is in The Tails: Skewness and Upregulation in scRNA-seq. UCLA. ProQuest ID: Gilshteyn_ucla_0031N_23577. Merritt ID: ark:/13030/m5b96p3g. Retrieved from https://escholarship.org/uc/item/8p46s4ck  
**Profile:** [ninagilshteyn.github.io](https://ninagilshteyn.github.io) · [github.com/NinaGilshteyn](https://github.com/NinaGilshteyn)

---

## Overview

CellSparse is a parallelized statistical algorithm for identifying upregulated genes in single-cell RNA sequencing (scRNA-seq) data. Rather than relying on conventional differential expression methods, CellSparse uses a **sparse null simulation approach** grounded in the distributional properties of scRNA-seq count data — specifically, the skewness between the mean and median of gene expression distributions.

The CellSparse algorithm tests the null hypothesis that outliers cells in a gene's histogram are caused by cells that are generally more transcriptionally active or 'heavy'. This algorithm simulated this question by randomly redistributing each cell's total counts to it's endogenously expressed genes. This test is appropriate for the question because heavier cells would have a higher chance of increasing the amount of counts for a gene; thus, the histogram for said gene would become right skewed by random chance. We ran 20,000 simulations using the Hoffman2 Cluster and check how many genes have significantly large weights relative to this null hypothesis. 

The pipeline is designed to run on an HPC cluster and leverages Python's `multiprocessing` library to parallelize null simulations across CPUs — enabling large-scale in-silico gene simulations that would otherwise be RAM-constrained.

---

## Pipeline overview

```
Input: scRNA-seq expression matrix (.csv)
         │
         ▼
Step 1 — Identify number and location of expressed genes per cell
         │
         ▼
Step 2 — Calculate total transcript count per cell
         │
         ▼
Step 3 — Run parallelized sparse size null simulation (HPC)
         │
         ▼
Step 4a — Process and flatten null simulation output files
Step 4b — Calculate observed mean–median weights from real data
Step 4c — Compute empirical two-sided p-values
         │
         ▼
Step 5 — Multiple testing correction (FDR/BH) + export results
         │
         ▼
Output: significant_genes.csv, corrected p-values, histograms
```

---

## File descriptions

### `Step1_NumberandLocation_of_Expressed_Genes.py`
Reads the input scRNA-seq expression matrix, transposes it so that genes are rows and cells are columns, and identifies for each cell which gene indices have non-zero expression. Outputs:
- `AmountExpressedGenes.csv` — number of expressed genes per cell
- `Coordinates_of_Expressed_Genes.csv` — row indices of expressed genes per cell

### `Step2_Number_of_Transcripts_in_eachcell.py`
Calculates the total RNA transcript count per cell by summing across all genes. Outputs:
- `Total_RNA_per_Cell.csv`

### `Step3_Parallelized_Sparse_Size_Null_Sim.py`
The core simulation engine. For each iteration of the null simulation, distributes cells across CPUs using `multiprocessing.Pool`, runs the sparse size null model, collects results, assembles them into a gene-indexed linked list structure, and computes summary statistics (mean, median) per gene per simulation. Writes batches of null statistics to disk at defined checkpoints. Outputs:
- `jobnumber_{n}_simulationnumber_{n}_null_stats.csv` (one per checkpoint)

Depends on helper modules: `SparseSizeNull`, `statscalc`, `jag_arr_to_linked_list`, `compressedwriter`.

### `Step4a_process_files.py`
Reads the combined null simulation output (`total_combined.csv`), flattens nested list structures, and computes the mean–median weight for each gene across all null simulations. Outputs:
- `correct_null_weights.csv`

### `Step4b_calc_observed_values.py`
Computes the observed mean–median weight for each expressed gene in the real data. Filters out genes with zero total expression before calculating weights. Outputs:
- `obs_weight.csv`

### `Step4c_p-val_calc.py`
Computes empirical two-sided p-values by comparing each gene's observed weight against its null weight distribution. For each gene: counts how many null values are ≥ the observed weight (upper tail) and how many are ≤ its negative (lower tail), then divides by twice the total number of simulations. Outputs:
- `newnewpvals_for_gene_index_as_column_index.csv`

### `Step5_correct_pvals.py`
Applies Benjamini–Hochberg FDR correction to all p-values using `statsmodels`. Generates histograms of the raw and corrected p-value distributions, then exports results sorted by significance and observed weight. Outputs:
- `fulldf_corrected_pval.csv`
- `significant_genes.csv`
- `not_significant_genes.csv`
- `histogram_of_nonzerop_vals.png`
- `corrected_pvals.png`

---

## Requirements

```
Python 3.8+
pandas
numpy
seaborn
statsmodels
multiprocessing (standard library)
```

Install dependencies:
```bash
pip install pandas numpy seaborn statsmodels
```

HPC environment (UCLA Hoffman2 or compatible SLURM cluster) is recommended for Step 3 due to the computational demands of the null simulation.

---

## Usage

### 1. Prepare your input
Provide a `.csv` expression matrix where rows are cells and columns are genes. A `labels` column will be dropped automatically.

### 2. Run Steps 1 and 2
```bash
python Step1_NumberandLocation_of_Expressed_Genes.py
python Step2_Number_of_Transcripts_in_eachcell.py
```

### 3. Run the null simulation (HPC)
Submit `Step3_Parallelized_Sparse_Size_Null_Sim.py` via your cluster job scheduler, passing the simulation batch number as an argument:
```bash
python Step3_Parallelized_Sparse_Size_Null_Sim.py <job_number>
```
Combine all output CSVs into `total_combined.csv` before proceeding. Linux script recommended. 

### 4. Compute weights and p-values
```bash
python Step4a_process_files.py
python Step4b_calc_observed_values.py
python Step4c_p-val_calc.py
```

### 5. Correct p-values and export results
```bash
python Step5_correct_pvals.py
```

---

## Output files

| File | Description |
|---|---|
| `significant_genes.csv` | Genes passing FDR correction, sorted by observed weight |
| `not_significant_genes.csv` | Non-significant genes |
| `fulldf_corrected_pval.csv` | All genes with corrected p-values and observed weights |
| `corrected_pvals.png` | Histogram of BH-corrected p-value distribution |
| `histogram_of_nonzerop_vals.png` | Histogram of non-zero raw p-values |

---

## Citation

If you use CellSparse in your research, please cite:

> Gilshteyn, N. & Deeds, E.J. (in preparation). *The Biology is in The Tails: Skewness and Upregulation in scRNA-seq Data.*

---

## License

This software is released under the **CellSparse Commercial Use License**. Free for academic and non-commercial research use. Commercial use requires a separate written license agreement. See `LICENSE` for full terms.
