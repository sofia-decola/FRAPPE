# FRAPPE
### Framework for Rhythmic Analysis of Phasic Proteomics Experiments
 
FRAPPE is a Python-based computational framework for analyzing high-dimensional proteomics time-series datasets. It was designed to help biologists generate hypotheses about how experimental conditions impact the stability and rhythmicity of the proteome, particularly in the context of circadian biology.
 
---
 
## Overview
 
Proteomics time-series data captures dynamic protein expression over time, but most existing tools are built for cross-sectional analysis. FRAPPE fills this gap by providing a suite of methods to detect phase shifts, period changes, and co-expression patterns across experimental groups.
 
 
Each module contains tutorials for use as well as the code itself for accessing the methods
---
 
## Modules 
 
### 1. Time-Lagged Cross Correlation (TLCC)
Compares two time-series proteomics datasets by circularly shifting one series across time and computing correlation at each lag. This identifies phase delays between a control and experimental group at the individual protein and proteome-wide level.
 
- Identifies the lag at which peak correlation occurs for each protein
- Can plot a histogram of phase delay distribution across the proteome
### 2. Rhythmic Period Analysis (RPA)
Evaluates rhythmic stability of protein expression across experimental groups by analyzing changes in oscillatory period length.
 
- Identifies proteins that gain, lose, or maintain rhythmicity between conditions 
- Calculates the magnitude of period shifts per protein
- Default visualization plots the distribution of period changes across the proteome, with an Anderson-Darling test for normality reported alongside the mean shift
### 3. Adjacency Matrix & Percent Correlation Threshold Analysis
Constructs a correlation matrix of protein expression values at time point 0 and applies a user-defined threshold to assess co-expression maintenance across experimental groups.
 
- Default correlation threshold: **0.8**
- Categorizes proteins into four groups:
  - `0 → 0`: Never above threshold
  - `1 → 1`: Maintains co-expression across all conditions
  - `0 → 1`: Gains co-expression
  - `1 → 0`: Loses co-expression
- Calculates percent of the proteome in each category for large-scale comparisons
### 4. Community Detection & StringDB Analysis
Clusters proteins by similar temporal expression patterns using network-based community detection, then queries the STRING protein-protein interaction (PPI) database to identify functional relationships within clusters.
 
- Uses the **Leiden algorithm** via `leidenalg` for community detection
- StringDB queries support **custom model organism** (scientific name or NCBI Taxon ID) and **custom functional annotation databases** (e.g., KEGG, Gene Ontology)
- Can be applied to individual proteins or entire expression clusters
---
**Key dependencies:**
- `pandas`, `numpy`, `scipy`
- `matplotlib`
- `scikit-learn`
- `python-igraph`, `leidenalg`
- `networkx`
- `requests`
- `kneed`
- `alive-progress`
- Python 3.10+
---
 
## Important -- Input Data Formatting
 
FRAPPE accepts any preprocessed proteomics time-series data. The following preprocessing pipeline is preferred: 
 
**ECHO** (Extended Circadian Harmonic Oscillator) for curve fitting, period estimation, and oscillation characterization
    Protein name (sometimes the column is 'Gene Name'), fitted timepoints and period values are necessary from ECHO outputs
 
Data should be sampled at regular intervals (e.g., every 3 hours) across a full circadian cycle or longer.
If data cannot be preprocessed into rhythms, the modules will work on any data in the following format
 
---
## Data Format Requirements - TLCC
Data should be a .csv (comma-separated values) with the first column being the protein names/labels, and the rest being numeric columns with expression data, ordered by time point. 
 
Certain methods allow for 'protein_col' parameter, which allow the user to specify the name of this column (Default is ‘Gene Name’, as in ECHO preprocessed tables) 
 
An example of this formatting is the following:
 
| Protein |	 Timepoint 1| 	Timepoint 2	| Timepoint 3| 	Timepoint 4| 	Timepoint 5| Timepoint 6 |
| ------------- |-------------|-------------|-------------|-------------|-------------|-------------|
| Sample 1 |	1.633117905| 1.4448585	| 	1.513810213| 	1.309553546 | 	1.302488129| 1.8384729 	|	
| Sample 2 | 	-0.630319173| -0.4228347	| 	-0.510500938| -0.5228194 	| 	-0.543457041| 	-0.448383157|		
| Sample 3	| -0.780221402| 0.17483084	| 0.238884	| 	0.178429468| 	0.306513019| 	1.376226634|
 
 
## Usage → Every module has a Tutorial .ipynb in src/tutorials folder
 
 
## License
 
This project is licensed under the MIT License. See `LICENSE` for details.