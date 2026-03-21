# Social Disconnection in the Brain: Loneliness and Age across Networks using Graph Theory

**Author:** Yen-Wen Chen

**Preprint:** https://doi.org/10.64898/2026.02.03.703621 | Under review at *Social Cognitive and Affective Neuroscience (SCAN)*


## Overview

This project examines the relationship between loneliness, age, and resting-state brain network organization using the HCP-Aging dataset and graph theory. While loneliness showed no independent main effects on neural measures, a significant loneliness × age interaction emerged for local connectivity in the Default Mode and Frontoparietal networks. Older age was independently associated with lower functional connectivity and less modular brain organization.


## Repository Structure

```
scripts/
├── Behavioral
│   ├── behav_preproc_recode_variables.Rmd        # Load and recode HCP-A behavioral/demographic variables
│   └── behav_stats_loneliness.Rmd                # Descriptive stats and regression models for loneliness
│
├── fMRI — Preprocessing & QC
│   ├── fmri_qc_motion_fd.Rmd                     # QC: framewise displacement per subject/run
│   └── fmri_preproc_demean_parcellate_cormat.slurm  # Demean → parcellate → cormat → average
│
├── fMRI — Network Analysis
│   ├── fmri_network_analysis_absedges.py          # Graph metrics: absolute edges (main text)
│   ├── fmri_network_analysis_posedges.py          # Graph metrics: positive edges only (supplementary)
│   └── fmri_network_analysis_negedges.py          # Graph metrics: negative edges only (supplementary)
│
├── fMRI — Statistics
│   ├── fmri_preproc_compile_nodal_palm.Rmd        # Compile nodal data and write PALM input files (absolute edges)
│   ├── fmri_stats_palm.m                          # FSL PALM permutation tests (N = 512)
│   ├── fmri_stats_palm_summarize.Rmd              # Summarize PALM output; scatter and interaction plots (absolute edges)
│   ├── run_fmri_network_analysis.slurm            # SLURM runner for network analysis
│   └── run_fmri_stats_palm.slurm                  # SLURM runner for PALM
```


## Pipeline Overview

```
Behavioral:
  behav_preproc_recode_variables.Rmd → datainfo.csv
  behav_stats_loneliness.Rmd         → regression models, descriptive stats

fMRI:
  Raw HCA resting-state fMRI
    → fmri_preproc_demean_parcellate_cormat.slurm
        Demean → Parcellate (Cole-Anticevic atlas) → Correlation matrix → Average across runs

    → fmri_network_analysis_{absedges|posedges|negedges}.py
        Build weighted graph → Compute nodal and global network metrics

    → fmri_preproc_compile_nodal_palm.Rmd
        Merge behavioral/demographic covariates → Write PALM input CSVs

    → fmri_stats_palm.m  (via run_fmri_stats_palm.slurm)
        FSL PALM permutation testing (N = 512, n = 5000 permutations)

    → fmri_stats_palm_summarize.Rmd
        Summarize PALM outputs → Transpose for Workbench → Scatter and interaction plots

QC (run separately):
    → fmri_qc_motion_fd.Rmd
        Framewise displacement timeseries and motion exclusion flags
```


## Post-processing: Visualizing Nodal Results in wb_view

PALM outputs (t-stats, FWE-p) can be converted back to CIFTI format for
visualization in Connectome Workbench (`wb_view`):

```bash
wb_command -cifti-convert -from-text -reset-scalars <text-in> <cifti-template> <cifti-out>
```

Steps:
1. Transpose the PALM output so that rows = parcels, columns = subjects (`fmri_stats_palm_summarize.Rmd`)
2. Use `cifti_label_import_from_CAplabel.plabel.nii` as the CIFTI template
3. Run `wb_command -cifti-convert -from-text` to convert all subjects in one go
4. Open the resulting CIFTI file in `wb_view`


## Data

This study used the **HCP-Aging Lifespan 2.0 Release** from the [Lifespan Human Connectome Project Aging (HCP-Aging)](https://www.humanconnectome.org/study/hcp-lifespan-aging). Data are available through the [NIMH Data Archive](https://nda.nih.gov/).


## Requirements

**Python** (network analysis)
- `networkx` 2.4
- `community` (python-louvain)
- `pandas`, `numpy`

**R** v4.0.3 (behavioral analysis, QC, PALM post-processing)
- `tidyverse`, `ggplot2`, `gridExtra`
- `broom`, `rstatix`

**MATLAB** (permutation statistics)
- [FSL PALM](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM)

**Connectome Workbench** v1.5.0 (preprocessing and visualization)
