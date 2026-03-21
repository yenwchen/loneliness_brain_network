#!/bin/bash
#
# ==============================================================================
# run_fmri_network_analysis.slurm
#
# Purpose:
#   SLURM job script to run the fMRI graph-theory network analysis on an
#   HPC cluster. Submits fmri_network_analysis_absedges.py, which computes
#   nodal and global network metrics from subject-level FC matrices using
#   absolute edge weights (dissertation version; main text results).
#
#   To run the positive-only or negative-only sensitivity analyses instead,
#   replace the python script call at the bottom with:
#     python fmri_network_analysis_posedges.py   # positive edges only
#     python fmri_network_analysis_negedges.py   # negative edges only
#   The pipeline and outputs are otherwise identical.
#
# Dependencies:
#   - pandas/1.1.4 (or later)
#   - See fmri_network_analysis_absedges.py for full Python dependencies
#
# Usage:
#   sbatch run_fmri_network_analysis.slurm
# ==============================================================================

#SBATCH --job-name=run_network_py
#SBATCH --time=12:00:00
#SBATCH --nodes=4
#SBATCH --output=run_network_py.txt
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=your.email@institution.edu
#SBATCH -p long-40core

# Note: 4 nodes used due to per-user node limit of 32 on this cluster


# load modules
module load shared
module load pandas/1.1.4


python fmri_network_analysis_absedges.py
