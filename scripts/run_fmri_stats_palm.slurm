#!/bin/bash
#
# ==============================================================================
# run_fmri_stats_palm.slurm
#
# Purpose:
#   SLURM job script to run FSL PALM permutation tests on graph-theory network
#   metrics derived from absolute edge weights (dissertation version; N = 512).
#   Executes fmri_stats_palm.m via MATLAB in no-display mode.
#
# Dependencies:
#   - MATLAB 2021a
#   - FSL PALM (must be accessible from within MATLAB)
#   - fmri_stats_palm.m and all input CSVs in the working directory
#
# Usage:
#   sbatch run_fmri_stats_palm.slurm
# ==============================================================================

#SBATCH --job-name=run_palm
#SBATCH --time=04:00:00
#SBATCH --ntasks-per-node=40
#SBATCH --output=out_run_palm.txt
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=your.email@institution.edu
#SBATCH -p short-40core


# ==============================================================================
# USER CONFIGURATION — update this path before running
# ==============================================================================

# Directory containing fmri_stats_palm.m and all PALM input CSVs
palm_dir="/path/to/HCA_stat/palm_models/abs_model"

# ==============================================================================


# load modules
module load matlab/2021a


# Create PALM output subdirectories (one per network metric)
mkdir -p "$palm_dir"/palm_output/betweencent
mkdir -p "$palm_dir"/palm_output/closecent
mkdir -p "$palm_dir"/palm_output/strength
mkdir -p "$palm_dir"/palm_output/normstrength
mkdir -p "$palm_dir"/palm_output/eigencent
mkdir -p "$palm_dir"/palm_output/pgcent
mkdir -p "$palm_dir"/palm_output/clustercoef
mkdir -p "$palm_dir"/palm_output/partcoef
mkdir -p "$palm_dir"/palm_output/modularity
mkdir -p "$palm_dir"/palm_output/average_shortest_path


# Run PALM via MATLAB in no-display mode from the palm_dir
# -nodisplay: suppress GUI; -r: run script without .m extension
cd "$palm_dir"
matlab -nodisplay -r fmri_stats_palm
