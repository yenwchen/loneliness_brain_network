#!/bin/bash
#
# ==============================================================================
# Demean_Parcel_Cormat_v1.slurm
#
# Purpose:
#   Preprocesses HCP-Aging (HCA) resting-state fMRI data into subject-level
#   functional connectivity correlation matrices. Runs on an HPC cluster via
#   SLURM using Connectome Workbench (wb_command / wb_shortcuts).
#
# Pipeline steps (per subject):
#   1. Demean and normalize dtseries (wb_shortcuts -cifti-demean -normalize)
#   2. Parcellate dtseries using the Cole-Anticevic (CA) network atlas
#      (wb_command -cifti-parcellate)
#   3. Compute per-run Fisher-z correlation matrices
#      (wb_command -cifti-correlation -fisher-z)
#   4. Average correlation matrices across 4 resting-state runs
#      (wb_command -cifti-average)
#
# Inputs:
#   - HCA resting-state dtseries:
#       <hca_dir>/HCA<subj>_V1_MR/MNINonLinear/Results/
#           rfMRI_REST{1,2}_{AP,PA}_Atlas_MSMAll_hp0_clean.dtseries.nii
#   - Cole-Anticevic parcellation atlas (dlabel.nii):
#       <atlas_dir>/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii
#   - Subject list CSV 
#
# Outputs (written to <output_dir>/):
#   - subj_demean/HCA<subj>/    : demeaned/normalized dtseries per run
#   - subj_parcel/HCA<subj>/    : parcellated ptseries per run
#   - subj_cormatPerRun/HCA<subj>/ : per-run Fisher-z correlation matrices
#   - subj_Avecormat/           : run-averaged correlation matrix per subject
#                                  (<num_id>_rest_CormatAveRuns.pconn.nii)
#
# Dependencies:
#   - Connectome Workbench >= 1.4.2  (module: connectome-workbench/1.4.2)
#
# Usage:
#   Edit the USER CONFIGURATION section below, then submit with:
#     sbatch Demean_Parcel_Cormat_v1.slurm
# ==============================================================================

#SBATCH --job-name=demean_parcel_cormat
#SBATCH --time=12:00:00
#SBATCH --nodes=8
#SBATCH --output=demean_parcel_cormat.txt
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=your.email@institution.edu
#SBATCH -p medium-28core


# ==============================================================================
# USER CONFIGURATION — update these paths before running
# ==============================================================================

# Root directory containing downloaded HCA fmriresults01 data
hca_dir="/path/to/HCA_data/fmriresults01"

# Root directory for all pipeline outputs
output_dir="/path/to/HCA_stat"

# Directory containing the Cole-Anticevic parcellation atlas
atlas_dir="/path/to/atlas"

# CSV file listing subjects to process (must include a numeric ID column)
subj_list="/path/to/HCA_subject_list.csv"

# Row range (0-based) of subjects to process in this job
# Adjust start/end to distribute subjects across SLURM array jobs
subj_start=128
subj_end=138   # exclusive upper bound


# ==============================================================================
# OUTPUT SUBDIRECTORIES (relative to output_dir — do not edit unless needed)
# ==============================================================================

demean_dir="subj_demean"
parcel_dir="subj_parcel"
cormat_dir="subj_cormatPerRun"
avecormat_dir="subj_Avecormat"


# ==============================================================================
# SETUP
# ==============================================================================

module load shared
module load connectome-workbench/1.4.2

# Locate the Cole-Anticevic parcellation dlabel file
atlas_dlabel=$(find "$atlas_dir" -name '*ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii')

# Read the numeric subject IDs from column 2 of the CSV (skip header row)
subjects=($(awk -F, 'NR>1 {print $2}' "$subj_list"))


# ==============================================================================
# MAIN LOOP — process each subject in the specified row range
# ==============================================================================

for (( i=subj_start; i<subj_end; i++ )); do

    # ------------------------------------------------------------------
    # Step 1: Demean and normalize dtseries
    # ------------------------------------------------------------------

    subj_out_demean="$output_dir/$demean_dir/HCA${subjects[$i]}"
    mkdir -p "$subj_out_demean"

    # Source resting-state dtseries directory for this subject
    subj_dtseries="$hca_dir/HCA${subjects[$i]}_V1_MR/MNINonLinear/Results"

    for filename in $(find "$subj_dtseries" -name 'rfMRI*REST*_*_Atlas_MSMAll*clean.dtseries.nii'); do
        f=$(basename "$filename")
        wb_shortcuts -cifti-demean -normalize \
            "$subj_out_demean/demean_normalize_HCA${subjects[$i]}_$f" \
            "$filename"
    done


    # ------------------------------------------------------------------
    # Step 2: Parcellate demeaned dtseries using the CA atlas
    # ------------------------------------------------------------------

    subj_out_parcel="$output_dir/$parcel_dir/HCA${subjects[$i]}"
    mkdir -p "$subj_out_parcel"

    for demeanfile in $(find "$subj_out_demean" -name 'demean_normalize*.dtseries.nii'); do
        demeanf=$(basename "$demeanfile")

        # Rename: demean_normalize_... -> parcel_CA_... (.dtseries -> .ptseries)
        outfile="${demeanf/demean_normalize/parcel_CA}"
        outfile="${outfile/dtseries/ptseries}"

        # -cifti-parcellate: parcellate along COLUMN direction
        wb_command -cifti-parcellate \
            "$demeanfile" \
            "$atlas_dlabel" \
            COLUMN \
            "$subj_out_parcel/$outfile"
    done


    # ------------------------------------------------------------------
    # Step 3: Compute per-run Fisher-z correlation matrices
    # ------------------------------------------------------------------

    subj_out_cormat="$output_dir/$cormat_dir/HCA${subjects[$i]}"
    mkdir -p "$subj_out_cormat"

    for parcelfile in $(find "$subj_out_parcel" -name 'parcel_CA*.ptseries.nii'); do
        ciftiin=$(basename "$parcelfile")

        # Rename: parcel_CA_... -> cormat_z_CA_... (.ptseries -> .pconn)
        outfile="${ciftiin/parcel/cormat_z}"
        outfile="${outfile/ptseries/pconn}"

        # -cifti-correlation with Fisher-z transform
        wb_command -cifti-correlation \
            "$parcelfile" \
            "$subj_out_cormat/$outfile" \
            -fisher-z
    done


    # ------------------------------------------------------------------
    # Step 4: Average correlation matrices across all 4 resting-state runs
    # ------------------------------------------------------------------

    out_avecormat="$output_dir/$avecormat_dir"
    mkdir -p "$out_avecormat"

    # Collect all per-run pconn files for this subject
    subj_cormat_files=()
    for cormatfile in $(find "$subj_out_cormat" -name 'cormat_z*.pconn.nii'); do
        subj_cormat_files+=("-cifti" "$cormatfile")
    done

    # Average across runs; output named by numeric subject ID
    wb_command -cifti-average \
        "$out_avecormat/${subjects[$i]}_rest_CormatAveRuns.pconn.nii" \
        "${subj_cormat_files[@]}"

done
