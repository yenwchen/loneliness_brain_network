% fmri_stats_palm.m
%
% Purpose:
%   Runs FSL PALM (Permutation Analysis of Linear Models) permutation tests
%   on subject-level graph-theory network metrics derived from ABSOLUTE edge
%   weights (dissertation version; main text results). Final sample N = 512.
%
%   PALM is called once per network metric. Each call runs 5000 permutations
%   with family-wise error correction across contrasts, two-tailed testing,
%   log p-values, and demeaning. Site is modeled as an exchangeability block.
%
% Model parameters:
%   Main predictors (in order): loneliness, age, loneliness x age interaction,
%                                sex, household income, marital status
%   Covariates: total GMV, mean FD, MoCA, site (4 dummy variables) + intercept
%
% Inputs (all paths relative to working directory set by run_fmri_stats_palm.slurm):
%   Nodal metric CSVs (subjects x parcels, n = 512):
%     abs_nodal_<metric>_n512.csv
%   Global metric CSVs (subjects x 1, n = 512):
%     abs_global_<metric>_n512.csv
%   Design matrix:          abs_design_matrix_n512.csv
%   T-contrast file:        abs_Tcontrast.csv
%   F-contrast file:        abs_Fcontrast.csv
%   Exchangeability blocks: abs_site_block_n512.csv
%
% Outputs (written to palm_output/<metric>/):
%   PALM output files: palm_abs_<metric>_*  (t-stats, FWE-p, etc.)
%
% Note:
%   Output directories (palm_output/<metric>/) must exist before running.
%   They are created automatically by run_fmri_stats_palm.slurm.
%
% Usage:
%   Run via SLURM: sbatch run_fmri_stats_palm.slurm
%   Or directly in MATLAB (from the correct working directory):
%     run('fmri_stats_palm.m')


% ------------------------------------------------------------------------------
% PALM flags used in all calls:
%   -vg auto      : automatically infer variance groups
%   -corrcon      : FWE correction across contrasts
%   -twotail      : two-tailed testing
%   -n 5000       : number of permutations
%   -logp         : output log10(p) instead of raw p
%   -demean       : demean the input data
%   -nouncorrected: suppress uncorrected p-value output
%   -eb           : exchangeability blocks (site)
% ------------------------------------------------------------------------------


% Betweenness centrality (nodal)
palm -i abs_nodal_betweencent_n512.csv -d abs_design_matrix_n512.csv -t abs_Tcontrast.csv -f abs_Fcontrast.csv -eb abs_site_block_n512.csv -vg auto -corrcon -twotail -n 5000 -logp -demean -nouncorrected -o palm_output/betweencent/palm_abs_betweencent

% Closeness centrality (nodal)
palm -i abs_nodal_closecent_n512.csv -d abs_design_matrix_n512.csv -t abs_Tcontrast.csv -f abs_Fcontrast.csv -eb abs_site_block_n512.csv -vg auto -corrcon -twotail -n 5000 -logp -demean -nouncorrected -o palm_output/closecent/palm_abs_closecent

% Strength (nodal)
palm -i abs_nodal_strength_n512.csv -d abs_design_matrix_n512.csv -t abs_Tcontrast.csv -f abs_Fcontrast.csv -eb abs_site_block_n512.csv -vg auto -corrcon -twotail -n 5000 -logp -demean -nouncorrected -o palm_output/strength/palm_abs_strength

% Normalized strength (nodal)
palm -i abs_nodal_normstrength_n512.csv -d abs_design_matrix_n512.csv -t abs_Tcontrast.csv -f abs_Fcontrast.csv -eb abs_site_block_n512.csv -vg auto -corrcon -twotail -n 5000 -logp -demean -nouncorrected -o palm_output/normstrength/palm_abs_normstrength

% Eigenvector centrality (nodal)
palm -i abs_nodal_eigencent_n512.csv -d abs_design_matrix_n512.csv -t abs_Tcontrast.csv -f abs_Fcontrast.csv -eb abs_site_block_n512.csv -vg auto -corrcon -twotail -n 5000 -logp -demean -nouncorrected -o palm_output/eigencent/palm_abs_eigencent

% PageRank centrality (nodal)
palm -i abs_nodal_pgcent_n512.csv -d abs_design_matrix_n512.csv -t abs_Tcontrast.csv -f abs_Fcontrast.csv -eb abs_site_block_n512.csv -vg auto -corrcon -twotail -n 5000 -logp -demean -nouncorrected -o palm_output/pgcent/palm_abs_pgcent

% Clustering coefficient (nodal)
palm -i abs_nodal_clustercoef_n512.csv -d abs_design_matrix_n512.csv -t abs_Tcontrast.csv -f abs_Fcontrast.csv -eb abs_site_block_n512.csv -vg auto -corrcon -twotail -n 5000 -logp -demean -nouncorrected -o palm_output/clustercoef/palm_abs_clustercoef

% Participation coefficient — CA network partition (nodal)
palm -i abs_nodal_partcoef_n512.csv -d abs_design_matrix_n512.csv -t abs_Tcontrast.csv -f abs_Fcontrast.csv -eb abs_site_block_n512.csv -vg auto -corrcon -twotail -n 5000 -logp -demean -nouncorrected -o palm_output/partcoef/palm_abs_partcoef

% Modularity — CA network partition (global, 1 value per subject)
palm -i abs_global_modularity_n512.csv -d abs_design_matrix_n512.csv -t abs_Tcontrast.csv -f abs_Fcontrast.csv -eb abs_site_block_n512.csv -vg auto -corrcon -twotail -n 5000 -logp -demean -nouncorrected -o palm_output/modularity/palm_abs_modularity

% Average shortest path length (global, 1 value per subject)
palm -i abs_global_average_shortest_path_n512.csv -d abs_design_matrix_n512.csv -t abs_Tcontrast.csv -f abs_Fcontrast.csv -eb abs_site_block_n512.csv -vg auto -corrcon -twotail -n 5000 -logp -demean -nouncorrected -o palm_output/average_shortest_path/palm_abs_average_shortest_path
