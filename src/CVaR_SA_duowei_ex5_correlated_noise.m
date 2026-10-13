% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

% Correlated-noise experiment: corrected CVaR and Gcv pseudo-RMSE.
% Paper settings: alpha=.999, alpha coefficient=.8, theta0=0, k0=2.
% 100 independent replications at N=1e3,1e4,1e5,1e6, for BOTH Y laws.
% Per user instruction, H retains the original Z^2 recursion, no truncation.
% See correlated_paper_settings.m for every parameter and step-size constant.
% The optimization recursions use the original state-N endpoint convention.
% Output files are written beside this script; prior outputs are backed up.
clearvars;
clc;
outputDir = fileparts(mfilename('fullpath'));
addpath(outputDir);
cfg = correlated_paper_settings();
results = run_correlated_paper_tables(cfg, outputDir);
