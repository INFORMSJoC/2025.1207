% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

function cfg = correlated_paper_settings()
% Paper: PDF pp.23,43-45 (EC.11-13), 2026-09-14T231419.819.
% User decision: retain the ORIGINAL squared-normal, untruncated H model.
cfg.d = 1;
cfg.phi = 0.9; cfg.CIalpha = 0.05; cfg.CIscale = 1/sqrt(2);
cfg.zcritical = 1.959963984540054; % Phi^(-1)(0.975)
cfg.k0 = 2;
cfg.alpha = 0.999; cfg.beta = 0.998;
cfg.gamma = 0.997; cfg.gammau = 0.997;
cfg.varphi = 0.996; cfg.eta = 0.995; cfg.varsigma = 1/6;
cfg.alphanc = 0.8; cfg.betanc = 4000;
cfg.gammanc = 15; cfg.gammaunc = 15;
cfg.varphinc = 12000; cfg.etanc = 100; cfg.cnc = 2;
% Total run budget in shifted step sizes. Fixed across the four table columns.
cfg.N1 = 1000000;
cfg.Nvalues = [1000 10000 100000 1000000];
cfg.xunhuancishu = 100;
cfg.distributions = [1 2]; % exponential, normal
cfg.theta0 = 0; cfg.v0 = 1; cfg.cv0 = 1; cfg.cu0 = 1;
cfg.Gv0 = 1; cfg.Gcv0 = 1; cfg.hn0 = 0;
cfg.omega = 1; cfg.betah = 1/3; cfg.alphah = 1/3;
cfg.LB = 0; cfg.UB = 2; cfg.thetatrue = 1;
cfg.noiseModel = 'original_garch';
cfg.baseSeed = 20260915;
cfg.workers = 4;
cfg.referenceOptions = struct('noiseModel','original_garch', ...
    'M0',Inf,'samplesPower',20,'numScrambles',16, ...
    'depth',80,'seed',420260915,'gridSize',81, ...
    'hGridStep',0.002,'hGridMaximum',128,'sampleBlockPower',18);
end
