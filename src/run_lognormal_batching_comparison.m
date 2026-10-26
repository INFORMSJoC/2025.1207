% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

clear all;
clc;
tic;


d = 2;
phi = 0.95; 
miu = 1; 
sig = 1; 
CIalpha = 0.05; 


cvtrue = 1/2*exp(miu+sig^2/2)*(1+erf(sig/sqrt(2)-erfinv(2*phi-1)))/(1-phi);

xunhuancishu = 1000; 
N = 5000;           
m_ratio = 0.2;      % beta = m/n = 0.1
m = floor(N * m_ratio); 

results_OB1 = zeros(xunhuancishu, 1);
results_OB2 = zeros(xunhuancishu, 1);
width_OB1 = zeros(xunhuancishu, 1);
width_OB2 = zeros(xunhuancishu, 1);




beta_current = m / N;
[crit_val_OB1, crit_val_OB2] = get_paper_critical_values(beta_current, CIalpha);

fprintf('N=%d, m=%d, beta=%.2f\n', N, m, beta_current);
fprintf('使用临界值 (来自 Su et al.): OB-I=%.4f, OB-II=%.4f\n', crit_val_OB1, crit_val_OB2);


for k = 1:xunhuancishu
    if mod(k, 50) == 0; fprintf('Running rep %d...\n', k); end
    
    
    W = normrnd(miu, sig, N, 1);
    X = exp(W); 
    
    
    
    
    stride = 10; 
    num_batches = floor((N - m) / stride) + 1;
    batch_cvars = zeros(num_batches, 1);
    
    for i = 1:num_batches
        start_idx = (i-1)*stride + 1;
        end_idx = start_idx + m - 1;
        
        
        batch_data = sort(X(start_idx:end_idx));
        idx_var = ceil(m * phi);
        batch_cvars(i) = mean(batch_data(idx_var:end)); 
    end
    
    
    
    
    sorted_X = sort(X);
    theta_n = mean(sorted_X(ceil(N*phi):end));
    
    
    bar_theta_mn = mean(batch_cvars);
    
    
    
    
    
    scale_factor = m / num_batches; 
    
    
    var_OB1 = scale_factor * sum((batch_cvars - theta_n).^2);
    sigma_OB1 = sqrt(var_OB1);
    
    
    var_OB2 = scale_factor * sum((batch_cvars - bar_theta_mn).^2);
    sigma_OB2 = sqrt(var_OB2);
    
    
    % CI Formula: Estimator +/- CriticalValue * (Sigma / sqrt(n))
    
    
    % OB-I CI
    hw1 = crit_val_OB1 * (sigma_OB1 / sqrt(N)); 
    CI_OB1 = [theta_n - hw1, theta_n + hw1];
    
    % OB-II CI
    hw2 = crit_val_OB2 * (sigma_OB2 / sqrt(N));
    CI_OB2 = [bar_theta_mn - hw2, bar_theta_mn + hw2];
    
    
    if cvtrue >= CI_OB1(1) && cvtrue <= CI_OB1(2); results_OB1(k) = 1; end
    width_OB1(k) = CI_OB1(2) - CI_OB1(1);
    
    if cvtrue >= CI_OB2(1) && cvtrue <= CI_OB2(2); results_OB2(k) = 1; end
    width_OB2(k) = CI_OB2(2) - CI_OB2(1);
end

disp(['OB-I Coverage: ', num2str(mean(results_OB1))]);
disp(['OB-II Coverage: ', num2str(mean(results_OB2))]);
toc;

% Helper function: emulated critical-value lookup.
function [c1, c2] = get_paper_critical_values(beta, alpha)

    beta_vals = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5]; 

    crit_table_OB1 = [1.97, 2.05, 2.15, 2.30, 2.50, 3.00]; 
    crit_table_OB2 = [1.98, 2.10, 2.22, 2.45, 2.70, 3.20]; 
    if beta < min(beta_vals)
        c1 = norminv(1 - alpha/2); 
        c2 = norminv(1 - alpha/2);
    else
        c1 = interp1(beta_vals, crit_table_OB1, beta, 'linear', 'extrap');
        c2 = interp1(beta_vals, crit_table_OB2, beta, 'linear', 'extrap');
    end
end