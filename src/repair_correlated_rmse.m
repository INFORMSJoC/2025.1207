% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

function R = repair_correlated_rmse(thetaguji,cvguji,Gcvguji,phi,fenbumoshi,cvtrue,opts)
%REPAIR_CORRELATED_RMSE Repair moving-target CVaR and scaled-gradient RMSE.
% Rows are independent optimization repetitions; optional columns are budgets.
% Paired theta, CVaR and gradient arrays must contain the same entries.
% The optimization recursion is not rerun or changed by this function.
% cvtrue is used ONLY for the separately labelled fixed-target ordinary RMSE.
if nargin<7, opts=struct(); end
validateattributes(thetaguji,{'numeric'},{'real','finite','nonempty'});
validateattributes(cvguji,{'numeric'},{'real','finite','nonempty'});
validateattributes(Gcvguji,{'numeric'},{'real','finite','nonempty'});
validateattributes(cvtrue,{'numeric'},{'real','finite','scalar'});
if numel(thetaguji)~=numel(cvguji) || numel(cvguji)~=numel(Gcvguji)
    error('Repair:Pairing','theta, CVaR and gradient must have equal paired sizes.');
end
if isvector(thetaguji)
    thetaguji=thetaguji(:); cvguji=cvguji(:); Gcvguji=Gcvguji(:);
elseif ~isequal(size(thetaguji),size(cvguji),size(Gcvguji))
    error('Repair:Pairing','For multiple budgets, arrays must have the same matrix shape.');
end
Ref=correlated_cvar_reference(thetaguji,phi,fenbumoshi,opts);
R=struct();
R.squaredErrorsPseudo=(cvguji-Ref.cvar).^2;
R.squaredErrorsGradient=(Gcvguji-Ref.scaledGradient).^2;
R.pseudoRMSE=sqrt(mean(R.squaredErrorsPseudo,1));
R.gradientPseudoRMSE=sqrt(mean(R.squaredErrorsGradient,1));
R.optimumRMSE=sqrt(mean((cvguji-cvtrue).^2,1));
R.independentOptimumRMSE=sqrt(mean((cvguji-Ref.optimumCVaR).^2,1));
R.fixedBenchmarkInput=cvtrue;
R.referenceCVaR=Ref.cvar;
R.referenceScaledGradient=Ref.scaledGradient;
R.referenceCVaRSE=Ref.cvarSE;
R.referenceScaledGradientSE=Ref.scaledGradientSE;
R.referenceOptimumCVaR=Ref.optimumCVaR;
R.referenceOptimumCVaRSE=Ref.optimumCVaRSE;
R.theta=thetaguji; R.cvEstimate=cvguji; R.scaledGradientEstimate=Gcvguji;
R.repetitions=size(thetaguji,1);
R.phi=phi; R.fenbumoshi=fenbumoshi;
B=Ref.settings.numScrambles;
rmseB=zeros(size(thetaguji,2),B);
grmseB=rmseB;
for b=1:B
    cb=reshape(Ref.batchCVaR(:,b),size(thetaguji));
    gb=reshape(Ref.batchScaledGradient(:,b),size(thetaguji));
    rmseB(:,b)=sqrt(mean((cvguji-cb).^2,1))';
    grmseB(:,b)=sqrt(mean((Gcvguji-gb).^2,1))';
end
R.pseudoRMSEReferenceSE=std(rmseB,0,2)'/sqrt(B);
R.gradientPseudoRMSEReferenceSE=std(grmseB,0,2)'/sqrt(B);
R.reference=Ref;
R.definition='sqrt(mean((C_hat-C_phi(theta))^2)); gradient target is (1-phi)*C_phi''(theta).';
R.uncertaintyNote='Reference SEs describe independent reference integration only; optimization sampling error and interpolation are separate.';
end
