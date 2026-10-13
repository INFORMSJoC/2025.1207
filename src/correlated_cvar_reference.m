% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

function Ref = correlated_cvar_reference(theta, phi, fenbumoshi, opts)
%CORRELATED_CVAR_REFERENCE Independent reference for the ORIGINAL Z^2 model.
% X(theta)=((theta-1)^2+2)*Y+H, H_next=1+(1+Z^2)*H/3.
% Y is Exp(1) (fenbumoshi=1) or N(0,1) (fenbumoshi=2); Z is N(0,1).
% H has its stationary marginal and is independent of Y. No noise truncation.
% Ref.scaledGradient is (1-phi)*C_phi'(theta), the algorithm's scale.
%
% Independent randomized Sobol replicates approximate H; Y is integrated out
% analytically. The identity E[H]=3 is used as a control variate, so an
% empirical heavy-tail mean is NOT used to estimate CVaR. Standard errors are
% across independent scrambles. Grid interpolation and finite-depth effects
% are separately reported and are not included in those standard errors.
%
% opts fields (all optional): samplesPower=16, numScrambles=12, depth=80,
% seed=20260915, gridSize=81, cacheFile='', verbose=true,
% noiseModel='original_garch', M0=Inf, hGridStep=.002, hGridMaximum=128,
% sampleBlockPower=18, stationaryCacheFile=''. The stationary cache is shared
% between the two Y distributions; it contains compressed weighted samples.
% Requires Statistics and Machine Learning Toolbox (sobolset and scramble).

if nargin < 4 || isempty(opts), opts=struct(); end
opts=fill_defaults(opts);
validateattributes(theta,{'numeric'},{'real','finite','nonempty'});
validateattributes(phi,{'numeric'},{'real','scalar','>',0,'<',1});
validateattributes(fenbumoshi,{'numeric'},{'scalar','integer','>=',1,'<=',2});
if ~strcmpi(opts.noiseModel,'original_garch') || ~isinf(opts.M0)
    error('Reference:Model','This reference implements the original Z^2 model without clipping.');
end
validateattributes(opts.samplesPower,{'numeric'},{'scalar','integer','>=',8,'<=',23});
validateattributes(opts.numScrambles,{'numeric'},{'scalar','integer','>=',2});
validateattributes(opts.depth,{'numeric'},{'scalar','integer','>=',16});
validateattributes(opts.gridSize,{'numeric'},{'scalar','integer','>=',5});
validateattributes(opts.hGridStep,{'numeric'},{'scalar','real','positive','finite'});
validateattributes(opts.hGridMaximum,{'numeric'},{'scalar','real','>',3,'finite'});
validateattributes(opts.sampleBlockPower,{'numeric'},{'scalar','integer','>=',8,'<=',20});
if mod(opts.gridSize,2)~=1, error('Reference:Grid','gridSize must be odd.'); end
if exist('sobolset','file')~=2
    error('Reference:Toolbox','Statistics and Machine Learning Toolbox is required.');
end
inputSize=size(theta);
thetaVector=theta(:);
scale=(thetaVector-1).^2+2;
settings=opts;
settings=rmfield(settings,{'cacheFile','verbose','stationaryCacheFile'});
settings.phi=phi;
settings.fenbumoshi=fenbumoshi;
settings.scaleMinimum=2;
settings.scaleMaximum=max(3,max(scale));
settings.referenceVersion=2;
cacheUsed=false;
if ~isempty(opts.cacheFile) && isfile(opts.cacheFile)
    saved=load(opts.cacheFile,'Grid','settings');
    if isfield(saved,'settings') && isequaln(saved.settings,settings) && isfield(saved,'Grid')
        Grid=saved.Grid;
        cacheUsed=true;
    end
end
if ~cacheUsed
    Grid=build_grid(settings,opts.verbose,opts.stationaryCacheFile);
    if ~isempty(opts.cacheFile)
        cacheParent=fileparts(opts.cacheFile);
        if ~isempty(cacheParent) && ~isfolder(cacheParent), mkdir(cacheParent); end
        save(opts.cacheFile,'Grid','settings','-v7.3');
    end
end
B=opts.numScrambles;
cvarB=zeros(numel(scale),B);
qB=cvarB;
tailMomentB=cvarB;
coarseCvarB=cvarB;
coarseTailB=cvarB;
coarse=1:2:numel(Grid.scale);
for b=1:B
    cvarB(:,b)=interp1(Grid.scale,Grid.cvar(:,b),scale,'pchip');
    qB(:,b)=interp1(Grid.scale,Grid.quantile(:,b),scale,'pchip');
    tailMomentB(:,b)=interp1(Grid.scale,Grid.tailYMoment(:,b),scale,'pchip');
    coarseCvarB(:,b)=interp1(Grid.scale(coarse),Grid.cvar(coarse,b),scale,'pchip');
    coarseTailB(:,b)=interp1(Grid.scale(coarse),Grid.tailYMoment(coarse,b),scale,'pchip');
end
gradientB=2*(thetaVector-1).*tailMomentB;
Ref=struct();
Ref.theta=theta;
Ref.cvar=reshape(mean(cvarB,2),inputSize);
Ref.quantile=reshape(mean(qB,2),inputSize);
Ref.scaledGradient=reshape(mean(gradientB,2),inputSize);
Ref.cvarSE=reshape(std(cvarB,0,2)/sqrt(B),inputSize);
Ref.quantileSE=reshape(std(qB,0,2)/sqrt(B),inputSize);
Ref.scaledGradientSE=reshape(std(gradientB,0,2)/sqrt(B),inputSize);
Ref.cvarInterpolationDiagnostic=reshape(abs(mean(cvarB-coarseCvarB,2)),inputSize);
Ref.gradientInterpolationDiagnostic=reshape(abs(2*(thetaVector-1).*mean(tailMomentB-coarseTailB,2)),inputSize);
% With stationary Var(H)=6, coupling H_0=3 gives this rigorous CVaR depth bound.
Ref.cvarFiniteDepthBound=sqrt(6)*(2/3)^opts.depth/(1-phi);
% Linear interpolation of the conditional lower partial moment has curvature
% at most 1/s (normal density is smaller). CVaR is the minimum over q of this
% objective, so the same uniform bound applies after solving the quantile.
Ref.cvarHDiscretizationBound=opts.hGridStep^2/(8*min(scale)*(1-phi));
Ref.optimumCVaR=mean(Grid.cvar(1,:));
Ref.optimumCVaRSE=std(Grid.cvar(1,:),0,2)/sqrt(B);
Ref.batchCVaR=cvarB;
Ref.batchScaledGradient=gradientB;
Ref.Grid=Grid;
Ref.settings=settings;
Ref.cacheUsed=cacheUsed;
Ref.method='Conditional analytic Y integration; randomized Sobol stationary H; exact E[H]=3 control variate.';
Ref.standardErrorMeaning='Reference numerical uncertainty across independent scrambles, excluding interpolation; not the optimization experiment standard deviation.';
Ref.interpolationDiagnosticMeaning='Absolute difference of mean fine-grid and nested coarse-grid pchip values; a convergence diagnostic, not a rigorous bound.';
end

function opts=fill_defaults(opts)
defaults=struct('samplesPower',16,'numScrambles',12,'depth',80,...
    'seed',20260915,'gridSize',81,'cacheFile','','verbose',true,...
    'noiseModel','original_garch','M0',Inf,'hGridStep',.002,...
    'hGridMaximum',128,'sampleBlockPower',18,'stationaryCacheFile','');
names=fieldnames(defaults);
for j=1:numel(names)
    if ~isfield(opts,names{j}) || isempty(opts.(names{j}))
        opts.(names{j})=defaults.(names{j});
    end
end
unknown=setdiff(fieldnames(opts),names);
if ~isempty(unknown), error('Reference:Option','Unknown option: %s',unknown{1}); end
end

function Grid=build_grid(settings,verbose,stationaryCacheFile)
previousRng=rng;
cleanup=onCleanup(@() rng(previousRng)); %#ok<NASGU>
rng(settings.seed,'twister');
S=linspace(settings.scaleMinimum,settings.scaleMaximum,settings.gridSize)';
B=settings.numScrambles;
N=2^settings.samplesPower;
Grid=struct('scale',S,'cvar',zeros(numel(S),B),...
    'quantile',zeros(numel(S),B),'tailYMoment',zeros(numel(S),B),...
    'samplesPerScramble',N,'numScrambles',B);
ticReference=tic;
stationarySettings=rmfield(settings,{'phi','fenbumoshi','scaleMinimum','scaleMaximum','gridSize'});
stationaryCacheUsed=false;
if ~isempty(stationaryCacheFile) && isfile(stationaryCacheFile)
    savedStationary=load(stationaryCacheFile,'stationarySettings','stationarySamples');
    if isfield(savedStationary,'stationarySettings') && isequaln(savedStationary.stationarySettings,stationarySettings)
        stationarySamples=savedStationary.stationarySamples;
        stationaryCacheUsed=true;
    end
end
if ~stationaryCacheUsed
    stationarySamples=cell(B,1);
end
for b=1:B
    if ~stationaryCacheUsed
        stationarySamples{b}=stationary_sample(settings,N);
    end
    H=stationarySamples{b}.nodes;
    weights=stationarySamples{b}.weights;
    previousQ=NaN;
    for j=1:numel(S)
        s=S(j);
        if isnan(previousQ)
            bracket=[-20*s,3+30*s];
        else
            bracket=[previousQ-2,previousQ+2];
        end
        fun=@(q) tail_probability(q,H,weights,s,settings.fenbumoshi)-(1-settings.phi);
        while fun(bracket(1))<0, bracket(1)=bracket(1)-10*s; end
        while fun(bracket(2))>0, bracket(2)=bracket(2)+10*s; end
        q=fzero(fun,bracket,optimset('TolX',1e-10,'Display','off'));
        [lowerPartial,tailYMoment]=conditional_moments(q,H,s,settings.fenbumoshi);
        meanY=double(settings.fenbumoshi==1);
        Grid.cvar(j,b)=q+(3+s*meanY-q+sum(weights.*lowerPartial))/(1-settings.phi);
        Grid.quantile(j,b)=q;
        Grid.tailYMoment(j,b)=sum(weights.*tailYMoment);
        previousQ=q;
    end
    if verbose
        fprintf('Reference distribution %d, scramble %d/%d: C(1)=%.9f, elapsed %.1fs\n',...
            settings.fenbumoshi,b,B,Grid.cvar(1,b),toc(ticReference));
    end
end
if ~stationaryCacheUsed && ~isempty(stationaryCacheFile)
    cacheParent=fileparts(stationaryCacheFile);
    if ~isempty(cacheParent) && ~isfolder(cacheParent), mkdir(cacheParent); end
    save(stationaryCacheFile,'stationarySettings','stationarySamples','-v7.3');
end
Grid.elapsedSeconds=toc(ticReference);
Grid.stationaryCacheUsed=stationaryCacheUsed;
Grid.meanPreservationMaxAbs=max(cellfun(@(z) z.meanPreservationError,stationarySamples));
Grid.cvarMean=mean(Grid.cvar,2);
Grid.cvarSE=std(Grid.cvar,0,2)/sqrt(B);
Grid.tailYMomentMean=mean(Grid.tailYMoment,2);
Grid.tailYMomentSE=std(Grid.tailYMoment,0,2)/sqrt(B);
end

function compressed=stationary_sample(settings,N)
pointSet=scramble(sobolset(settings.depth),'MatousekAffineOwen');
gridNodes=(1.5:settings.hGridStep:settings.hGridMaximum)';
counts=zeros(size(gridNodes));
tailNodes=[];
meanH=0;
blockSize=2^min(settings.samplesPower,settings.sampleBlockPower);
for first=1:blockSize:N
    last=min(first+blockSize-1,N);
    U=pointSet(first:last,:);
    % |Z| inverse CDF removes an unnecessary symmetry from the integration.
    H=3*ones(last-first+1,1);
    for d=settings.depth:-1:1
        ud=min(max(U(:,d),eps),1-eps);
        H=1+((1+2*erfinv(ud).^2)/3).*H;
    end
    meanH=meanH+sum(H)/N;
    above=H>gridNodes(end);
    tailNodes=[tailNodes;H(above)]; %#ok<AGROW>
    H=H(~above);
    position=(H-1.5)/settings.hGridStep;
    left=min(floor(position)+1,numel(gridNodes)-1);
    fraction=position-(left-1);
    counts=counts+accumarray(left,1-fraction,size(counts))+...
        accumarray(left+1,fraction,size(counts));
end
occupied=counts>0;
compressed.nodes=[gridNodes(occupied);tailNodes];
compressed.weights=[counts(occupied);ones(size(tailNodes))]/N;
compressed.meanPreservationError=abs(sum(compressed.nodes.*compressed.weights)-meanH);
assert(abs(sum(compressed.weights)-1)<1e-10,'Reference weights do not sum to one.');
assert(compressed.meanPreservationError<1e-9,'Reference compression did not preserve the sample mean.');
end

function tail=tail_probability(q,H,weights,s,mode)
a=(q-H)/s;
if mode==1
    tail=sum(weights.*exp(-max(a,0)));
else
    tail=sum(weights.*(0.5*erfc(a/sqrt(2))));
end
end

function [lowerPartial,tailYMoment]=conditional_moments(q,H,s,mode)
a=(q-H)/s;
if mode==1
    ap=max(a,0);
    ea=exp(-ap);
    lowerPartial=s*(ap+expm1(-ap));
    tailYMoment=(ap+1).*ea;
else
    density=exp(-0.5*a.^2)/sqrt(2*pi);
    lowerPartial=s*(a.*(0.5*erfc(-a/sqrt(2)))+density);
    tailYMoment=density;
end
end
