% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

function result = simulate_correlated_case(cfg, N, seed)
%SIMULATE_CORRELATED_CASE One scalar replication of correlated-noise CVaR SA.
%   R = SIMULATE_CORRELATED_CASE(CFG,N,SEED) implements the recursions in
%   CVaR_SA_duowei_ex5_correlated_noise.m without allocating trajectory arrays.
%   N updates are performed, but the reported estimates use state N, after
%   N-1 updates, exactly as the original post-processing did. R.nextState is
%   state N+1. R.H is H after N-1 updates; R.HLastGenerated is H after N.
%
%   Required CFG fields (same names as the original script):
%     phi, CIalpha, fenbumoshi, k0, N1,
%     alpha,beta,gamma,gammau,varphi,eta,varsigma,
%     alphanc,betanc,gammanc,gammaunc,varphinc,etanc,cnc.
%   Optional scalar fields: theta0,v0,cv0,cu0,Gv0,Gcv0 (all default 1),
%     hn0 (0), omega (1), betah,alphah (1/3), LB (0), UB (2), cvtrue,
%     CIscale (1); use CIscale=1/sqrt(2) for the article's interval width.
%     zcritical (computed from CIalpha if absent), a precomputed normal
%     quantile; for CIalpha=.05 use 1.959963984540054.
%   CFG.noiseModel, when specified, must be 'original_garch'; this retains
%   the user's explicitly selected squared-Gaussian H recursion.
%   Optional CFG.steps contains arrays alphan,betan,gamman,gammaun,
%     varphin,etan,cn of length at least N, enabling reuse across replications.
%   Optional CFG.randomInputs contains zH,direction,Y,Yplus,Yminus (length N)
%     for deterministic equivalence testing. These are the ACTUAL innovations
%     and actual distribution draws, not uniform transform inputs.
%   Optional CFG.captureTrace=true returns all states in R.trace, columns
%     [theta,v,cv,cu,Gv,Gcv,H]. It is false in production to limit memory use.
%
%   Three independent Y draws are used at each iteration, including the plus
%   and minus perturbations. Bulk random-number generation has the same law
%   as the original draws, but is not intended to reproduce its RNG ordering.
%   A private RandStream leaves the caller's RNG state unchanged.
%   Confidence-interval variance and the absolute-value operation reproduce
%   the original script; this function does not validate that variance model.

timer = tic;
validateattributes(N, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(seed, {'numeric'}, {'scalar','integer','nonnegative','finite'});
required = {'phi','CIalpha','fenbumoshi','k0','N1', ...
    'alpha','beta','gamma','gammau','varphi','eta','varsigma', ...
    'alphanc','betanc','gammanc','gammaunc','varphinc','etanc','cnc'};
for j = 1:numel(required)
    if ~isfield(cfg, required{j})
        error('CVaR:MissingConfiguration','Missing cfg.%s.', required{j});
    end
    validateattributes(cfg.(required{j}), {'numeric'}, {'scalar','finite'});
end
assert(cfg.phi > 0 && cfg.phi < 1, 'CVaR:Phi','cfg.phi must lie in (0,1).');
assert(cfg.CIalpha > 0 && cfg.CIalpha < 1, 'CVaR:CIalpha','cfg.CIalpha must lie in (0,1).');
assert(any(cfg.fenbumoshi == [1 2]), 'CVaR:Distribution','fenbumoshi must be 1 or 2.');
validateattributes(cfg.k0, {'numeric'}, {'scalar','integer','nonnegative'});
assert(cfg.N1 >= 0 && cfg.cnc > 0, 'CVaR:Steps','N1 must be nonnegative and cnc positive.');

theta = optional(cfg,'theta0',1); v = optional(cfg,'v0',1);
cv = optional(cfg,'cv0',1); cu = optional(cfg,'cu0',1);
Gv = optional(cfg,'Gv0',1); Gcv = optional(cfg,'Gcv0',1);
H = optional(cfg,'hn0',0);
omega = optional(cfg,'omega',1);
betah = optional(cfg,'betah',1/3); alphah = optional(cfg,'alphah',1/3);
LB = optional(cfg,'LB',0); UB = optional(cfg,'UB',2);
assert(LB <= UB && H >= 0 && omega >= 0 && betah >= 0 && alphah >= 0, ...
    'CVaR:Model','Invalid bounds or nonnegative volatility-model parameters.');
phi = cfg.phi; tailProbability = 1-phi; k0 = cfg.k0;
CIscale = optional(cfg,'CIscale',1);
assert(CIscale > 0, 'CVaR:CIscale','CIscale must be positive.');
if isfield(cfg,'noiseModel')
    assert(strcmp(cfg.noiseModel,'original_garch'), 'CVaR:NoiseModel', ...
        'This kernel preserves the original squared-Gaussian H recursion.');
end
captureTrace = isfield(cfg,'captureTrace') && cfg.captureTrace;
if captureTrace
    trace = zeros(N+1,7);
end

if isfield(cfg,'steps')
    steps = cfg.steps;
    stepNames = {'alphan','betan','gamman','gammaun','varphin','etan','cn'};
    for j = 1:numel(stepNames)
        assert(isfield(steps,stepNames{j}) && numel(steps.(stepNames{j})) >= N, ...
            'CVaR:Steps','Each cfg.steps array must contain at least N entries.');
    end
    alphan = steps.alphan; betan = steps.betan; gamman = steps.gamman;
    gammaun = steps.gammaun; varphin = steps.varphin;
    etan = steps.etan; cn = steps.cn;
else
    ngrid = (1:N)'; shifted = ngrid+cfg.N1;
    alphan = cfg.alphanc ./ ngrid.^cfg.alpha;
    betan = cfg.betanc ./ shifted.^cfg.beta;
    gamman = cfg.gammanc ./ shifted.^cfg.gamma;
    gammaun = cfg.gammaunc ./ shifted.^cfg.gammau;
    varphin = cfg.varphinc ./ shifted.^cfg.varphi;
    etan = cfg.etanc ./ shifted.^cfg.eta;
    cn = cfg.cnc ./ shifted.^cfg.varsigma;
end

if isfield(cfg,'randomInputs')
    randomInputs = cfg.randomInputs;
    inputNames = {'zH','direction','Y','Yplus','Yminus'};
    for j = 1:numel(inputNames)
        assert(isfield(randomInputs,inputNames{j}) && numel(randomInputs.(inputNames{j})) >= N, ...
            'CVaR:RandomInputs','Each cfg.randomInputs array must contain at least N entries.');
    end
    zH = randomInputs.zH; direction = randomInputs.direction;
    Y = randomInputs.Y; Yplus = randomInputs.Yplus; Yminus = randomInputs.Yminus;
    assert(all(abs(direction(1:N)) == 1), 'CVaR:Direction','Directions must be +1 or -1.');
else
    stream = RandStream('mt19937ar','Seed',seed);
    zH = randn(stream,N,1);
    direction = 2*(rand(stream,N,1) >= 0.5)-1;
    if cfg.fenbumoshi == 1
        Y = -log(1-rand(stream,N,1));
        Yplus = -log(1-rand(stream,N,1));
        Yminus = -log(1-rand(stream,N,1));
    else
        Y = randn(stream,N,1); Yplus = randn(stream,N,1); Yminus = randn(stream,N,1);
    end
end

tailHistory = zeros(k0,1);
tailPrevious1 = 0; tailPrevious2 = 0;
for n = 1:N
    if captureTrace
        trace(n,:) = [theta,v,cv,cu,Gv,Gcv,H];
    end
    if n == N
        terminal = [theta,v,cv,cu,Gv,Gcv,H];
    end

    % Keep the original H update, including its operation order.
    un = H^(1/2)*zH(n);
    H = omega+betah*H+alphah*un^2;
    barcn = cn(n)/max([1,abs(Gv),abs(Gcv)]);
    perturbation = barcn*direction(n);
    X = Y(n)*((theta-1)*(theta-1)+2)+H;
    Xplus = Yplus(n)*((theta+perturbation-1)*(theta+perturbation-1)+2)+H;
    Xminus = Yminus(n)*((theta-perturbation-1)*(theta-perturbation-1)+2)+H;
    vplus = v+barcn*Gv*direction(n);
    vminus = v-barcn*Gv*direction(n);

    if X <= v
        vnew = v+betan(n)*(phi-1);
    else
        vnew = v+betan(n)*phi;
    end
    if X >= v
        tail = X;
    else
        tail = 0;
    end
    cvnew = cv+gamman(n)*(tail/tailProbability-cv);
    if k0 == 2
        varianceSample = tail^2+2*tail*(tailPrevious1+tailPrevious2);
        tailPrevious2 = tailPrevious1; tailPrevious1 = tail;
    else
        varianceSample = tail^2+2*tail*sum(tailHistory);
        for lag = k0:-1:2
            tailHistory(lag) = tailHistory(lag-1);
        end
        if k0 > 0
            tailHistory(1) = tail;
        end
    end
    cunew = cu+gammaun(n)*(varianceSample-cu);
    Gvnew = Gv-varphin(n)*((Xplus <= vplus)-(Xminus <= vminus))*direction(n)/(2*barcn);
    Gcvnew = Gcv+etan(n)*((Xplus*(Xplus >= vplus)-Xminus*(Xminus >= vminus))*direction(n)/(2*barcn)-Gcv);
    thetaCandidate = theta-alphan(n)*Gcv/tailProbability;
    if thetaCandidate >= LB && thetaCandidate <= UB
        thetanew = thetaCandidate;
    elseif thetaCandidate < LB
        thetanew = LB;
    else
        thetanew = UB;
    end

    theta = thetanew; v = vnew; cv = cvnew; cu = cunew;
    Gv = Gvnew; Gcv = Gcvnew;
end

result = stateStruct(terminal);
result.N = N; result.seed = seed; result.updates = N-1;
result.gamman = gamman(N);
result.varianceRaw = result.cu-(1+2*k0)*(result.cv*tailProbability)^2;
result.sigma = sqrt(abs(result.varianceRaw));
result.standardError = sqrt(result.gamman)*result.sigma/tailProbability;
result.CIstandardError = CIscale*result.standardError;
result.CIscale = CIscale;
if isfield(cfg,'zcritical')
    zcritical = optional(cfg,'zcritical',NaN);
    assert(zcritical > 0, 'CVaR:CriticalValue','zcritical must be positive.');
else
    zcritical = sqrt(2)*erfcinv(cfg.CIalpha);
end
result.CI = result.cv+[-1 1]*zcritical*result.CIstandardError;
if isfield(cfg,'cvtrue')
    result.CIN = (result.cv-cfg.cvtrue)/sqrt(result.gamman*abs(result.varianceRaw))*tailProbability/CIscale;
else
    result.CIN = NaN;
end
result.nextState = stateStruct([theta,v,cv,cu,Gv,Gcv,H]);
result.nextState.updates = N;
result.HLastGenerated = H;
result.elapsedSeconds = toc(timer);
if captureTrace
    trace(N+1,:) = [theta,v,cv,cu,Gv,Gcv,H];
    result.trace = trace;
    result.traceColumns = {'theta','v','cv','cu','Gv','Gcv','H'};
end
end

function value = optional(cfg,name,defaultValue)
if isfield(cfg,name)
    value = cfg.(name);
else
    value = defaultValue;
end
validateattributes(value, {'numeric'}, {'scalar','finite'}, '', ['cfg.' name]);
end

function state = stateStruct(values)
state = struct('theta',values(1),'v',values(2),'cv',values(3),'cu',values(4), ...
    'Gv',values(5),'Gcv',values(6),'H',values(7));
end
