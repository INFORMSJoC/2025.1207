% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

function results = run_correlated_paper_tables(cfg, outputDir)
%RUN_CORRELATED_PAPER_TABLES Run all 8 cases with paired endpoints.
% Each N column has 100 independent replications. N1 is the fixed maximum
% budget. References use independent randomness.
if ~isfolder(outputDir), mkdir(outputDir); end
started = tic;
tag = 'CVaR_correlated_k0_2_R100_corrected';
rawPath = fullfile(outputDir,[tag '_raw.mat']);
resume=[];
checkpointPath=fullfile(outputDir,[tag '_checkpoint.mat']);
if isfield(cfg,'resume') && cfg.resume && isfile(checkpointPath)
    resume=load(checkpointPath);
    checkNames={'N1','Nvalues','xunhuancishu','baseSeed','distributions', ...
        'alpha','beta','gamma','gammau','varphi','eta','varsigma', ...
        'alphanc','betanc','gammanc','gammaunc','varphinc','etanc','cnc', ...
        'theta0','v0','cv0','cu0','Gv0','Gcv0','hn0','k0','phi','CIalpha', ...
        'CIscale','noiseModel','omega','betah','alphah','LB','UB'};
    for z=1:numel(checkNames)
        assert(isequaln(cfg.(checkNames{z}),resume.cfg.(checkNames{z})), ...
            'Checkpoint settings differ: %s.',checkNames{z});
    end
end
% Back up existing output files before starting a new run.
if isfile(rawPath) && isempty(resume)
    backup = fullfile(outputDir,['backup_' tag '_' datestr(now,'yyyymmdd_HHMMSS')]);
    mkdir(backup);
    oldFiles = dir(fullfile(outputDir,[tag '*']));
    for f=1:numel(oldFiles)
        if ~oldFiles(f).isdir
            copyfile(fullfile(outputDir,oldFiles(f).name),backup);
        end
    end
end
fprintf('MODEL: original Z^2 H recursion, no truncation.\n');
fprintf('k0=%d; replications=%d; N1=%d; alpha=%.3f; alphanc=%.1f; theta0=%.1f\n', ...
    cfg.k0,cfg.xunhuancishu,cfg.N1,cfg.alpha,cfg.alphanc,cfg.theta0);
fprintf('Columns N: %s\n',mat2str(cfg.Nvalues));
pool = [];
if cfg.workers > 1 && license('test','Distrib_Computing_Toolbox')
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('Threads',cfg.workers); end
end
cases = struct([]);
if ~isempty(resume), cases=resume.cases; end
for distribution = cfg.distributions
    localCfg = cfg;
    localCfg.fenbumoshi = distribution;
    if distribution==1
        name = 'exponential'; localCfg.cvtrue = 11.6089;
    else
        name = 'normal'; localCfg.cvtrue = 8.8803;
    end
    for column=1:numel(cfg.Nvalues)
        N = cfg.Nvalues(column);
        if ~isempty(cases) && any([cases.fenbumoshi]==distribution & [cases.N]==N)
            continue
        end
        localCfg.steps = makeSteps(localCfg,N);
        replications = cell(cfg.xunhuancishu,1);
        fprintf('START %s N=%d R=%d\n',name,N,cfg.xunhuancishu);
        caseTimer=tic;
        firstPending=1;
        previousElapsed=0;
        if ~isempty(resume) && resume.distribution==distribution && resume.N==N
            replications=resume.replications;
            firstPending=resume.last+1;
            if isfield(resume,'caseElapsed'), previousElapsed=resume.caseElapsed; end
            fprintf('RESUME %s N=%d at replication %d\n',name,N,firstPending);
        end
        % Batches provide durable checkpoints and readable progress.
        for first=firstPending:10:cfg.xunhuancishu
            last=min(first+9,cfg.xunhuancishu);
            batch=cell(last-first+1,1);
            if ~isempty(pool)
                parfor j=1:numel(batch)
                    r=first+j-1;
                    seed=cfg.baseSeed+distribution*1000000+column*10000+r;
                    batch{j}=simulate_correlated_case(localCfg,N,seed);
                end
            else
                for j=1:numel(batch)
                    r=first+j-1;
                    seed=cfg.baseSeed+distribution*1000000+column*10000+r;
                    batch{j}=simulate_correlated_case(localCfg,N,seed);
                end
            end
            replications(first:last)=batch;
            caseElapsed=previousElapsed+toc(caseTimer);
            save(checkpointPath, ...
                'cfg','cases','replications','distribution','N','last','caseElapsed','-v7');
            fprintf('PROGRESS %s N=%d %d/%d elapsed=%.1fs\n', ...
                name,N,last,cfg.xunhuancishu,toc(caseTimer));
        end
        c=struct('distribution',name,'fenbumoshi',distribution,'N',N, ...
            'cvtrueOriginal',localCfg.cvtrue,'replications',[replications{:}], ...
            'simulationWallSeconds',previousElapsed+toc(caseTimer));
        if isempty(cases), cases=c;
        else, cases(end+1)=c; end %#ok<AGROW>
        save(rawPath,'cfg','cases','-v7');
    end
end
fprintf('ALL_OPTIMIZATION_COMPLETE %.1fs\n',toc(started));
if isfield(cfg,'deferReference') && cfg.deferReference
    results=struct('cfg',cfg,'cases',cases,'simulationWallSeconds',toc(started));
    return
end
results = finalize_correlated_tables(cfg,cases,outputDir);
results.totalWallSeconds=toc(started);
save(fullfile(outputDir,[tag '_results.mat']),'results','cfg','cases','-v7');
export_correlated_matlab_excel(results,outputDir);
fprintf('ALL_COMPLETE %.1fs; saved in %s\n',toc(started),outputDir);
end

function steps=makeSteps(c,N)
n=(1:N)'; shifted=n+c.N1;
steps=struct('alphan',c.alphanc./n.^c.alpha, ...
    'betan',c.betanc./shifted.^c.beta,'gamman',c.gammanc./shifted.^c.gamma, ...
    'gammaun',c.gammaunc./shifted.^c.gammau,'varphin',c.varphinc./shifted.^c.varphi, ...
    'etan',c.etanc./shifted.^c.eta,'cn',c.cnc./shifted.^c.varsigma);
end
