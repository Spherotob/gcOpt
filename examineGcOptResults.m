%% Examine gcOpt results struct for GC mutants
% Maximal growth rate
% Minimal production rate at maximal growth (Ymin)
% Ratio between pRmin and maximal theoretical production rate
% Minimal productivity at maximal growth
% Ratio between Pmin and maximal theoretical productivity
% GCS
% Growth state according to MiMBl

function ex = examineGcOptResults(results, results_wt, filename_fluxData)

if nargin < 3
    filename_fluxData    = [];
end

% if gcOpt results are invalid stop examination
if results.validity ~= 1
    ex{1}.gcs   = -2;
    return
end

rxnNumBM                = find(ismember(results.model.rxns,results.model.bmRxn));
rxnNumSubs              = find(ismember(results.model.rxns,results.model.subsRxn));
rxnNumTarget            = find(ismember(results.model.rxns,results.model.targetRxn));

yRange_wt      = results_wt.yieldR;
mRange_wt      = results_wt.muR;
pRange_wt      = results_wt.prodR;

yRange_mut       = results.printData.yieldR;
mRange_mut       = results.printData.muR;
pRange_mut       = results.printData.prodR;


%% Wild type
[th_maxY, ~]    = max(yRange_wt);
[th_maxP, ~]    = max(pRange_wt);


ex  = cell(2,1);
compSave    = [];
%% Maximal growth rate
[ex{1}.maxMu, ~]   = max(mRange_mut);
compSave           = [compSave,ex{1}.maxMu];

%% Minimal production rate at maximal growth (minPR)
isMax               = yRange_mut(mRange_mut==ex{1}.maxMu);
ex{1}.minPR_maxMu    = min(isMax);
compSave            = [compSave,ex{1}.minPR_maxMu];

%% Ratio between PRmin and maximal theoretical production rate [%]
ex{1}.ratio_minPR   = (ex{1}.minPR_maxMu/th_maxY)*100;
compSave           = [compSave,ex{1}.ratio_minPR];

%% Minimal productivity at maximal growth (Pmin)
isMax           = pRange_mut(mRange_mut==ex{1}.maxMu);
ex{1}.minP_maxMu   = min(isMax);
compSave           = [compSave,ex{1}.minP_maxMu];

%% Ratio between Pmin and maximal theoretical productivity [%]
ex{1}.ratio_minP   = (ex{1}.minP_maxMu/th_maxP)*100;
compSave           = [compSave,ex{1}.ratio_minP];

%% Growth coupling strength
ex{1}.gcs  = calcGCS(results, results_wt);
compSave           = [compSave,ex{1}.gcs];

%% Yield at maximal growth
model   = results.model;
% apply KOs
model   = changeRxnBounds(model,model.rxns(results.KORxnNum(:,1)),0,'b');
% FBA
model   = changeObjective(model,model.bmRxn);
sol     = optimizeCbModel(model,'max','one');
ex{1}.Y_maxMu   = sol.x(rxnNumTarget)/(-sol.x(rxnNumSubs));

%% Growth state according to MiMBl

model                   = results.model;
% Calculate reference flux distribution from data
opt.filename    = filename_fluxData;
opt.fluxFac     = 1;
if isempty(filename_fluxData)
    [model,~,~]     = createRefFD(model,[],0,opt);
else
    [model,~,~]     = createRefFD(model,[],1,opt);
end
refFluxDist         = model.fd_ref;
ex{3}.refFluxDist   = refFluxDist;
model_mut   = model;
% apply reaction deletions
model_mut   = changeRxnBounds(model_mut,model_mut.rxns(results.KORxnNum(:,1)),0,'b');
% MiMBl
solMiMBl = MiMBl(model_mut,refFluxDist,0);
ex{3}.mutFluxDist   = solMiMBl.x;
% extract yield/mu pair
ex{1}.mimbl_mu  = solMiMBl.x(rxnNumBM); 
ex{1}.mimbl_Y   = solMiMBl.x(rxnNumTarget)/-solMiMBl.x(rxnNumSubs);
ex{1}.mimbl_PR  = solMiMBl.x(rxnNumTarget);  
compSave           = [compSave,ex{1}.mimbl_mu, ex{1}.mimbl_Y, ex{1}.mimbl_PR];

% compact structure
ex{2}   = compSave;
end