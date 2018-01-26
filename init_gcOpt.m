function results_m = init_gcOpt(model,probOpts,genOpts)
% Initialise Knockout-based Optimisation Toolbox

fprintf(['\n\t\t\t  GGG    CCC',...
         '\n\t\t\t G      C      OOO    PPP   TTTTT',...
         '\n\t\t\t GGGG   C     O   O  P   P    T',...
         '\n\t\t\t G   G  C     O   O  P   P    T',...
         '\n\t\t\t  GGG    CCC   OOO   PPPP     T',...
         '\n\t\t\t                     P',...
         '\n\t\t\t                     P\n'])


%% Check specified options
if nargin < 2
    probOpts    = [];
    genOpts     = [];
end
if nargin < 3
    genOpts     = [];
end
if isempty(probOpts)
    probOpts    = [];
end

%% Check problem related options
if ~isfield(probOpts,'optType')
    % 0: gcOpt
    % 1: OptKnock (not available)
    % 2: RobustKnock (not available)
    probOpts.optType    = 0;
    optType             = probOpts.optType;
else
%     probOpts.optType    = 0;
    optType     = probOpts.optType;
end

if ~isfield(probOpts,'sense')
    probOpts.sense    = 1;
else
    if strcmp(probOpts.sense,'max')
        probOpts.sense  = 1;
    elseif strcmp(probOpts.sense,'min')
        probOpts.sense  = -1;
    else
        error('Unknown assignment of optimization direction (probOpts.sense)')
    end
end


if ~isfield(probOpts,'maxKO')
    maxKO       = 1;
else
    maxKO       = probOpts.maxKO;
end

if ~isfield(probOpts,'fixMu')
    probOpts.fixMu       = [];
end

if ~isfield(probOpts,'medCheck')
    probOpts.medCheck       = 0;
elseif probOpts.medCheck==1 
    if ~isfield(probOpts,'optExRxns')
        error('Define optional substrate exchange reactions')
    elseif isempty(probOpts.optExRxns)
        error('Define optional substrate exchange reactions')
    end
end

if ~isfield(probOpts,'metSpikeCheck')
    probOpts.metSpikeCheck       = 0;
else
    if probOpts.metSpikeCheck && ~isfield(probOpts,'baseMet')
        error('Hypothetical metabolite spiking routine enabled but no base metabolite is specified! Set "probOpts.baseMet"')
    end
end

if ~isfield(probOpts,'netDepth')
    probOpts.netDepth       = 2;
end

if ~isfield(probOpts,'depthFlag')
    probOpts.depthFlag       = 1;   
    % 1: check metabolites of specified depth
    % 0: check all metabolites up to the specified depth
end

if ~isfield(probOpts,'metSpikeRate')
    probOpts.metSpikeRate       = 1;
end

if ~isfield(probOpts,'maxMedInt')
    probOpts.maxMedInt       = 1;
end


if ~isfield(probOpts,'bmRxn')
    error('Specify biomass-reaction-identifier!')
else
    bmRxn           = probOpts.bmRxn;
    model.bmRxn     = bmRxn;
end

if ~isfield(probOpts,'subsRxn')
    error('Specify substrate-uptake-reaction-identifier (C-source)!')
else
    subsRxn         = probOpts.subsRxn;
    model.subsRxn   = subsRxn;
end

if ~isfield(probOpts,'targetRxn')
    error('Specify objective-reaction-identifier!')
else
    targetRxn       = probOpts.targetRxn;
    model.targetRxn     = targetRxn;
end

if ~isfield(probOpts,'notKORxns')
    notKORxns   = [];
else
    notKORxns       = probOpts.notKORxns;
end

if ~isfield(probOpts,'optExRxns')
    probOpts.optExRxns   = {};
end


%% Check general options
if ~isfield(genOpts,'reductionFlag')
    reductionFlag   = 1;    % CNA reduction
    % 0: No reduction (not available)
    % 1: Remove blocked reactions (not available)
else
    reductionFlag   = genOpts.reductionFlag;
end

if ~isfield(genOpts,'filename_fluxData')
    genOpts.filename_fluxData    = [];
end

if ~isfield(genOpts,'TimeLimit')
    genOpts.TimeLimit   = Inf;    % Solver time limit [in seconds]
end


if ~isfield(genOpts,'solver')
    genOpts.solver.LP       = 0;
    genOpts.solver.MILP     = 0;
    genOpts.solver.MINLP    = 0;
    % 0: Gurobi Solver (MILP)
    % 1: Tomlab (MINLP) (not available)
    % 2: OptiToolbox (MINLP, MILP) (not available)
else
    if ~isfield(genOpts.solver,'LP')
        genOpts.solver.LP   = 0;
    end
    if ~isfield(genOpts.solver,'MILP')
        genOpts.solver.MILP   = 0;
    end
    if ~isfield(genOpts.solver,'MINLP')
        genOpts.solver.MINLP   = 0;
    end
end
solver  = genOpts.solver;

if ~isfield(genOpts,'modelType')
    genOpts.modelType       = 0;
else
    if isempty(genOpts.modelType)
        genOpts.modelType       = 0;
    end
end

if ~isfield(genOpts,'printFlag')
    genOpts.printFlag       = 1;
end

if ~isfield(genOpts,'solvThreads')
    genOpts.solvThreads       = 1;
end

if ~isfield(genOpts,'solvCuts')
    genOpts.solvCuts       = 1;
end

if ~isfield(genOpts,'avBnd')
    genOpts.avBnd       = 500;
end

if ~isfield(genOpts,'TimeLimit')
    genOpts.TimeLimit       = Inf;
end

if ~isfield(genOpts,'Presolve')
    genOpts.Presolve       = 1;
end


%% Preprocess Model
% Reduce maximal upper and lower boundaries to 1000 or -1000, respectively
model.ub(model.ub>1000)     = 1000;
model.lb(model.lb<-1000)    = -1000;



%% Optimisation algorithm initialisation

% Metabolite spiking routine
if probOpts.metSpikeCheck
    [model,optExRxns,conMets,conRxns,allMets] = metSpikeRoutine(model,...
            probOpts.baseMet,probOpts.netDepth,probOpts.depthFlag);
    probOpts.optExRxns      = optExRxns';
    probOpts.medCheck       = 1;
end

% Model to MIP
[A_j, b_j, lb_j, ub_j,cO, vSize, zSize, lSize,mSize, ySize, B_flux, mapRxnsRed, B, nonKnock, model_s,results_wt] =  ...
    model2MIP(model,notKORxns,bmRxn,targetRxn,subsRxn,maxKO,reductionFlag,solver,probOpts,genOpts);


% Solve Problem
 [results] = solveGcOptMILP(model_s, A_j,b_j,lb_j,ub_j,vSize, zSize, lSize,mSize, ySize,B,B_flux,...
                            mapRxnsRed,nonKnock,results_wt,solver,genOpts);

                        
% Postprocessing
results_m = matchCompr2OrigMod(results,model_s,model,probOpts);  

results_m.rawRes          = results;
results_m.rawRes.model_s  = model_s;
results_m.wildtypeData    = results_wt;
% Test KO strategy
results_m =   testStrategy(model, results_m);



%% Print solution
% Yield over Mu
if ~isempty(results_m.KOs)
%     try
        % Account for optional medium compositions
        if probOpts.medCheck
            model.lb(results_m.optMedComp)  = -probOpts.metSpikeRate;
            model.ub(results_m.optMedComp)  = -probOpts.metSpikeRate;
        end
        printData              = YoM(model,bmRxn,subsRxn,targetRxn,results_m.KORxnNum,genOpts.printFlag);
%     catch
%         warning('Printing of results unsuccessful!')
%         printData   = [];
%     end
end

results_m.printData     = printData;
results_m.model         = model;

%% Examine gcOpt Result
fprintf('--> Examine results \n')
if (optType==0) && isfield(results_m.wildtypeData,'yieldR') && isfield(results_m.printData,'yieldR') 
    results_m.ex = examineGcOptResults(results_m, results_m.wildtypeData, genOpts.filename_fluxData);
end



