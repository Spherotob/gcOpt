%% gcOpt application Example

% gcOpt requires gurobi MILP solver, COBRA Toolbox. Make sure that
% all are installed and/or initialized before running the example or gcOpt in general.


%% Load and setup Model

% iAF1260 core model
load('iAF1260Core.mat')

model = changeRxnBounds(model,'EX_o2_e',0,'l');    %  anaerobic condition
model = changeRxnBounds(model,'EX_glc__D_e',-20,'l'); % Constrain maximal glucose consumption rate

notknockable        = [];
notknockable        = [notknockable; find(not(cellfun('isempty',(strfind(model.rxnNames,'diffusion')))))];
notknockable        = [notknockable; find(not(cellfun('isempty',(strfind(model.rxnNames,'exchange')))))];
% ATPM
notknockable        = [notknockable; 11];
probOpts.notKORxns  = model.rxns(notknockable);

% assign biomass formation, substrate uptake and target reaction
probOpts.bmRxn              = 'BIOMASS_Ecoli_core_w_GAM';
probOpts.subsRxn            = 'EX_glc__D_e';
probOpts.targetRxn          = 'EX_succ_e';


%% Additional or optional properties/options regarding the problem formulation


% Optimization direction
% 'max': maximization of objective reaction
% 'min': minimization of objective reaction
probOpts.sense  = 'max';

% Maximal allowable reaction deletions
probOpts.maxKO      = 7;

% Manually fix growth rate at which gcOpt optimizes the minimal guaranteed
% production rate (leave empty if growth rate should be automatically fixed
% at maximal productivity)
probOpts.fixMu      = 0.15;


%% General options (Mostly regarding the gurobi solver)
% "genOpts" are optional, default values are assigned here
genOpts.solvThreads     = 7;    % Number of parallel threads
genOpts.solvCuts        = 1;    % MILP Cuts setting
genOpts.avBnd           = 500;  % Maximal upper bound of auxiliary variables (dualisation)
genOpts.TimeLimit       = Inf;  % Maximal runtime of the solver [sec]
genOpts.Presolve        = 1;    % Presolver of gurobi solver

% Network compression
genOpts.compressFlag   = 1;
% 0: No compression
% 1: Remove blocked reactions

% model type/ specifier type (0: BIGG database)
genOpts.modelType   = 0;

% Specifiy Excel File name including flux data for reference flux
% distribution
genOpts.filename_fluxData   = [];



%% Initialize gcOpt
results = init_gcOpt(model,probOpts,genOpts);


%% Save results to a ".mat" file
% save('results_example_gcOpt','results')

