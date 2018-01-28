%% gcOpt example employing the Escherichia coli iAF1260core model
%
%
% The metabolic model (iAF1260Core.mat) was provided and downloaded from the BiGG database, hosted and 
% maintained by the Systems Biology Research Group at the University of California, San Diego:
% "http://bigg.ucsd.edu/"
%
% Compatibility notes:
% - gcOpt requires gurobi MILP solver. Refer to "http://www.gurobi.com/"
%   for further information.
%       -> Academic licenses are freely available.
% 
%  - COBRA toolbox need to be installed. Cf.: "https://opencobra.github.io/cobratoolbox/stable/"
%
%
% Refer to init_gcOpt.m header for further information about the inputs/outputs of gcOpt
%
% AUTHOR:
%     - Tobias Alter      29th January 2018
%         Institute of Applied Microbiology, RWTH Aachen University

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
genOpts.solvThreads     = 4;    % Number of parallel threads
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


%% Initialize gcOpt
results = init_gcOpt(model,probOpts,genOpts);


%% Save results to a ".mat" file
% save('results_example_gcOpt','results')

