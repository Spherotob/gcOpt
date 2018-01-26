%% gcOpt application Example

% gcOpt requires gurobi MILP solver, COBRA Toolbox as well as CellNetAnalyzer. Make sure that
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

probOpts.bmRxn              = 'BIOMASS_Ecoli_core_w_GAM';
probOpts.subsRxn            = 'EX_glc__D_e';
probOpts.targetRxn          = 'EX_lac__D_e';


%% Additional or optional properties/options regarding the problem formulation


% Optimization direction
% 'max': maximization of objective reaction
% 'min': minimization of objective reaction
probOpts.sense  = 'max';

% Maximal allowable reaction deletions
probOpts.maxKO      = 2;

% Manually fix growth rate at which gcOpt optimizes the minimal guaranteed
% production rate (leave empty if growth rate should be automatically fixed
% at maximal productivity)
probOpts.fixMu      = 0.15;


%% Test for alternative, growth-coupling inducing medium compositions and hypothetical metabolite supply

% Alternative, growth-coupling inducing medium compositions
probOpts.medCheck   = 0;    % Set to "1" to activate function
% Specifiy optional exchange reactions
% probOpts.optExRxns  = model.rxns(find(not(cellfun('isempty',(strfind(model.rxns,'Ex_'))))));
probOpts.optExRxns      = model.rxns([71,76,78]);

% Hypothetical metabolite spiking to identify bottlenecks in metabolite supply 
probOpts.metSpikeCheck  = 0;            % Set to "1" to activate function
probOpts.baseMet        = '2-Butanone'; % Specifiy metabolite identifier of base metabolite
probOpts.depthFlag      = 1;            % 1: Only consider metabolites of the specified distance from the base metabolite
                                        % 0: Consider all metabolites within the specified distance from the base metabolite
probOpts.maxMedInt      = 1;            % Maximal number of parallel metabolite spiking reactions 

% Maximal metabolite or medium composition spiking rate [mmol/g/h]
probOpts.metSpikeRate   = 1;                                         
     

%% General options (Mostly regarding the gurobi solver)
% "genOpts" are optional, default values are assigned here
genOpts.solvThreads     = 3;    % Number of parallel threads
genOpts.solvCuts        = 1;    % MILP Cuts setting
genOpts.avBnd           = 500;  % Maximal upper bound of auxiliary variables (dualisation)
genOpts.TimeLimit       = Inf;  % Maximal runtime of the solver [sec]
genOpts.Presolve        = 1;    % Presolver of gurobi solver

% Network compression
genOpts.reductionFlag   = 1;
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

