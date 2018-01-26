function [results] = solveGcOptMILP(model_s,A_j,b_j,lb_j,ub_j,vSize, zSize, lSize,mSize, ySize,B,B_flux,mapRxnsRed,nonKnock,results_wt,solver,genOpts)


fprintf('--> Solve Optimization Problem \n\n')

d_vSize     = zSize+lSize;
cO_s        = model_s.cO_s;

cO_j    = [cO_s',zeros(1,d_vSize+ySize+mSize)];                % Outer objective  
numDualConstr   = size(b_j,1);


A_j         = [-A_j(1,:);A_j];
b_j         = [0;b_j];


switch solver.MILP
    case 0      % gurobi solver
        gurModel.A          = A_j;
        gurModel.obj        = cO_j;
        gurModel.rhs        = full(b_j);
        gurModel.lb         = lb_j;
        gurModel.ub         = ub_j;
        
        % sense type
        sense               = cell(size(A_j,1),1);
        sense(:)            = {'<'};
        gurModel.sense      = char(sense);
        
        % Variable Type
        vtype                               = cell(d_vSize+vSize+mSize+ySize,1);
        vtype(:)                            = {'C'};  
        vtype((end-(ySize+mSize)+1):end,:)  = {'B'};    
        gurModel.vtype                      = char(vtype);
        
        gurModel.modelsense = 'max';
        
        % Optimisation parameters (See gurobi documentation)
        params.Cutoff           = 1e-4;                 % results_wt.optObjFlux/genOpts.Cutoff;
        params.ObjScale         = 1e-03;                % Objective Value is divided by the max objective coeeficient to the power of ObjScale
        params.FeasibilityTol   = 1e-07;                % Tolerance for constraint violations
        params.MIPGap           = 1e-03;                % Relative difference between incumbent and best bound below which solver terminates
        params.IntFeasTol       = 1e-09;
        params.Presolve         = genOpts.Presolve;     % Presolve
        params.Threads          = genOpts.solvThreads;  % Number of parallel threads
        params.TimeLimit        = genOpts.TimeLimit;
        params.Heuristics       = 0.01;
        params.Cuts             = genOpts.solvCuts;     % Defines aggressiveness of cutting plane routine
        params.BranchDir        = -1;                   % Searches lower branch first

        % Initialise gurobi optimisation
        res     = gurobi(gurModel,params);
        
        results.solverOutput    = res;
%         save optim
        if isfield(res,'x') 
            resKO       = res.x(d_vSize+vSize+mSize+1:end,1);   % Extract booleans related to reaction deletions
            resMed      = res.x(d_vSize+vSize+1:(d_vSize+vSize+mSize),1);
            resFluxes   = res.x(1:vSize,1);
        else
            resKO       = ones(ySize,1);
            resFluxes   = zeros(vSize,1);
            resMed      = zeros(mSize,1);
        end            
end

%% Map KOs
KOmap               = mapRxnsRed;
KOmap(:,nonKnock)   = [];
Kos                 = KOmap*resKO;
Kos(nonKnock,1)     = 1;
results.KOs         = Kos;
results.KORxnNames  = model_s.model.rxnNames(Kos<1e-5);
results.KORxnNum    = find(Kos<1e-05);
results.nonKnock    = nonKnock;
results.mapRxnsRed  = mapRxnsRed;

%% Map optional media compositions
try
    results.resMed          = resMed;
catch
    results.resMed          = [];
end
    
%% Map fluxes
% Use "B" matrix of the compressed model
results.fluxes      = mapRxnsRed*(B_flux*resFluxes);


%% Map objective
results.objRxn      = model_s.model.rxnNames(find(mapRxnsRed(:,find(B(:,find(cO_s)),1))));

% save resultsSolver
end