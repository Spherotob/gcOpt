function [model,optExRxns,conMets,conRxns,allMets] = metSpikeRoutine(model,baseMet,netDepth,depthFlag)
% outputs network around a metabolite and adds exchange reactions to the
% specified metabolites

baseMetNum      = find(strcmp(model.mets,baseMet));
S               = model.S;

conMets     = cell(netDepth,1);
conRxns     = cell(netDepth,1);

% metabolites to be disregarded
disMets     = {'H+';'H2O';'ATP';'Orthophosphate';'ADP';'NAD+';'NADH';...
                'Pyrophosphate';'CO2';'H+[e]';'NADP+';'NADPH';'NH4+';...
                'Oxygen';'CoA';'AMP';'FAD';'FADH2';'GMP';'UTP';'HCO3';'CTP';...
                'CDP';'CMP';'O2';'GTP';'GDP';'Phosphate'};
disMets     = find(ismember(model.metNames,disMets));           
disMets     = [disMets;find(~cellfun('isempty',strfind(model.metNames,'RNA')))];

jumpMets    = {'(E)-2-methylbutanal oxime';'(Z)-2-methylbutanal oxime';'(2R)-Hydroxy-2-methylbutanenitrile'};

%% Determine metabolites for each depth level
actLvlMets  = baseMetNum;
allConRxns  = [];
allMets     = [];
for l=1:netDepth
    conRxns{l}      = find(any(S(actLvlMets,:),1));
    conRxns{l}      = conRxns{l}(~ismember(conRxns{l},allConRxns));
    allConRxns      = [allConRxns,conRxns{l}];
    
    
    actLvlMetsAll   = find(any(S(:,conRxns{l}),2));
    actLvlMets      = actLvlMetsAll(~ismember(actLvlMetsAll,[actLvlMets;disMets]));
    conMets{l}      = actLvlMets;
    allMets         = unique([allMets;actLvlMets]);
end


%% Add exchange reactions
if depthFlag
    % only add exchange reaction for the highest depth
    mets     = conMets{netDepth};   
else
    % add exchange reactions for all metabolites up to the maximum depth
    mets     = allMets;
end

jumpMetsNum     = find(ismember(model.mets,jumpMets));
mets            = mets(~ismember(mets,jumpMetsNum));
numMets         = size(mets,1);
optExRxnsNum   = zeros(numMets,1);
for i=1:numMets
    % if metabolite is external choose existing exchange reaction
    exRxn           = find((model.S(mets(i),:)~=0).*(sum(model.S~=0,1)==1));
    if isempty(exRxn)
        model           = addExchangeRxn(model,model.mets(mets(i)),0,0);
        optExRxnsNum(i)    = size(model.rxns,1);
    else
        optExRxnsNum(i)    = exRxn;
    end
end

optExRxns   = model.rxns(optExRxnsNum);

end
