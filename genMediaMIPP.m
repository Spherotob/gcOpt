function [Av,Ay,b,numKORxns,mSize,nonKnockTot,nonKnockTot_s,optExRxns] = genMediaMIPP(nRxns,nRxns_s,...
                    S_s,ub_s,lb_s,nMets,notknockables,optExRxnsSplit,...
                    optExRxnsIrr,mapRxnsRed,B)


Av      = [];
Ay      = [];
Am      = [];
b       = [];
I       = eye(nRxns_s);

optExRxns   = [optExRxnsSplit;optExRxnsIrr];

% Match non-knockables and optional external fluxes to split network
[~,nonKnock]        = find(mapRxnsRed(notknockables,:));
[~,nonKnock_s]      = find(B(nonKnock,:));

nonKnock_s_orig     = nonKnock_s;
nonKnock_orig       = nonKnock;


nonKnock_s          = nonKnock_s(~ismember(nonKnock_s,optExRxns));


nonKnockTot_s   = [nonKnock_s;optExRxns];
nonKnockTot     = find(any(B(:,nonKnockTot_s),2));

% knockable_s     = 1:nRxns_s;
% knockable_s     = knockable_s(~ismember(knockable_s,nonKnockTot_s))';

% knockable       = 1:nRxns;
% knockable       = knockable(~ismember(knockable,nonKnockTot))';


% numKORxns   = nRxns-numNonKO-length(find(any(B(:,optExRxns))));     % y Size
numKORxns   = nRxns-length(nonKnockTot);

mSize       = length(optExRxns);                                    % m Size

% Special constraints (substrate uptake, ATP maintenance etc.) -> include i
% general constraints!

% Constraints non-knockable reactions
Av      = [Av;I(nonKnock_s,:);-I(nonKnock_s,:)];
Ay      = [Ay;zeros(2*size(nonKnock_s,1),numKORxns)];
Am      = [Am;zeros(2*size(nonKnock_s,1),mSize)];
b       = [b;ub_s(nonKnock_s);-lb_s(nonKnock_s)];

% Equality constraints S*v=0
Av      = [Av;S_s;-S_s];
Ay      = [Ay;zeros(nMets,numKORxns);zeros(nMets,numKORxns)];
Am      = [Am;zeros(nMets,mSize);zeros(nMets,mSize)];
b       = [b;zeros(nMets,1);zeros(nMets,1)];

% Add media composition booleans
Av      = [Av;I(optExRxns,:);-I(optExRxns,:)];
Ay      = [Ay;zeros(2*mSize,numKORxns)];
Am      = [Am;-diag(ub_s(optExRxns));diag(lb_s(optExRxns))];
b       = [b;zeros(2*mSize,1)];


% KO assignments + flux boundaries (implying positive and knockable fluxes)

Iknockable                  = I;
Iknockable(nonKnockTot_s,:) = [];
ubKnockable                 = ub_s;
lbKnockable                 = lb_s;
ubKnockable(nonKnockTot_s)  = [];
lbKnockable(nonKnockTot_s)  = [];
Bknockable                  = B;
Bknockable(:,nonKnockTot_s) = [];
% Bknockable                      = Bknockable(any(Bknockable,2),:);
Bknockable(nonKnockTot,:)   = [];

Ay1     = Bknockable*diag(ubKnockable);

Ay2     = Bknockable*diag(lbKnockable);

Av      = [Av;Iknockable;-Iknockable];          
Ay      = [Ay;-Ay1';Ay2'];
Am      = [Am;zeros(2*size(Iknockable,1),mSize)];
b       = [b; zeros(2*(nRxns_s-size(nonKnockTot_s,1)),1)];

% Av      = [Av;I(knockable_s,:);-I(knockable_s,:)];          
% Ay      = [Ay;-diag(ub_s(knockable_s));diag(lb_s(knockable_s))];
% Am      = [Am;zeros(2*size(knockable_s,1),mSize)];
% b       = [b; zeros(2*size(knockable_s,1),1)];


% merge knockable matrices part
Ay      = [Am,Ay];


end
