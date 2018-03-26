function [minFlux,maxFlux,Vmin,Vmax] = manualFVA(model,slctRxns,parOpt)

numRxns     = length(model.rxns);   % Total number of reactions 

if nargin<2
    slctRxns_idx     = 1:numRxns;   % Total number of reactions 
    slctRxns_num     = numRxns;
    parOpt           = 1;   % activate parallel computation
elseif nargin<3
    % analyze selected rxns only
    if ~isempty(slctRxns)
        slctRxns_idx     = zeros(length(slctRxns),1);
        for i=1:length(slctRxns)
            slctRxns_idx(i)  = find(ismember(model.rxns,slctRxns{i}));
        end
        slctRxns_num    = length(slctRxns_idx);
    else
        slctRxns_idx     = 1:numRxns;
        slctRxns_num     = numRxns;
    end
    parOpt           = 1;   % activate parallel computation
else
    % analyze selected rxns only
    if ~isempty(slctRxns)
        slctRxns_idx     = zeros(length(slctRxns),1);
        for i=1:length(slctRxns)
            slctRxns_idx(i)  = find(ismember(model.rxns,slctRxns{i}));
        end
        slctRxns_num    = length(slctRxns_idx);
    else
        slctRxns_idx     = 1:numRxns;
        slctRxns_num     = numRxns;
    end
end


%% Parallel Computing
v=ver;
PCT='Parallel Computing Toolbox';
if  any(strcmp(PCT,{v.Name}))% &&license('test',PCT)    
    p = gcp('nocreate');
    if isempty(p)
        poolsize = 0;
    else
        poolsize = p.NumWorkers;
    end
    PCT_status=1;
else
     PCT_status=0;  % Parallel Computing Toolbox not found.    
end


%% Loop through reactions
[minFlux,maxFlux]   = deal(zeros(slctRxns_num,1));
[Vmin,Vmax]         = deal(zeros(numRxns,slctRxns_num));     % Store flux distributions for each LP

if PCT_status && parOpt
    % Parallel computing


    parfor i=1:slctRxns_num
        % workaround for using COBRA optimization tool 
        globalSolverVariable('gurobi');
        
        rxnNum      = slctRxns_idx(i);
        passModel   = model;
        passModel   = changeObjective(passModel,passModel.rxns(rxnNum));
        % Maximize flux
        FBAsol      = optimizeCbModel(passModel,'max');
        if strcmp(FBAsol.origStat,'INFEASIBLE')
            maxFlux(i)  = 0;
            warning('INFEASIBLE MAX MODEL')
            Vmax(:,i)   = zeros(numRxns,1);
        else
            maxFlux(i)  = FBAsol.x(rxnNum);
            Vmax(:,i)   = FBAsol.x;
        end
        % Minimize flux
        FBAsol      = optimizeCbModel(passModel,'min');
        if strcmp(FBAsol.origStat,'INFEASIBLE')
            warning('INFEASIBLE MIN MODEL')
            minFlux(i)  = 0;
            Vmin(:,i)   = zeros(numRxns,1);
        else
            minFlux(i)  = FBAsol.x(rxnNum);
            Vmin(:,i)   = FBAsol.x;
        end
    end
      
else
    % Serial Computing
    for i=1:slctRxns_num
%         changeCobraSolver('gurobi6','LP');
        rxnNum      = slctRxns_idx(i);
        passModel   = model;
        passModel   = changeObjective(passModel,passModel.rxns(rxnNum));
        % Maximize flux
        FBAsol      = optimizeCbModel(passModel,'max');
        if strcmp(FBAsol.origStat,'INFEASIBLE')
            maxFlux(i)  = 0;
            Vmax(:,i)   = 0;
        else
            maxFlux(i)  = FBAsol.x(rxnNum);
            Vmax(:,i)   = FBAsol.x;
        end
        % Minimize flux
        FBAsol      = optimizeCbModel(passModel,'min');
        if strcmp(FBAsol.origStat,'INFEASIBLE')
            minFlux(i)  = 0;
            Vmin(:,i)   = 0;
        else
            minFlux(i)  = FBAsol.x(rxnNum);
            Vmin(:,i)   = FBAsol.x;
        end
    end
      
end


end