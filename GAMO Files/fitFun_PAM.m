function [fitVal,fluxDist] = fitFun_PAM(model,targets,targetBounds,opt_fitFun)
% fitness function employing protein allocation models. Flux distributions
% are calculated by simply maximizing the biomass formation reaction.

% INPUT
% model:        Reversible metabolic Model in COBRA format 
% targets:      Reaction numbers of targets for the respective intervention
%               strategy.
% targetBounds: Lower and upper bounds for specific targets
% opt_fitFun:   fitness function options (Here MiMBL options)

% OUTPUT

%% General Parameter


%% apply targets
% allocate mutant model
model_m     = model;
% lower bounds
model_m.lb(targets)      = targetBounds(:,1);
% upper bounds
model_m.ub(targets)      = targetBounds(:,2);


%% optimize model (determine mutant substrate uptake)
[fluxDist,~,~] = determineMutantPhenotype(model,model_m,opt_fitFun.subsGrid,opt_fitFun.initGridNum,...
                        opt_fitFun.maxOEELimit,opt_fitFun.wildtypeGridSolutions);


% check if valid solutions exists
if fluxDist==-1
    fitVal  = -1;
    return
end

% calculate fitness value from flux distribution
mu      = fluxDist(model.bmRxnNum);
uptake  = fluxDist(model.subsRxnNum);
if mu<opt_fitFun.minGrowth || uptake==0
    fitVal  = 0;
else
    % choose fitness function paramter
    switch opt_fitFun.fitParam
        case 0
            % Biomass-product coupled yield
            fitVal  = (fluxDist(model.targetRxnNum)*mu)/uptake;
        case 1
            % yield
            fitVal  = fluxDist(model.targetRxnNum)/uptake; 
        case 2
            % production rate
            fitVal  = fluxDist(model.targetRxnNum); 
    end
end

end