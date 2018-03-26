function [fitVal,fluxDist] = fitFun_MiMBl(model,model_i,targets_i,targetBounds_i,opt_fitFun)
% fitness function using MiMBl and a simple growth model to evaluate the
% fitness of a chromosome (intervention strategy)

% INPUT
% model:        Reversible metabolic Model in COBRA format 
% model_i:      Irreversible metabolic Model in COBRA format containg reference flux
%               distribution (refFluxDist).
% targets:      Reaction numbers of targets for the respective intervention
%               strategy.
% targetBounds: Lower and upper bounds for specific targets
% opt_fitFun:   fitness function options (Here MiMBL options)

% OUTPUT

%% General Parameter
% X0  = 0.1;     % initial biomass


%% apply targets (irreversible)
% lower bounds
model_i.lb(targets_i)      = targetBounds_i(:,1);
% upper bounds
model_i.ub(targets_i)      = targetBounds_i(:,2);
% unconstrain substrate uptake 
model_i.ub(model_i.subsRxnNum)  = 1000;



% conduct MiMBl
sol     = MiMBl(model_i,model_i.fd_ref,1,opt_fitFun.excl_rxns_i);
if isempty(sol.x)
    % infeasible model
    fitVal      = -1;
    fluxDist    = -1;
    return
end
fluxDist    = sol.x;
% calculate productivity [mmol] (from timepoint 0 to 1) using a simple growth model
mu      = sol.x(model.bmRxnNum);
uptake  = sol.x(model.subsRxnNum);
if mu<=0 || uptake==0
    fitVal  = 0;
else
    fitVal  = (sol.x(model.targetRxnNum)*mu)/-uptake; 
end

end