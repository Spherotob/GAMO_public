function [fitVal,fluxDist] = multiObj_MiMBl(model,model_i,opt_fitFun)
% Evaluates MiMBl solution for multiobjective fitness function

% INPUT
% model_i:      irreversible metabolic model
% opt_fitFun:   fitness function options

% OUTPUT
% fitVal:       MiMBl fitness function value
% fluxDist:     flux distribution of MiMBl solution


%%
% conduct MiMBl
sol     = MiMBl(model_i,model_i.fd_ref,1,opt_fitFun.excl_rxns_i);
if isempty(sol.x)
    % infeasible model
    fitVal      = -1;
    fluxDist    = -1;
else
    % calculate fitness
    fitVal     = (sol.x(model.bmRxnNum)*sol.x(model.targetRxnNum)...
                            /-sol.x(model.subsRxnNum))/opt_fitFun.maxMiMBl;
    fluxDist    = sol.x;
end

end