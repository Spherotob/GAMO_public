function [fitVal,fluxDist,maxGrowth] = multiObj_OptKnock(model,gurProb,gurProb_1norm,opt_fitFun)
% Evaluates OptKnock solution for multiobjective fitness function

% INPUT
% model:        metabolic model
% gurobi:       gurobi optimization problem
% opt_fitFun:   fitness function options

% OUTPUT
% fitVal:       OptKnock fitness function value
% fluxDist:     flux distribution of MiMBl solution
% maxGrowth     Maximal growth rate

%%
gurProb.obj         = gurProb.obj_BM;
gurProb.modelsense  = 'max';
% solve LP, maximize for biomass
sol     = gurobi(gurProb,opt_fitFun.gurParams);
if strcmp(sol.status,'INFEASIBLE')
    % infeasible model
    fitVal      = -1;
    fluxDist    = -1;
    maxGrowth   = -1;
else
    % solve 1-norm minimization problem to minimize flux values
    gurProb_1norm.rhs(end)  = sol.objval;    % constrain to optimal objective function value of previous LP
    sol_1norm   = gurobi(gurProb_1norm,opt_fitFun.gurParams);
    if strcmp(sol_1norm.status,'INFEASIBLE')
        % infeasible model
        fitVal      = -1;
        fluxDist    = -1;
        maxGrowth   = -1;
    else
        fitVal      = sol_1norm.x(model.targetRxnNum)/opt_fitFun.maxOptKnock;
        fluxDist    = sol_1norm.x(1:opt_fitFun.nRxns);
        maxGrowth   = sol_1norm.x(model.bmRxnNum);
    end
end


end