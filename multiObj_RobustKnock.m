function [fitVal,fluxDist,minTarget] = multiObj_RobustKnock(model,gurProb,opt_fitFun,maxGrowth)
% Evaluates RobustKnock solution for multiobjective fitness function

% INPUT
% model:        metabolic model
% gurobi:       gurobi optimization problem
% opt_fitFun:   fitness function options

% OUTPUT
% fitVal:       RobustKnock fitness function value
% fluxDist:     flux distribution of MiMBl solution

%%
gurProb.obj         = gurProb.obj_P;
gurProb.modelsense  = 'min';
% set biomass formation constraints
gurProb.lb(model.bmRxnNum)   = maxGrowth*0.99;
gurProb.ub(model.bmRxnNum)   = maxGrowth;
% solve LP
sol     = gurobi(gurProb,opt_fitFun.gurParams);
if strcmp(sol.status,'INFEASIBLE')
    % infeasible model
    fitVal      = -1;
    fluxDist    = -1;
    minTarget   = -1;
else
    fitVal      = sol.x(model.targetRxnNum)/opt_fitFun.maxRbstKnock;
    fluxDist    = sol.x;
    minTarget   = sol.x(model.targetRxnNum);
end


end