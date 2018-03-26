function [fitVal,fluxDist] = multiObj_gcOpt(model,gurProb,opt_fitFun,maxGrowth,minTarget)
% Evaluates gcOpt solution for multiobjective fitness function

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
growth_red                   = maxGrowth*0.9;
gurProb.lb(model.bmRxnNum)   = growth_red;
gurProb.ub(model.bmRxnNum)   = growth_red;
% solve LP
sol     = gurobi(gurProb,opt_fitFun.gurParams);
if strcmp(sol.status,'INFEASIBLE')
    % infeasible model
    fitVal      = -1;
    fluxDist    = -1;
else
    % calculate simplified GCS, is seen as the worst estimation
    % calculate area below the minimally guaranteed production yield
    b           = ((minTarget/maxGrowth)-(sol.x(model.targetRxnNum)/growth_red))/...
                        ((1/maxGrowth)-(1/growth_red));    % intercept
    a           = (minTarget-b)/maxGrowth;  % slope
    if b<=0
        % weak growth coupling
        fitVal  = (minTarget*(maxGrowth-(-b/a)))/2;
    else
        % strong growth coupling
        fitVal  = (maxGrowth*minTarget)-(((minTarget-b)*maxGrowth)/2);
    end
    % add area beyond the maximal growth rate
    fitVal  = fitVal+((minTarget*(opt_fitFun.maxMu-maxGrowth))/2);
    fitVal  = fitVal/opt_fitFun.maxgcOpt;
    
    fluxDist    = sol.x;
end


end