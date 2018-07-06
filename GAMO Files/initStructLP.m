function [gurProb,gurProb_1norm,gurParams] = initStructLP(model)
% initilization of structures needed to solve linear programm for maximizing biomass formation
% using gurobi (6 and later)

% INPUTS
% model:        Metabolic model COBRA format

% OUTPUTS
% gurProb:      gurobi problem structure
% gurProb_1norm gurobi problem structure of 1-norm minimization problem
% gurParams:    gurobi options


%% general parameter
[nMets,nRxns]   = size(model.S);


%% setup gurobi LP structure
% lefthand side of constraints
gurProb.A       = sparse(model.S);
% righthand side of linear constraints                
gurProb.rhs     = zeros(nMets,1);
% sense of constraints
gurProb.sense(1:nMets)  = '=';               
% objective function (biomass formation, set as default)               
gurProb.obj_BM      = +strcmp(model.rxns,model.bmRxn);
gurProb.obj         = gurProb.obj_BM;
% objective function (product formation)               
gurProb.obj_P       = +strcmp(model.rxns,model.targetRxn);
% sense of objective function (maximization)
gurProb.modelsense  = 'max';
% variable type
gurProb.vtype(1:nRxns)      = 'C';
% lower variable bounds
gurProb.lb      = model.lb;
% upper variable bounds
gurProb.ub      = model.ub;


%% setup gurobi LP for minimizing 1-norm 
% lefthand side of constraints
gurProb_1norm.A         = [model.S,sparse(nMets,nRxns);...
                            speye(nRxns),speye(nRxns);...
                            -speye(nRxns),speye(nRxns);...
                            gurProb.obj_BM',sparse(1,nRxns)];
% righthand side of linear constraints                
gurProb_1norm.rhs       = zeros(nMets+1+(2*nRxns),1);   % last value is the optimal objective function value of the corresponding maximization for biomass LP
% sense of constraints
gurProb_1norm.sense(1:nMets)                        = '='; 
gurProb_1norm.sense((nMets+1):(nMets+(2*nRxns)))    = '>'; 
gurProb_1norm.sense(nMets+1+(2*nRxns))              = '='; 
% objective function
gurProb_1norm.obj   = [zeros(nRxns,1);ones(nRxns,1)];
% sense of objective function (maximization)
gurProb_1norm.modelsense  = 'min';
% variable type
gurProb_1norm.vtype(1:(2*nRxns))    = 'C';
% lower variable bounds
gurProb_1norm.lb      = [model.lb;zeros(nRxns,1)];
% upper variable bounds
gurProb_1norm.ub      = [model.ub;ones(nRxns,1).*Inf];


%% gurobi Options
gurParams.OutputFlag    = 0;
gurParams.Presolve      = 1;
gurParams.Threads       = 1;

end