function solMiMBl = MiMBl(model,fluxDistWT,irrevFlag,excl_rxns)
% Conduction of Minimization of Metabolite Balance according to Brochado
% 2012

 
if nargin<4
    excl_rxns     = [];  
end

if irrevFlag
    modelIrrev          = model;
    irrevFluxDistWT     = fluxDistWT;
else
    modelIrrev        = rev2irr(model);
    irrevFluxDistWT   = fd_rev2irr(model,modelIrrev,fluxDistWT);
    excl_rxns         =[modelIrrev.rev2irrev{excl_rxns}]';
%     [modelIrrev,irrevFluxDistWT] = conv2Irr(model,fluxDistWT);
%     [~,excl_rxns]     = find(modelIrrev.mapIrr2Rev(excl_rxns,:)); 
end

% Parameter
[nMets,nRxns]   = size(modelIrrev.S);

% Strictly positive stoichiometric coefficients
alpha   = abs(modelIrrev.S);    

% Exclude specified fluxes from metabolite balance calculations
alpha(:,excl_rxns)    = 0;

% Metabolite turnover of wildtype flux distribution
tIrrevFluxDistWT    = alpha*irrevFluxDistWT;


%% First Optimization: Minimization of metabolite balances
% Consider absolute differences of wildtype and mutant metabolite balances
% Introduce help variables and rewrite LP:
% min sum(x)
%     s.t.
%     t_i,wt - t_i = x_i
%     and additional constraints for fluxes
         
fac     = Inf;
clearvars gurProb gurParams
% % No absolute objective function
% gurProb.rhs                     = [zeros(nMets,1);...
%                                     tIrrevFluxDistWT;...
%                                     -tIrrevFluxDistWT];
% gurProb.A                       = [sparse(nMets,nMets),sparse(modelIrrev.S);...
%                                     speye(nMets),sparse(alpha);...
%                                     speye(nMets),-sparse(alpha)];
% gurProb.obj                     = [ones(1,nMets),zeros(1,nRxns)];
% 
% gurProb.sense(1:nMets,1)        = '=';
% gurProb.sense((nMets+1):(3*nMets),1)    = '>';
% 
% gurProb.lb                      = [-fac.*ones(nMets,1);modelIrrev.lb];
% gurProb.ub                      = [fac.*ones(nMets,1);modelIrrev.ub];
% gurProb.vtype(1:(nMets+nRxns),1)    = 'C';

% Absolute objective funtion
gurProb.rhs                     = [zeros(nMets,1);...
                                    tIrrevFluxDistWT;...
                                    zeros(2*nMets,1)];
                                
gurProb.A                       = [sparse(nMets,nMets),sparse(modelIrrev.S),sparse(nMets,nMets);...
                                    speye(nMets),sparse(alpha),sparse(nMets,nMets);...
                                    speye(nMets),sparse(nMets,nRxns),-speye(nMets);...
                                    -speye(nMets),sparse(nMets,nRxns),-speye(nMets)];                             

gurProb.obj                     = [zeros(1,nMets),zeros(1,nRxns),ones(1,nMets)];


gurProb.sense(1:(2*nMets),1)                = '=';
gurProb.sense(((2*nMets)+1):(4*nMets),1)    = '<';

gurProb.lb                              = [-fac.*ones(nMets,1);modelIrrev.lb;zeros(nMets,1)];
gurProb.ub                              = [fac.*ones(nMets,1);modelIrrev.ub;fac.*ones(nMets,1)];
gurProb.vtype(1:((2*nMets)+nRxns),1)    = 'C';


gurProb.modelsense              = 'min';
% gurProb.objcon                  = sum(tmWT);


% gurobi parameters
gurParams.OutputFlag    = 0;
gurParams.Presolve      = 1;
gurParams.Threads       = 1;
gurParams.Method        = 1;    % use dual simplex for better convergense
% gurParams.Cutoff        = 0;


foSol       = gurobi(gurProb,gurParams);

% Cancel function if first optimization model is infeasible
if strcmp(foSol.status,'INFEASIBLE')
%     solMiMBl.x              = NaN.*ones(size(modelIrrev.mapIrr2Rev,1),1);
    solMiMBl.x              = [];
    solMiMBl.x_i            = [];
    solMiMBl.xWT            = modelIrrev.mapIrr2Rev*irrevFluxDistWT;
    solMiMBl.firstOptSol    = foSol;
    return;
else
%     trueobjval  = [ones(1,nMets),zeros(1,nMets+nRxns)]*foSol.x;     % objective value for non-absolute metabolite balances
    trueobjval  = foSol.objval;     % objective value for absolute metabolite balances
end
    
%% Second Optimization: Minimization of fluxes
% Here, consider optimal objective of previous LP
% As for the first optimization problem, absolute values have to be
% considered in the objective

% Only consider reactions with non-zero fluxes in the wildtype flux
% distribution

irrevFluxDistWT(irrevFluxDistWT<1e-7)   = 0;    % remove too small values
obj                     = ones(1,nRxns)./abs(irrevFluxDistWT');
obj(irrevFluxDistWT==0) = 0;

% exclude specified reactions from objective
obj(excl_rxns)        = 0;


clearvars gurProb 
% % No absolute objective function
% gurProb.obj                     = [obj,zeros(1,nRxns+nMets)];
% 
% gurProb.rhs                     = [zeros(nMets,1);...
%                                     foSol.objval;...
%                                     irrevFluxDistWT;...
%                                     -irrevFluxDistWT;...
%                                     tIrrevFluxDistWT;...
%                                     -tIrrevFluxDistWT];
%                                 
% gurProb.A                       = [sparse(nMets,nMets+nRxns),sparse(modelIrrev.S);...
%                                     sparse(1,nRxns),ones(1,nMets),sparse(1,nRxns);...
%                                     speye(nRxns),sparse(nRxns,nMets),speye(nRxns);...
%                                     speye(nRxns),sparse(nRxns,nMets),-speye(nRxns);...
%                                     sparse(nMets,nRxns),speye(nMets),sparse(alpha);...
%                                     sparse(nMets,nRxns),speye(nMets),-sparse(alpha)];
%                                 
% gurProb.sense(1:(nMets+1),1)                        = '=';
% gurProb.sense((nMets+2):((3*nMets)+(2*nRxns)+1),1)  = '>';
% gurProb.lb                      = [0.*ones(nMets+nRxns,1);modelIrrev.lb];
% gurProb.ub                      = [fac.*ones(nMets+nRxns,1);modelIrrev.ub];
% gurProb.vtype(1:(nMets+(2*nRxns)),1)    = 'C';


% % Absolute objective function values
% gurProb.obj                     = [zeros(1,(2*nRxns)+nMets),obj];
% 
% gurProb.rhs                     = [zeros(nMets,1);...
%                                     %foSol.objval;...
%                                     trueobjval;...
%                                     irrevFluxDistWT;...
%                                     tIrrevFluxDistWT;...
%                                     zeros(2*nRxns,1)];
%                                 
% gurProb.A                       = [sparse(nMets,nMets+nRxns),sparse(modelIrrev.S),sparse(nMets,nRxns);...
%                                     sparse(1,nRxns),ones(1,nMets),sparse(1,nRxns),sparse(1,nRxns);...
%                                     speye(nRxns),sparse(nRxns,nMets),speye(nRxns),sparse(nRxns,nRxns);...
%                                     sparse(nMets,nRxns),speye(nMets),sparse(alpha),sparse(nMets,nRxns);...
%                                     speye(nRxns),sparse(nRxns,nMets+nRxns),-speye(nRxns);...
%                                     -speye(nRxns),sparse(nRxns,nMets+nRxns),-speye(nRxns)];
%                             
%                                 
% gurProb.sense(1:((2*nMets)+1+nRxns),1)                          = '=';
% gurProb.sense(((2*nMets)+2+nRxns):((2*nMets)+(3*nRxns)+1),1)    = '<';
% 
% gurProb.lb                      = [-fac.*ones(nMets+nRxns,1);modelIrrev.lb;zeros(nRxns,1)];
% gurProb.ub                      = [fac.*ones(nMets+nRxns,1);modelIrrev.ub;fac.*ones(nRxns,1)];
% 
% gurProb.vtype(1:(nMets+(3*nRxns)),1)    = 'C';


% Absolute objective function values and absolute metabolite balances
gurProb.obj                     = [zeros(1,(2*nRxns)+nMets),obj,zeros(1,nMets)];

gurProb.rhs                     = [zeros(nMets,1);...
                                    %foSol.objval;...
                                    trueobjval;...
                                    irrevFluxDistWT;...
                                    tIrrevFluxDistWT;...
                                    zeros(2*nRxns,1);...
                                    zeros(2*nMets,1)];
                                
                                       % v_i_WT - v_i  , t_i_WT - t_i  ,  v_i  ,   |v_i_WT - v_i|
gurProb.A                       = [sparse(nMets,nMets+nRxns),sparse(modelIrrev.S),sparse(nMets,nRxns),sparse(nMets,nMets);...
                                    sparse(1,nRxns),sparse(1,nMets),sparse(1,nRxns),sparse(1,nRxns),ones(1,nMets);...
                                    speye(nRxns),sparse(nRxns,nMets),speye(nRxns),sparse(nRxns,nRxns),sparse(nRxns,nMets);...
                                    sparse(nMets,nRxns),speye(nMets),sparse(alpha),sparse(nMets,nRxns),sparse(nMets,nMets);...
                                    speye(nRxns),sparse(nRxns,nMets+nRxns),-speye(nRxns),sparse(nRxns,nMets);...
                                    -speye(nRxns),sparse(nRxns,nMets+nRxns),-speye(nRxns),sparse(nRxns,nMets);...
                                    sparse(nMets,nRxns),speye(nMets),sparse(nMets,2*nRxns),-speye(nMets);...
                                    sparse(nMets,nRxns),-speye(nMets),sparse(nMets,2*nRxns),-speye(nMets)];
%                             
                                
gurProb.sense(1:((2*nMets)+1+nRxns),1)                          = '=';
gurProb.sense(((2*nMets)+2+nRxns):((4*nMets)+(3*nRxns)+1),1)    = '<';

gurProb.lb                      = [-fac.*ones(nMets+nRxns,1);modelIrrev.lb;zeros(nRxns,1);zeros(nMets,1)];
gurProb.ub                      = [fac.*ones(nMets+nRxns,1);modelIrrev.ub;fac.*ones(nRxns,1);fac.*ones(nMets,1)];

gurProb.vtype(1:((2*nMets)+(3*nRxns)),1)    = 'C';


gurProb.modelsense              = 'min';

soSol       = gurobi(gurProb,gurParams);


solMiMBl.firstOptSol    = foSol;
solMiMBl.secondOptSol   = soSol;

% truesecondObj   = [zeros(1,nRxns),ones(1,nMets),zeros(1,2*nRxns)]*soSol.x;

% save irreversible solution


% Output flux solutions in reversible format
if strcmp(soSol.status,'OPTIMAL') 
    solMiMBl.x_i        = soSol.x(nMets+nRxns+1:(end-nRxns-nMets));
    solMiMBl.x          = modelIrrev.mapIrr2Rev*soSol.x(nMets+nRxns+1:(end-nRxns-nMets));
    solMiMBl.rc_so      = modelIrrev.mapIrr2Rev*soSol.rc(nMets+nRxns+1:(end-nRxns-nMets));
    solMiMBl.rc_fo      = modelIrrev.mapIrr2Rev*foSol.rc((nMets+1):(nMets+nRxns));
%     solMiMBl.x              = modelIrrev.mapIrr2Rev*soSol.x(nMets+nRxns+1:end);
else
%     solMiMBl.x              = NaN.*ones(size(modelIrrev.mapIrr2Rev,1),1);
    solMiMBl.x              = [];
    solMiMBl.x_i            = [];
end
    
solMiMBl.xWT            = modelIrrev.mapIrr2Rev*irrevFluxDistWT;

% save('foSol','foSol')
end