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
    % extract phenotype parameter
    mu      = sol.x(model.bmRxnNum);
    uptake  = sol.x(model.subsRxnNum);
    if mu<opt_fitFun.minGrowth || uptake==0
        fitVal  = 0;
    else
        % calculate fitness and choose fitness parameter
        switch opt_fitFun.fitParam
            case 0
                % biomass-product coupled yield
                fitVal     = (sol.x(model.bmRxnNum)*sol.x(model.targetRxnNum)...
                                        /-sol.x(model.subsRxnNum))/opt_fitFun.maxMiMBl;
            case 1
                % yield
                fitVal     = (sol.x(model.targetRxnNum)/-sol.x(model.subsRxnNum))/opt_fitFun.maxMiMBl;
            case 2
                % production rate
                fitVal     = sol.x(model.targetRxnNum)/opt_fitFun.maxMiMBl;                                
        end
    end
    fluxDist    = sol.x;
end

end