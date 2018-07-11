function [popFit,popFitObj,popFD,popObjVal] = evalFitness(pop,pop_Tbin,model,model_i,Np,targets,fitFun_type,opt_fitFun)
% calculate fitness of the population's chromosomes. 

%% INPUT
% pop:              Population containing intervention strategies as
%                   chromosomes
% pop_Tbin          
% chr_map_obj:      Container map holding fitness of all unique
%                   chromosomes ever created. Keys are strings of the
%                   binary intervention format
% model:            Stoichiometric metabolic model
% model_i:          Irreversible stoichiometric metabolic model
% Np:               Number of chromosomes in the population
% targets:        struct containing target reaction numbers, scores and
%                 bounds
% fitFun_type:      Type of fitness function
% opt_fitFun:       fitness function options

%% OUTPUT

%% General parameters


%% calculate fitness
popFit      = zeros(Np,1);
popFD       = cell(Np,1);
popObjVal   = popFD;
for i=1:Np
    % calculate fitness
    % choose fitness function
    switch fitFun_type
        case 0
            % MiMBl fitness function
            [~,act_targets_i,~,act_targetBounds_i] = evalTargets(pop(i,:),targets);
            
%             act_targets_i       = targets.rxnNum_i(pop(i,:));    % actual intervention strategy (reaction numbers)
%             act_targetBounds_i  = targets.bound_i(pop(i,:),:);     % respective target bounds

            % MiMBl and simple growth model
            [popFit(i),popFD{i}]   = fitFun_MiMBl(model,model_i,act_targets_i,act_targetBounds_i,opt_fitFun);

        case 1
            % multiobjective fitness function
            
             [act_targets,act_targets_i,act_targetBounds,act_targetBounds_i] = evalTargets(pop(i,:),targets);

%             % get targets
%             act_targets_i       = targets.rxnNum_i(pop(i,:));
%             act_targets         = targets.rxnNum(pop(i,:));
%             % get target bounds
%             act_targetBounds_i  = targets.bound_i(pop(i,:),:);
%             act_targetBounds    = targets.bound(pop(i,:),:);

            % multiobjective fitness function
            [popFit(i),popFD{i},popObjVal{i}]    = fitFun_multiObj(model,model_i,act_targets_i,act_targets,...
                                               act_targetBounds_i,act_targetBounds,opt_fitFun);  

        otherwise
            error('Unknown fitness function')
    end 
end

%% correct fitness to account for the minimization of interventions

popFitObj   = popFit;   % save non-transformed fitnesses
if opt_fitFun.minInt
    popFit      = transformFitness(popFit,pop_Tbin,opt_fitFun);
end

end