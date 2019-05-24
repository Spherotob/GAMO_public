function [popFit,popFitObj,popFD,popObjVal,chr_map] = evalFitness_mem(pop,pop_Tbin,chr_map,model,model_i,Np,targets,fitFun_type,opt_fitFun)
% calculate fitness of the population's chromosomes. Check if fitness has
% already been calculated and saved before.

%% INPUT
% pop:              Population containing intervention strategies as
%                   chromosomes
% pop_Tbin:         Population in binary intervention format
% chr_map:          Container map holding fitness of all unique
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
% popFitObj         Non-transformed fitnesses

%% General parameters

%% check if fitness already exist for the chromosomes

keySet      = cellstr(num2str(sort(pop,2))); 
isAvail     = isKey(chr_map,keySet);


%% calculate fitness
popFit      = zeros(Np,1);
popFD       = cell(Np,1);
popObjVal   = cell(Np,1);

for i=1:Np
    if isAvail(i)
        % load fitness and objective values from container
        keyValue        = chr_map(keySet{i});
        popFit(i)       = keyValue(1);    % First entry is the total fitness
    else
        % calculate fitness
            % choose fitness function
        switch fitFun_type
            case 0
                % MiMBl fitness function
                [~,act_targets_i,~,act_targetBounds_i] = evalTargets(pop(i,:),targets);

%                 act_targets_i       = targets.rxnNum_i(pop(i,:));    % actual intervention strategy (reaction numbers)
%                 act_targetBounds_i  = targets.bound_i(pop(i,:),:);     % respective target bounds

                % MiMBl and simple growth model
                [popFit(i),popFD{i}]   = fitFun_MiMBl(model,model_i,act_targets_i,act_targetBounds_i,opt_fitFun);

                % Update container map of fitness values
                chr_map(keySet{i})      = popFit(i);
                
            case 1
                % multiobjective fitness function
                [act_targets,act_targets_i,act_targetBounds,act_targetBounds_i] = evalTargets(pop(i,:),targets);

                
                
%                 % get targets
%                 act_targets_i       = targets.rxnNum_i(pop(i,:));
%                 act_targets         = targets.rxnNum(pop(i,:));
%                 % get target bounds
%                 act_targetBounds_i  = targets.bound_i(pop(i,:),:);
%                 act_targetBounds    = targets.bound(pop(i,:),:);

                % multiobjective fitness function
                [popFit(i),popFD{i},popObjVal{i}]    = fitFun_multiObj(model,model_i,act_targets_i,act_targets,...
                                                   act_targetBounds_i,act_targetBounds,opt_fitFun);  
             
                                               
                % Update container map of fitness values
                chr_map(keySet{i})      = [popFit(i),popObjVal{i}];                               
                                               
                  
            case 2
                % protein allocation model optimization
                [act_targets,act_targets_i,act_targetBounds,act_targetBounds_i] = evalTargets(pop(i,:),targets);
                % optimize model
                [popFit(i),popFD{i}]     = fitFun_PAM(model,act_targets,act_targetBounds,opt_fitFun);
            otherwise
                error('Unknown fitness function')
        end 
    end
end

%% correct fitness to account for the minimization of interventions

popFitObj   = popFit;   % save non-transformed fitnesses
if opt_fitFun.minInt
    popFit      = transformFitness(popFit,pop_Tbin,opt_fitFun);
end

end