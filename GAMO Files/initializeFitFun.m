function [popFit,popFitObj,chr_map,targets,model,model_i,opt_fitFun] = initializeFitFun(model,targets,pop,pop_Tbin,fitFun_type,opt_fitFun,optFocus)
% calculate fitness of initial population and build storage container map
% additionally transform reversible model an irreversible model if
% necessary

% INPUT
% model:          Metabolic model in COBRA format
% targets:        struct containing target reaction numbers, scores and
%                 bounds
% pop:            initial population in intervention number format
% fitFun_type:    Type of fitness function
% opt_fitFun:     Options for fitness function evaluation

% OUTPUT
% popFit:         Fitness of the initial population
% targets:        Struct containing all possible targets
% model:          Updated metabolic model in COBRA format
% model_i:        Irreversible model


%% General parameter
if ispc
    slash   = '\';
else
    slash   = '/';
end
[popSize,K_tot] = size(pop);
Nt_tot          = length(targets.rxnNum);   % number of possible deletion and insertion targets (reactions!)
% weight_objVal   = 1/7;  % weighting factor of slave objective values

% save maximal number of interventions
opt_fitFun.numIntMaxFit     = K_tot;

% translate optimization focus to boolean
% (0): gene targets
% (1): reaction targets
opt_fitFun.optFocus     = strcmp(optFocus,'rxns');




%% create irreversible model (if necessary for fitness function)
% MiMBl and simple growth model
%         [model_i,refFluxDist_i] = conv2Irr(model,model.refFluxDist);
model_i                 = rev2irr(model);
model_i.fd_ref          = fd_rev2irr(model,model_i,model.fd_ref);
model_i.subsRxnNum      = find(ismember(model_i.rxns,model_i.subsRxn));
model_i.hetRxnNum       = model_i.rev2irrev(targets.rxnNum_hri);
% adapt targets to irreversible model
% distinguish between genes and reaction targets
if opt_fitFun.optFocus
    % reaction targets
    targets.rxnNum_i        = model_i.rev2irrev(targets.rxnNum);
    bound_i                 = cell(Nt_tot,2);
    for i=1:length(targets.rxnNum)
        if length(targets.rxnNum_i{i}) > 1
            % reversible reaction
            if ~any(targets.bound(i,:))
                % reaction deletion
                bound_i{i,1}     = targets.bound(i,:);
                bound_i{i,2}     = targets.bound(i,:);
            elseif all(targets.bound(i,:)>=0)
                % knockdown or irreversible insertion target, positive flux values
                bound_i{i,1}     = [targets.bound(i,1),0];
                bound_i{i,2}     = [targets.bound(i,2),0];
            elseif all(targets.bound(i,:)<=0) 
                % knockdown or irreversible insertion target, negative flux values
                bound_i{i,1}     = [0,-targets.bound(i,2)];
                bound_i{i,2}     = [0,-targets.bound(i,1)];
            elseif (targets.bound(i,1)* targets.bound(i,2))<0
                % reversible insertion target
                bound_i{i,1}     = [0,0];
                bound_i{i,2}     = [targets.bound(i,2),-targets.bound(i,1)];
            end
        else
            % irreversible reaction
            bound_i{i,1}     = targets.bound(i,1);
            bound_i{i,2}     = targets.bound(i,2);
        end
    end
    targets.bound_i     = bound_i;
%     targets.rxnNum      = num2cell(targets.rxnNum);
else
  
    % gene targets
    B_t     = model_i.B(:,targets.rxnNum);  % mapping matrix for reaction targets only
    
    % determine direction of reference flux (mapIrr2Rev)
    B_t_bound   = B_t;
    for i=1:Nt_tot
        back_rxn    = model_i.mapIrr2Rev(targets.rxnNum(i),:)==-1;  % position of backwards flux
        forw_rxn    = model_i.mapIrr2Rev(targets.rxnNum(i),:)==1;  % position of forward flux
        if any(back_rxn) 
            % reversible reaction
            if model.fd_ref(targets.rxnNum(i))<0
                B_t_bound(back_rxn,i)     = 1;
                B_t_bound(forw_rxn,i)     = 0;
            else
                B_t_bound(back_rxn,i)     = 0;
                B_t_bound(forw_rxn,i)     = 1;
            end
        else
            % irreversible reaction
        end
    end
    
    % determine irreversible bounds for insertion targets
    rxnNum_i                = model_i.rev2irrev(targets.rxnNum);
    bound_i                 = cell(Nt_tot,2);
    for i=1:Nt_tot
        rxnPos              = i;
        if length(rxnNum_i{rxnPos}) > 1
            % reversible reaction
            if ~any(targets.bound(rxnPos,:))
                % reaction deletion
                bound_i{i,1}     = targets.bound(rxnPos,:);
                bound_i{i,2}     = targets.bound(rxnPos,:);
            elseif all(targets.bound(rxnPos,:)>=0)
                % knockdown or irreversible insertion target, positive flux values
                bound_i{i,1}     = [targets.bound(rxnPos,1),0];
                bound_i{i,2}     = [targets.bound(rxnPos,2),0];
            elseif all(targets.bound(rxnPos,:)<=0) 
                % knockdown or irreversible insertion target, negative flux values
                bound_i{i,1}     = [0,-targets.bound(rxnPos,2)];
                bound_i{i,2}     = [0,-targets.bound(rxnPos,1)];
            elseif (targets.bound(rxnPos,1)* targets.bound(rxnPos,2))<0
                % reversible insertion target
                bound_i{i,1}     = [0,0];
                bound_i{i,2}     = [targets.bound(rxnPos,2),-targets.bound(rxnPos,1)];
            end
        else
            % irreversible reaction
            bound_i{i,1}     = targets.bound(rxnPos,1);
            bound_i{i,2}     = targets.bound(rxnPos,2);
        end
    end
    targets.bound_i     = bound_i;
    targets.rxnNum_i    = rxnNum_i;
    
    
    targets.B_t         = B_t;
    targets.B_t_bound   = B_t_bound;
    targets.rev2irrev   = model_i.rev2irrev;
    

%     % create dependence matrix of gene bounds and irreversible reactions
%     % upper bound
%     nR_i    = length(model_i.rxns);     % number of irreversible reactions
%     
% 
%     targets.rxnNum_i    = cell(Nt_tot,1);
%     bound_i             = cell(Nt_tot,2);
%     for i=1:Nt_tot
%         rxnNum_i    = model_i.rev2irrev(targets.rxnNum(i));
%         numRxns     = length(rxnNum_i);
%         for j=1:numRxns
%             bound   = [targets.bound{i,1}(j),targets.bound{i,2}(j)];    % load bound of the specific target reaction
%             if length(rxnNum_i{j}) > 1
%                 % reversible reaction
%                 if ~any(bound)
%                     bound_i{i,1}        = [bound_i{i,1},bound];
%                     bound_i{i,2}        = [bound_i{i,2},bound];
%                 elseif all(bound>=0)
%                     % knockdown, positive flux
%                     bound_i{i,1}        = [bound_i{i,1},bound(1),0];
%                     bound_i{i,2}        = [bound_i{i,2},bound(2),0];
%                 elseif all(bound<=0)
%                     % knockdown, negative flux
%                     bound_i{i,1}        = [bound_i{i,1},0,-bound(1)];
%                     bound_i{i,2}        = [bound_i{i,2},0,-bound(2)];
%                 else
%                     warning('Target bounds suggest a reversible reaction')
%                 end
%             else
%                 % irreversible reaction
%                 bound_i{i,1}        = [bound_i{i,1},bound(1)];
%                 bound_i{i,2}        = [bound_i{i,2},bound(2)];
%             end
%         end
%         targets.rxnNum_i{i}     = [rxnNum_i{:}];
%     end 
%     targets.bound_i     = bound_i;
end


%% prepare model
% set bounds of herterologous reaction to zero
model.lb(model.hetRxnNum)     = 0;
model.ub(model.hetRxnNum)     = 0;
for i=1:length(model_i.hetRxnNum)
    model_i.lb(model_i.hetRxnNum{i})  = 0;
    model_i.ub(model_i.hetRxnNum{i})  = 0;
end      


%% special options for particular fitness function
switch fitFun_type
    case 0
        % MiMBl
        [~,opt_fitFun.excl_rxns_i]  = find(model_i.mapIrr2Rev(opt_fitFun.excl_rxns,:));
     
    case 1
        % Multiobjective fitness function        
        
        % determine active and leading objectives
        % (G):  gcOpt
        % (R):  RobustKnock
        % (O):  Optknock    
        % (M):  MiMBl
        switch opt_fitFun.leadObj
            case 'G'
                opt_fitFun.leadObj  = 4;
            case 'R'
                opt_fitFun.leadObj  = 3;
            case 'O'
                opt_fitFun.leadObj  = 2;
            case 'M'
                opt_fitFun.leadObj  = 1;
            otherwise
                error('No valid master objective specified!')
        end
        
        isObj                      = zeros(4,1);
        isObj(opt_fitFun.leadObj)  = 1;
        if ~isempty(opt_fitFun.isObj)
            for i=1:length(opt_fitFun.isObj)
                switch opt_fitFun.isObj(i)
                    case 'G'
                        isObj(4)  = 1;
                    case 'R'
                        isObj(3)  = 1;
                    case 'O'
                        isObj(2)  = 1;
                    case 'M'
                        isObj(1)  = 1;
                    otherwise
                        error('No valid objective specified!')
                end
            end  
        end
        % switch notation
        opt_fitFun.isObj   = isObj;
        
        % MiMBl
        [~,opt_fitFun.excl_rxns_i]  = find(model_i.mapIrr2Rev(opt_fitFun.excl_rxns,:));
        
        % write tailor-made fitFun_multiObj m-file
        fileID  = fopen(['AddFiles',slash,'fitFun_multiObj.m'],'w');
        % fix header
        fprintf(fileID,['function [fitVal,fluxDist,objVal] = fitFun_multiObj(model,model_i,targets_i,targets,targetBounds_i,targetBounds,opt_fitFun)\n',...
                        '%% General parameters\n',...   
                        'objVal = -ones(1,4);\n',...
                        'fluxDist = cell(1,4);\n',...
                         '%% Apply targets\n',...
                         'model_i.lb(targets_i) = targetBounds_i(:,1);\n',...
                         'model_i.ub(targets_i) = targetBounds_i(:,2);\n',...
                         'model_i.ub(model_i.subsRxnNum) = 1000;\n',...
                         '%% FBA problem formulation\n',...
                         'gurProb = opt_fitFun.gurProb;\n',... 
                         'gurProb.lb(targets) = targetBounds(:,1);\n',...
                         'gurProb.ub(targets) = targetBounds(:,2);\n',...
                         'gurProb_1norm = opt_fitFun.gurProb_1norm;\n',...
                         'gurProb_1norm.lb(targets) = targetBounds(:,1);\n',...
                         'gurProb_1norm.ub(targets) = targetBounds(:,2);\n',...
                         '%% Determine objectives\n']);
        
        % variable objective evaluation
        % MiMBl ?
        if isObj(1)
            fprintf(fileID,'[objVal(1),fluxDist{1}] = multiObj_MiMBl(model,model_i,opt_fitFun);\n');
        end
        % FBA problems?
        if any(isObj(2:end))
            % Optknock solution is always required
            fprintf(fileID,'[objVal(2),fluxDist{2},maxGrowth] = multiObj_OptKnock(model,gurProb,gurProb_1norm,opt_fitFun);\n');
        end
        if isObj(3) || isObj(4)
            % RobustKnock or gcOpt
            fprintf(fileID,['if maxGrowth>0\n\t',...
                '[objVal(3),fluxDist{3},minTarget] = multiObj_RobustKnock(model,gurProb,opt_fitFun,maxGrowth);\n']);
            if isObj(4)
                % also include gcOpt
                fprintf(fileID,['\tif minTarget>0\n\t\t',...
                    '[objVal(4),fluxDist{4}] = multiObj_gcOpt(model,gurProb,opt_fitFun,maxGrowth,minTarget);\n\t',...
                    'end\n end\n']);
            else
                fprintf(fileID,'end\n');
            end
        end
        % Caculate total fitness
        eval_Obj                        = isObj;
        eval_Obj(opt_fitFun.leadObj)    = 0;
        eval_Obj                        = find(eval_Obj);
        
        fprintf(fileID,['%% Calculate total fitness\n',...
                        'objVal_t = objVal;\n',...
                        'objVal_t(objVal==-1 | isnan(objVal)) = 0;\n',...
                        'fitVal = ((objVal_t(',num2str(opt_fitFun.leadObj),')',...
                         '*',num2str(opt_fitFun.weighting(opt_fitFun.leadObj)),')']);
        if isempty(eval_Obj)
            fprintf(fileID,[')*(objVal_t(',num2str(opt_fitFun.leadObj),')~=0);\n']); % if lead objective is zero then the total fitness is zero
        else
            for i=1:length(eval_Obj)
                fprintf(fileID,['+(objVal_t(',num2str(eval_Obj(i)),')']);
                fprintf(fileID,['*',num2str(opt_fitFun.weighting(eval_Obj(i))),')']);
            end
%             fprintf(fileID,';\n');
            fprintf(fileID,[')*(objVal_t(',num2str(opt_fitFun.leadObj),')~=0);\n']); % if lead objective is zero then the total fitness is zero
        end
        fprintf(fileID,[ 'fluxDist = fluxDist{',num2str(opt_fitFun.leadObj),'};\n','end']);
            
        fclose(fileID);            
        
         % set up standard LP
        [opt_fitFun.gurProb,opt_fitFun.gurProb_1norm,opt_fitFun.gurParams]   = initStructLP(model);
        opt_fitFun.nRxns    = length(model.rxns);   % save number of reactions (for 1-norm minimization)
       
        % calculate maximal theoretical production rate and maximal
        % objectives as normalization factors to be able to compare
        % different objectives in a multi-objective approach
        
        if sum(isObj)>1 || opt_fitFun.leadObj~=1
            % more 
            % maximal growth rate
            gurProb     = opt_fitFun.gurProb;
            gurProb.obj = gurProb.obj_BM;
            sol         = gurobi(gurProb,opt_fitFun.gurParams);
            maxMu       = sol.x(model.bmRxnNum);
            opt_fitFun.maxMu    = maxMu;
            % theoretical production rate
            gurProb.obj = gurProb.obj_P;
            gurProb.lb(model.bmRxnNum)  = 0;
            gurProb.ub(model.bmRxnNum)  = 0;
            sol                     = gurobi(gurProb,opt_fitFun.gurParams);
            opt_fitFun.theoProd     = sol.x(model.targetRxnNum);
            % assign optimal objectives
            opt_fitFun.maxOptKnock  = opt_fitFun.theoProd;
            opt_fitFun.maxRbstKnock = opt_fitFun.theoProd;
            opt_fitFun.maxgcOpt     = (opt_fitFun.theoProd*maxMu)/2;
            % optimal objective MiMBl, distinguish fitness parameter
            switch opt_fitFun.fitParam
                case 0
                    % Biomass-product coupled yield
                    grid        = 50;
                    step        = opt_fitFun.theoProd/grid;
                    prodRate    = 0:step:opt_fitFun.theoProd;
                    objMiMBl    = zeros(1+grid,1);
                    parfor i=1:(1+grid)
                        % load problem
                        gurProb         = opt_fitFun.gurProb;
                        gurProb.obj     = gurProb.obj_BM;
                        % set production rate
                        gurProb.lb(model.targetRxnNum)  = prodRate(i);
                        gurProb.ub(model.targetRxnNum)  = prodRate(i);
                        sol                             = gurobi(gurProb,opt_fitFun.gurParams);
                        % save objective
                        objMiMBl(i)     = prodRate(i)*sol.x(model.bmRxnNum)/-sol.x(model.subsRxnNum);
                    end
                    opt_fitFun.maxMiMBl = max(objMiMBl);
                    
                case 1
                    % yield
                    opt_fitFun.maxMiMBl     = opt_fitFun.theoProd/-sol.x(model.subsRxnNum);
                    
                case 2
                    % production rate
                    opt_fitFun.maxMiMBl     = opt_fitFun.theoProd;
            end
        else
            % only one objective function, dismiss normalization of
            % objective function value/fitness
            opt_fitFun.maxMiMBl     = 1;
            opt_fitFun.maxOptKnock  = 1;
            opt_fitFun.maxRbstKnock = 1;
            opt_fitFun.maxgcOpt     = 1;
        end
        
        % consider minimal allowable growth rate
        opt_fitFun.gurProb.lb(model.bmRxnNum)  = opt_fitFun.minGrowth;
        
        
    otherwise
        error('Unknown fitness function')   
        
end


%% Calculate Fitness of each chromosome
popFit      = zeros(popSize,1);     % fitness for each chromosome
popFD       = cell(popSize,1);
popObjVal   = popFD;
parfor i=1:popSize
    % choose fitness function
    switch fitFun_type
        case 0
%             act_targets_i       = targets.rxnNum_i(pop(i,:));    % actual intervention strategy (reaction numbers)
%             act_targetBounds_i  = targets.bound_i(pop(i,:),:);     % respective target bounds
              
            [~,act_targets_i,~,act_targetBounds_i] = evalTargets(pop(i,:),targets);     % m-file compatible
%             [~,act_targets_i,~,act_targetBounds_i] = evalTargets(pop(i,:),targets.bound,targets.bound_i,...
%                                                         targets.rxnNum,targets.rxnNum_i);   % mex file compatible

            % MiMBl and simple growth model
            [popFit(i),popFD{i}]   = fitFun_MiMBl(model,model_i,act_targets_i,act_targetBounds_i,opt_fitFun);
              
        case 1
 
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



%% initialize container for chromosomes
% use string format of pop
keySet      = cellstr(num2str(sort(pop,2)));

% use string of the binary intervention format of the population as keys
% keySet      = cellstr(num2str(pop_Tbin));
if ~isempty(popFit)
    switch fitFun_type
        case 0
            chr_map     = containers.Map(keySet,popFit);
        case 1
            chr_map     = containers.Map();
            for i=1:length(keySet)
                chr_map(keySet{i})      = [popFit(i),popObjVal{i}];
            end
    end
else
    chr_map     = containers.Map;
end


end