function [KO,KD,Gene_KO,Gene_KD,Ins] = translatePop(model,model_c,indv,optFocus,targets,K,K_hri)
% Translate individuals of a indvulation to intervention strategies

%%
% INPUT
% model:                original metabolic model
% model_c:              compressed metabolic model
% indv:                 target set individual
% optFocus:             individuals carry reaction or gene targets
% targets:              struct containing all necessary about the target space
% K:                    Number of deletion targets
% K_hri:                Number of insertion targets

% OUTPUT
% KO:                   Deletions targets
% KD:                   Knockdown targets
% Gene_KO               gene deletion targets
% Gene_KD               gene knockdown targets
% Ins:                  heterologous inseertion targets


%%

% get parameters
indv_del = indv(1:K);     % deletion targets of individual
indv_ins = indv((K+1):(K_hri+K));   % insertion targets

% % check if compressed model contains transformation vector and matrices
% if ~isfield(model_c,'comprMapVec')
%     model_c.comprMapVec     = [1:length(model_c.rxns)]';
% end


Gene_KO     = [];
Gene_KD     = [];
if strcmp(optFocus,'gene')
    % get reaction targets from gene targets
    [targetRxnNum,~,targetRxnBounds,~]  = evalTargets(indv,targets);
    % extract deletion gene targets
    gene_targets        = {targets.genes{indv_del}};
    gene_targets_num    = cell(K,1); 
    for j=1:K
        % get gene numbers in original model (decompress)
        if gene_targets{j}(1)=='a'  % artifical gene (cluster of multiple genes)
            gene_targets_num{j}     = find(strcmp(model_c.geneMatCompress,gene_targets{j}));
            gene_targets_dc         = model.genes(gene_targets_num{j});
        else
            gene_targets_num{j}     = find(strcmp(model_c.genes,gene_targets{j}));
            gene_targets_dc         = model_c.genes(gene_targets_num{j});
        end
        gene_targets_num{j}     = find(ismember(model.genes,gene_targets_dc));
        % distinguish KD /KO targets
        bound_gene  = targets.bound_gene(find(strcmp(targets.genes,gene_targets{j})),:);
        if all(bound_gene==0)
            Gene_KO     = [Gene_KO,gene_targets_num{j}'];
        else
            Gene_KD     = [Gene_KD,gene_targets_num{j}'];
        end
    end
    % calculate number of reaction targets
    K_rxns  = length(targetRxnNum)-K_hri;
    % analyze insertion targets
    targetRxnNum_hri    = targets.rxnNum(indv_ins+targets.shift);
    Ins                 = model_c.comprMapVec(targetRxnNum_hri);
    
elseif strcmp(optFocus,'rxns')
    % reaction targets
    targetRxnNum    = targets.rxnNum(indv_del);
    targetRxnBounds = targets.bound(indv_del,:);
    % extract gene targets (all possible candidates)
    for j=1:K
        gene_targets_num    = find(any(model_c.rxnGeneMat(targetRxnNum(j),:),1));
        gene_targets_dc     = model.genes(gene_targets_num);
        gene_targets_num    = find(ismember(model.genes,gene_targets_dc));

        % distinguish KD /KO targets
        if all(targetRxnBounds(j,:)==0)
            Gene_KO     = [Gene_KO,gene_targets_num'];
        else
            Gene_KD     = [Gene_KD,gene_targets_num'];
        end
    end
    % declare number of reaction deletion targets
    K_rxns  = K;
    % analyze insertion targets
    targetRxnNum_hri    = targets.rxnNum(indv_ins);
    Ins                 = model_c.comprMapVec(targetRxnNum_hri);
end

KO  = zeros(1,K_rxns);
KD  = zeros(1,K_rxns);
for j=1:K_rxns
    if all(targetRxnBounds(j,:)==0)
        % knockout target
        KO(j)   = unique(model_c.comprMapVec(targetRxnNum(j)));
    else
        % knockdown target
        KD(j)   = unique(model_c.comprMapVec(targetRxnNum(j)));
    end
end
KO  = KO(KO~=0);
KD  = KD(KD~=0);


    


end