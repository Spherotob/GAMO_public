function model_u = compressModel(model)
% compresses a metabolic model in terms of deleting blocked reactions
% INPUT
% model:    Stoichiometric metabolic model

% OUTPUT
% model_u:    Compressed model

c = fix(clock);
disp([num2str(c(4:end)),': Compress model ...'])

%% General variables
numRxns     = length(model.rxns);
numMets     = length(model.mets);

% %% create reference flux distribution
% model       = changeObjective(model,model.bmRxn);
% sol_ref     = optimizeCbModel(model,'max','one');
% fd_ref      = sol_ref.x;


%% rational analysis (see Gagneur 2004)
% conduct in an iterative manner
model_u         = model;    % reduced reversible model
keep_rxns_s     = eye(length(model.rxns));   % save reactions to be kept
redo            = 1;
while redo
    model_i         = rev2irr(model_u);
    % analyze model and get parameters
    S   = full(model_i.S);
    
%     % null space analysis of stoichiometric matrix
%     K           = null(S,'r');  % calculate kernel (null-space matrix) of S
%     null_rxns   = find(~any(K,2));  % exclude reactions with null flux at any steady state

    
    % setup gurobi problem of mass balance LP
    gurProb         = [];
    gurParams       = [];
    [nMets,nRxns]   = size(model_u.S);
    gurProb.A       = sparse(model_u.S);
    gurProb.rhs     = zeros(nMets,1);
    gurProb.sense(1:nMets)  = '=';
    gurProb.modelsense  = 'max';
    gurProb.vType(1:nRxns)  = 'C';
    % "unconstrain" flux variables !!!
    gurProb.lb      = ones(nRxns,1)*-Inf;
    gurProb.ub      = ones(nRxns,1)*Inf;
    % gurobi Options
    gurParams.OutputFlag    = 0;
    gurParams.Presolve      = 1;
    gurParams.Threads       = 1;
    gurParams.method        = 1;
    % calculate maximal and minimal values
    maxFlux     = zeros(nRxns,1);
    minFlux     = zeros(nRxns,1);
    parfor i=1:nRxns
        parGurProb          = gurProb;
        parGurProb.obj     = zeros(nRxns,1);
        parGurProb.obj(i)  = 1;
        % maximzation
        parGurProb.modelsense  = 'max';
        sol                 = gurobi(parGurProb,gurParams);
        if strcmp(sol.status,'OPTIMAL')
            maxFlux(i)  = sol.objval;
        else
            maxFlux(i)  = -1;
        end
        % minimization
        parGurProb.modelsense  = 'min';
        sol                 = gurobi(parGurProb,gurParams);
        if strcmp(sol.status,'OPTIMAL')
            minFlux(i)  = sol.objval;
        else
            minFlux(i)  = -1;
        end
    end
    null_rxns    = find(maxFlux==0 & minFlux==0);  

    % uni-directional reactions pointing to a sink/ from a source and
    % metabolites with only one reaction
    % determine sink and source metabolites
%     mult_b      = model_u.lb.*model_u.ub;   % multiply reaction bounds
%     mult_b      = mult_b';
%     en_rxns     = [model_u.lb~=0 | model_u.ub~=0]';     % enabled reactions
%     rxnsPerMet  = sum((S.*en_rxns)~=0,2);  % number of reactions incorporating a metabolite 
%     bi_rxns     = mult_b<0;     % determine bi-directional reactions
%     bi_S        = S.*bi_rxns;   % bi-directional reactions in S
%     bi_S        = bi_S~=0;
%     dir_S       = S.*mult_b;    % determine reaction directionalities in S (positive values imply uni-directionality)
%     ss_mets     = ((any(dir_S>=0,2) | any(dir_S<=0,2))  & any(dir_S,2) & ~any(bi_S,2)) | rxnsPerMet<2;    % determine sinks and sources
%     ss_rxns     = find(any(S(ss_mets,:),1));
    
    en_rxns     = [model_i.lb~=0 | model_i.ub~=0]';     % enabled reactions
    rxnsPerMet  = sum((S.*en_rxns)~=0,2);
    ss_mets     = ((all(S>=0,2) | all(S<=0,2)) & any(S,2)) | rxnsPerMet<2;    % determine sinks and sources
    ss_rxns_i   = find(any(S(ss_mets,:),1));
    ss_rxns     = [];
    for i=1:length(model_i.rev2irrev)
        if all(ismember(model_i.rev2irrev{i},ss_rxns_i))
            ss_rxns    = [ss_rxns;i];
        end  
    end

    
    % something
    
    % % reactions of constant flux ratios (row of S differing only by a scalar)
    % ident_rxns_s1    = {};
    % ident_rxns_s2    = [];
    % for r1=1:(numRxns-1)
    %     col_S           = S(:,r1);  % extract column from S
    %     entry_col_S     = find(col_S);  % find metabolites participating in reaction
    %     norm_S          = S./col_S(entry_col_S(1));     % normalize S by one stoichiometric coefficient
    %     norm_S_ident    = [ones(1,r1),sum(abs(norm_S(:,r1+1)-norm_S(:,r1)),1)];   % determing reaction with the same stoichiometrics
    %     ident_rxns      = find(norm_S_ident==0);
    %     if ~isempty(ident_rxns)
    %         ident_rxns_s1    = [ident_rxns_s1;ident_rxns];
    %         ident_rxns_s2    = [ident_rxns_s2;r1];        
    %     end  
    % end
    
    
    % delete reactions 
    del_rxns    = unique([null_rxns;ss_rxns]);
%     del_rxns_i  = [model_i.rev2irrev{del_rxns}];
    
%     del_rxns    = unique([ss_rxns]);
    if isempty(del_rxns)
        redo    = 0;
    else
%         model_i                 = deleteReactions(model_i,del_rxns);
        model_u                 = deleteReactions(model_u,del_rxns);
        keep_rxns_s(:,del_rxns) = [];
    end
    
    
end
%%
% % Analyze reactions to be deleted
keep_rxns   = find(any(keep_rxns_s,2));
blckd_rxns  = find(~any(keep_rxns_s,2));
% del_rxns    = [];
% for i=1:length(model.rxns)
%     if all(ismember(model_u.rev2irrev{i},del_rxns_i))
%         del_rxns    = [del_rxns;i];
%     end  
% end

%% create mapping matrix and vector
numBlck_rxns  = length(blckd_rxns);
comprMapMat   = zeros(numRxns,numBlck_rxns);
comprMapVec   = zeros(numBlck_rxns,1);
for j=1:numRxns
    rxnMap  = find(strcmp(model_u.rxns,model.rxns{j}));
    if ~isempty(rxnMap)
        comprMapMat(j,rxnMap)     = 1;
        comprMapVec(rxnMap)       = j;
    end
end
model_u.comprMapMat     = sparse(comprMapMat);
model_u.comprMapVec     = comprMapVec;


%% check consistency of compressed model



%% merge genes
geneMatCompress     = cell(length(model_u.genes),1);
genes_save          = model_u.genes;
genes_merged        = cell(length(model_u.genes),1);
for i=1:length(model_u.rxns)
    % check if multiple genes exclusively encode for the same reaction
    gene_inc            = find(model_u.rxnGeneMat(i,:));
    is_exclusive_bin    = sum(model_u.rxnGeneMat(:,gene_inc),1)==1;
    if sum(is_exclusive_bin)>1
        
        is_exclusive    = find(is_exclusive_bin);
        % create unifying gene name
        conti   = 1;
        while conti
            nums            = round(rand(4,1).*10);
            uni_gene        = 'a';
            for j=1:4
                uni_gene        = [uni_gene,num2str(nums(j))];
            end
            if ~any(strcmp(model_u.genes,uni_gene)) && length(uni_gene)==5
                conti   = 0;
            end
        end
%         % save which genes are incorporated in unifying gene
%         gene_inc_tot    = model_u.genes{gene_inc(is_exclusive(1))};
%         for j=1:length(is_exclusive)
%             gene_inc_tot    = [gene_inc_tot,'+',model_u.genes{gene_inc(is_exclusive(j))}];
%         end

        % rewrite gene rules            
        for j=1:length(is_exclusive)
            gene        = model_u.genes{gene_inc(is_exclusive(j))};
            lenGene     = length(gene);
            is_gene     = strfind(model_u.grRules,gene);
            pos_gene    = find(~cellfun('isempty',is_gene));
            for k=1:length(pos_gene)
                pos_ident   = is_gene{pos_gene(k)};
                grRule      = model_u.grRules{pos_gene(k)};
                lenRule     = length(grRule);
                for l=1:length(pos_ident)
                    pos     = strfind(model_u.grRules{pos_gene(k)},gene);
                    grRule  = model_u.grRules{pos_gene(k)};
                    if pos(1)==1
                        model_u.grRules{pos_gene(k)}    = uni_gene;
                    else
                        model_u.grRules{pos_gene(k)}    = [grRule(1:(pos(1)-1)),uni_gene];
                    end
                    if lenRule>lenGene
                        model_u.grRules{pos_gene(k)}    = [model_u.grRules{pos_gene(k)},grRule((pos(1)+lenGene):end)];
                    end    
                end
            end
        end
        % if all genes encoding for one reaction are exclusive simplify
        % gene rule
        if sum(is_exclusive_bin)==length(gene_inc)
            model_u.grRules{i}  = uni_gene;
        end
            
        % write relation matrix
        for j=1:length(is_exclusive)
            is_gene                     = strcmp(genes_save,model_u.genes{gene_inc(is_exclusive(j))});
            geneMatCompress{is_gene}    = uni_gene;
        end
        % rewrite gene name
        model_u.genes{gene_inc(is_exclusive(1))}    = uni_gene;
        keep_genes      = 1:length(model_u.genes);
        keep_genes(gene_inc(is_exclusive(2:end)))   = [];
        model_u.genes     = model_u.genes(keep_genes);
        % adapt gene reaction relation matrix
        model_u.rxnGeneMat    = model_u.rxnGeneMat(:,keep_genes);        
    end
end
model_u.geneMatCompress     = geneMatCompress;


end