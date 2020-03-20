function analyzeGAMOresults(model,results_file,target,varargin)
% Analyze GAMO results and extract relevant information from the simulated
% intervention strategies
%
% INPUT
% model:                    COBRA model
% results_file:             GAMO results file
% target:                   specifiy target type
%                               'rxns': reaction deletions
%                               'genes': gene deletions
%
% Optional
%     'fitMin':           Minimum relevant fitness
%     'numRelTarget':     Number of relevant targets
%     'plot':             Plot results (1):plot  (0): no plot
%     'targetNotation':   Choose type of target notation
%                           'name':    print target names (e.g. rxnNames)
%                           'ID'       print IDs (e.g. rxns)
%%


%% read input

% minimum relevant fitness
if any(strcmp(varargin,'fitMin'))
    fitMin  = varargin{find(strcmp(varargin,'fitMin'))+1};
else
    fitMin  = 1e-10;    % small positive number to exclude solutions with a fitness of zero
end
    
% plot results
if any(strcmp(varargin,'plot'))
    plotVar  = varargin{find(strcmp(varargin,'plot'))+1};
else
    plotVar  = 1;    % enable plot by default
end

% number of relevant targets
if any(strcmp(varargin,'targetNotation'))
    targetNotation  = varargin{find(strcmp(varargin,'targetNotation'))+1};
else
    targetNotation  = 'name';   
end

% number of relevant targets
if any(strcmp(varargin,'numRelTarget'))
    numRelTarget  = varargin{find(strcmp(varargin,'numRelTarget'))+1};
else
    numRelTarget  = 15;    % 
end

%% Parameter

%% extract results
fitness_best    = [results_file.res_best(:).fitness];
fit_best_pos    = find(fitness_best>=fitMin);    % choose relevant fitness increase
% choose target type
if strcmp(target,'rxns')
%     gene_del        = [results_file.res_best(fit_best_pos).KO];
    strategies_best = {results_file.res_best(fit_best_pos).KO}';
elseif strcmp(target,'genes')
%     gene_del        = [results_file.res_best(fit_best_pos).Gene_KO];
    strategies_best = {results_file.res_best(fit_best_pos).Gene_KO}';
else
    warning('Invalid target type! Genes are chosen as targets')
    target          = 'genes';
%     gene_del        = [results_file.res_best(fit_best_pos).Gene_KO];
    strategies_best = {results_file.res_best(fit_best_pos).Gene_KO}';
end
num_strat       = length(fit_best_pos);


%% determine unique solutions
unique_strategies_char   = {};
unique_strategies   = {};
for i=1:length(strategies_best)
    unique_strategies_char{end+1}    = num2str(sort(strategies_best{i}));
end
[unique_strategies_char,ia,~]   = unique(unique_strategies_char);
unique_strategies_fitness = fitness_best(fit_best_pos(ia));
for i=1:length(unique_strategies_char)
    unique_strategies{end+1}    = str2num(unique_strategies_char{i});
end
num_strat_unique    = length(unique_strategies);


%% count occurences of unique deletions in unique solutions
gene_del        = [unique_strategies{:}];
gene_del_unq    = unique(gene_del);
gene_del_count  = zeros(length(gene_del_unq),1);
for i=1:length(gene_del_unq)
    for j=1:num_strat_unique
        gene_del_count(i)   = gene_del_count(i)+any(ismember(unique_strategies{j},gene_del_unq(i)));        
    end
end

%% sort results
[gene_del_count_sort,P]   = sort(gene_del_count,'descend');
gene_del_unq_sort         = gene_del_unq(P);
% normalize with number of unique, best strategies
gene_del_count_sort_norm    = 100*(gene_del_count_sort/num_strat_unique);

%% sum of fitness for participating deletions
fitSum_gene_del_sort     = zeros(length(gene_del_unq_sort),1);
for i=1:length(unique_strategies)
    % find gene deletion matches
    match   = ismember(gene_del_unq_sort,unique_strategies{i});
    fitSum_gene_del_sort(match)  = fitSum_gene_del_sort(match)+unique_strategies_fitness(i);
end
fitMean_gene_del_sort   = fitSum_gene_del_sort./num_strat_unique;


%% save, plot and return

% check number of relevant targets
if length(gene_del_unq_sort)<numRelTarget
    numRelTarget    = length(gene_del_unq_sort);
end

 % rename targets
 if strcmp(target,'genes')
    % gene names
    target_names  = model.genes(gene_del_unq_sort(1:numRelTarget));
 elseif strcmp(target,'rxns')
     if strcmp(targetNotation,'name')
        % reaction names
        target_names  = model.rxnNames(gene_del_unq_sort(1:numRelTarget));
     elseif strcmp(targetNotation,'ID')
         % reaction IDs
         target_names  = model.rxns(gene_del_unq_sort(1:numRelTarget));
     end
 end
    
for i=1:length(target_names)
    target_names{i}   = [target_names{i}(1:2),'\',target_names{i}(3:end)];
end


if plotVar

    figure;
    hold on
    bar(gene_del_count_sort_norm(1:numRelTarget));
    set(gca,'XTick',1:numRelTarget)
    set(gca,'XTickLabel',target_names)
    xtickangle(45)
    ylabel('Occurence [%]')
    title(['Occurences in deletion strategies (N_{unique}=',num2str(num_strat_unique),')'])
    hold off

    figure;
    hold on
    bar(fitMean_gene_del_sort(1:numRelTarget));
    set(gca,'XTick',1:numRelTarget)
    set(gca,'XTickLabel',target_names)
    xtickangle(45)
    ylabel('Mean fitness')
    title('Mean fitness of designs including respective targets')
    hold off

    
   
end


end