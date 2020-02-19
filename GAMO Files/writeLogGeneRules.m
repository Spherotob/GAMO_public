function targets = writeLogGeneRules(model,targets,opt_ga)
% function writes m-File to evaluate gene rule logic for target reactions

% Input
% model:          metabolic model (will be altered, do not specify as an output)
% targetRxns:     Set of target reactions (numbers)
% targetGenes:    Set of target genes (identifier)
% KDID:           vector with entries showing if target is a KD target

% Output
% m-File to evaluate gene rule logic


%% general variables

targetRxns      = targets.rxnNum;
targetGenes     = targets.genes;
targetBounds    = targets.bound_gene;
KDID            = targets.KDID;
numTargetGenes  = length(targetGenes);

%% set up function 
%% Homogenize gene identifiers
targetsGenesRename  = cell(numTargetGenes,1);
for i=1:numTargetGenes
    targetsGenesRename{i}   = ['x',sprintf('%04d',i)];
    lenGene     = length(targetGenes{i});
    % change gene Rule identifier
    isGenePos   = strfind(model.grRules,targetGenes{i});
    isGene      = find(~cellfun('isempty',isGenePos));
    for j=1:length(isGene)
        genePos     = isGene(j);
        identPos    = isGenePos{genePos};
        lenRule     = length(model.grRules{genePos});
        for k=1:length(identPos)
            pos     = strfind(model.grRules{genePos},targetGenes{i});
            grRule  = model.grRules{genePos};
            if pos(1)==1
                model.grRules{genePos}  = [targetsGenesRename{i}];
            else
                model.grRules{genePos}  = [model.grRules{genePos}(1:(pos(1)-1)),targetsGenesRename{i}];
            end
            if lenRule>lenGene
                model.grRules{genePos}  = [model.grRules{genePos},grRule((pos(1)+lenGene):end)];
            end
        end
    end
    % change gene vector
    model.genes{strcmp(model.genes,targetGenes{i})}     = targetsGenesRename{i};
end
targetGenes     = targetsGenesRename;




%% write logical gene rules 
log_grRule  = cell(length(targetRxns) ,1); % string containing grRule in logical format
for i=1:length(targetRxns) 
    if isempty(model.grRules{targetRxns(i)})
%         warning('No gene annotated for reaction. No deletion possible')
        log_grRule{i}  = [log_grRule{i},'1'];  
        continue;
    end
    
    % distinguish between KD and KO target
    if KDID(i)
        num     = 2;
    else
        num     = 1;
    end
    
    for n=1:num       
            
        grRule          = strsplit(model.grRules{targetRxns(i)});   % get and split rule string
        log_grRule_s    = [];
        genes_rule      = [];   % involved genes in grRule
        for  j=1:length(grRule)
            % place logical operators
            if strcmp(grRule{j},'and')
                log_grRule_s  = [log_grRule_s,'&'];     % logical AND
            elseif strcmp(grRule{j},'or')
                log_grRule_s  = [log_grRule_s,'|'];     % logical OR
            else
                % set opening brackets
                ob  = grRule{j}=='(';
                for k=1:length(find(ob))
                    log_grRule_s  = [log_grRule_s,'('];     % opening bracket
                end
                grRule{j}   = grRule{j}(~ob);   % delete brackets

                % set gene variable
                cb          = grRule{j}==')';   % get position of closing brackets
                geneID      = grRule{j}(~cb);   % isolate gene ID
                
                if ~isempty(geneID)
                    % if only brackets were found skip gene search
                    if n==2
                        % search for knockdown gene target
                        geneID  = [geneID,'KD'];
                    end
                    genePos         = find(strcmp(targetGenes,geneID));     % position of gene in target vector
                    if isempty(genePos)
                        % gene not in target space so always functioning
                        log_grRule_s  = [log_grRule_s,'1'];  
                    else
                        log_grRule_s    = [log_grRule_s,['gene(',num2str(genePos),')']];
                        genes_rule      = [genes_rule,genePos];
                    end
                end
                  
                % set closing bracket
                for k=1:length(find(cb))
                    log_grRule_s  = [log_grRule_s,')'];     % closing bracket
                end
                
                grRule{j}   = grRule{j}(~cb);   % delete brackets      
            end 
        end
        % save KO rule
        if n==1
            log_grRule_KO   = log_grRule_s;
        elseif n==2
            log_grRule_KD   = log_grRule_s;
        end
        % process bounds of gene rule
        ub_rule     = zeros(1,length(genes_rule));  % upper flux bounds (absolute) of gene targets in grRule
        for j=1:length(genes_rule)
            if targetBounds(genes_rule(j),2)<0
                ub_rule(j)      = -targetBounds(genes_rule(j),1);
            else
                ub_rule(j)      = targetBounds(genes_rule(j),2);
            end
        end

            
        [targets.ub_rule{i},I]      = sort(ub_rule,'ascend');
            
        targets.genes_rule{i}       = genes_rule(I);
        
    end
    % merge, for KD targets, standard rule and KO rule
    if num==2
        log_grRule{i}   = ['(',log_grRule_KO,')<=(',log_grRule_KD,')'];
    elseif num==1
        log_grRule{i}   = log_grRule_KO;
    end
end


%% If only Knockouts are considered write reduced m-File
is_KD   = any(KDID);
switch is_KD
    case 0
        % Only KOs considered
        % write master function
        fileID  = fopen([opt_ga.AddFilesDir,opt_ga.slash,'evalTargets.m'],'w');
        fprintf(fileID,'function [targetRxnNum,targetRxnNum_i,targetRxnBounds,targetRxnBounds_i] = evalTargets(geneTargets,targets)\n');
        fprintf(fileID,['gene = ones(targets.Nt,1);\n',...
                        'gene(geneTargets(1:targets.K)) = 0;\n']);               
        fprintf(fileID,'%% Calculate reaction interventions from gene targets\n');   
%         fprintf(fileID,'tic\n');          
        fprintf(fileID,'targetRxnLog = evalGeneLog(gene);\n');
%         fprintf(fileID,'toc\n'); 
        
        
        if targets.Nt_hri>0
            fprintf(fileID,['%% Consider also heterologous insertions\n',...
                    'targetRxnLog = [targetRxnLog==0;false(',num2str(targets.Nt_hri),',1)];\n',...
                    'targetRxnLog(geneTargets(',num2str(targets.K+1),':end)+targets.shift) = true;\n']);
            fprintf(fileID,['targetRxnBounds_i_cell = targets.bound_i(targetRxnLog,:);\n',...
                    'targetRxnBounds_i = [[targetRxnBounds_i_cell{:,1}]'',[targetRxnBounds_i_cell{:,2}]''];\n',...
                    'targetRxnBounds = targets.bound(targetRxnLog,:);\n',...
                    'targetRxnNum = targets.rxnNum(targetRxnLog);\n',...
                    'targetRxnNum_i = [targets.rxnNum_i{targetRxnLog}]'';\n',...
                    'end\n']);     
        else
            fprintf(fileID,['if all(targetRxnLog)\n\t',...
                    'targetRxnNum = [];\n\t',...
                    'targetRxnNum_i = [];\n\t',...
                    'if length(geneTargets)<=targets.K\n\t\t',...
                    'targetRxnBounds = [0,0];\n\t\t',...
                    'targetRxnBounds_i = [0,0];\n\t',...
                    'else\n\t\t',...
                    'targetRxnBounds = [];\n\t\t',...
                    'targetRxnBounds_i = [];\n\t',...
                    'end\n',...
                    'else\n\t']);
            fprintf(fileID,'targetRxnLog = targetRxnLog==0;\n\t');
            fprintf(fileID,['targetRxnBounds_i_cell = targets.bound_i(targetRxnLog,:);\n\t',...
                    'targetRxnBounds_i = [[targetRxnBounds_i_cell{:,1}]'',[targetRxnBounds_i_cell{:,2}]''];\n\t',...
                    'targetRxnBounds = targets.bound(targetRxnLog,:);\n\t',...
                    'targetRxnNum = targets.rxnNum(targetRxnLog);\n\t',...
                    'targetRxnNum_i = [targets.rxnNum_i{targetRxnLog}]'';\n',...
                    'end\n']);   
            fprintf(fileID,'end');
        end
        fclose(fileID);
             
        % write gene logicals evaluation targets
        fileID  = fopen([opt_ga.AddFilesDir,opt_ga.slash,'evalGeneLog_MAT.m'],'w');
        fprintf(fileID,'function targetRxnLog = evalGeneLog(gene)\n');
        fprintf(fileID,'targetRxnLog = [');
        numRules    = length(log_grRule);
        for i=1:numRules
            if i~=numRules
                fprintf(fileID,[log_grRule{i},';...\n\t']);
            else
                fprintf(fileID,[log_grRule{i},'];\n']);
            end
        end
        fprintf(fileID,'end');
        fclose(fileID);
        % translate to a mex-file
        cd(opt_ga.AddFilesDir);
        codegen -o evalGeneLog ...
                evalGeneLog_MAT.m ...
                -args {zeros(length(targetGenes),1)}
        cd(opt_ga.homeDir)   
        % put m-file in the bin
%         delete(['AddFiles',slash,'evalGeneLog.m']);
        
        
    case 1
        % Consider also KDs

        % write rules for deriving KD bounds
        bound_rules     = cell(length(targetRxns) ,1);  % bound rules 
        pairs_s         = {};
        for i=1:length(targetRxns)
            if ~KDID(i)
                bound_rules{i}  = '0';
                continue;
            end

            % translate gene rules    
            grRule          = model.grRules{targetRxns(i)};   % get and split rule string
            % adding outer brackets
            grRule          = ['(',grRule,')'];
            split_grRule    = strsplit(grRule);
            % find logical operators and replace them
            f_op_and                = find(~cellfun('isempty',strfind(split_grRule,'and')));
            split_grRule(f_op_and)  = {','};
            f_op_or                 = find(~cellfun('isempty',strfind(split_grRule,'or')));
            split_grRule(f_op_or)  = {'*'};
            % merge string
            bound_rule_s    = strjoin(split_grRule);
            bound_rule_u    = bound_rule_s;
            % process AND logicals
            f_op_and    = strfind(bound_rule_s,',');
            for j=1:length(f_op_and)
                f_op_and    = strfind(bound_rule_s,',');
                % search for closest opening bracket
                n   = 1;    % number of closing brackets encountered
                for k=1:f_op_and(j)
                    pos     = f_op_and(j)-k;    % position in string
                    if strcmp(bound_rule_s(pos),')')
                        n   = n+1;
                    elseif strcmp(bound_rule_s(pos),'(')
                        n   = n-1;
                    elseif strcmp(bound_rule_s(pos),'[') && n==1
                        % min is already set
                        break;
                    end
                    if n==0
                        % found position of opening bracket
                        bound_rule_s    = [bound_rule_s(1:(pos-1)),'min([',bound_rule_s((pos+1):end)];
                        break;
                    end
                end

                % search for closest closing bracket
                f_op_and    = strfind(bound_rule_s,',');
                n   = 1;    % number of opening brackets encountered
                for k=1:(length(bound_rule_s)-f_op_and(j))
                    pos     = f_op_and(j)+k;    % position in string
                    if strcmp(bound_rule_s(pos),'(')
                        n   = n+1;
                    elseif strcmp(bound_rule_s(pos),')')
                        n   = n-1;
                    elseif strcmp(bound_rule_s(pos),']') && n==1
                        % min is already set
                        break;
                    end
                    if n==0
                        % found position of closing bracket
                        bound_rule_s    = [bound_rule_s(1:(pos-1)),'])',bound_rule_s((pos+1):end)];
                        break;
                    end
                end
            end

            % replace gene identifiers
            f_gene  = strfind(bound_rule_s,'x');
%             f_gene  = [f_gene,strfind(bound_rule_s,'s')];
            % extract gene ID
            geneID  = cell(length(f_gene),1);
            for j=1:length(f_gene)
                geneID{j}   = bound_rule_s(f_gene(j):(f_gene(j)+4));
            end
            % change with array name containing gene bounds
            for j=1:length(geneID)
                geneID_pos   = strfind(bound_rule_s,geneID{j});
                gene_pos     = find(strcmp(targetGenes,[geneID{j},'KD']));
                if isempty(gene_pos)
                    bound_rule_s    = [bound_rule_s(1:(geneID_pos-1)),'1',...
                                    bound_rule_s((geneID_pos+5):end)];
                else
                    bound_rule_s    = [bound_rule_s(1:(geneID_pos-1)),'ub_gene(',num2str(gene_pos),')',...
                                    bound_rule_s((geneID_pos+5):end)];
                end

            end
            bound_rules{i}  = bound_rule_s;


        end

        % Berücksichtige nur upper bounds von KD targets
        % integriere das in die zu erstellende m-File, berücksichtige irreversibles
        % Modell


        %% write m file
        fileID  = fopen(['AddFiles',slash,'evalTargets.m'],'w');
        fprintf(fileID,'function [targetRxnNum,targetRxnNum_i,targetRxnBounds,targetRxnBounds_i] = evalTargets(geneTargets,targets)\n');
        fprintf(fileID,['gene = ones(targets.Nt,1);\n',...
                        'gene(geneTargets) = 0;\n']);
        fprintf(fileID,'%% Calculate reaction interventions from gene targets\n');           
        fprintf(fileID,'targetRxnLog = [');
        numRules    = length(log_grRule);
        for i=1:numRules
            if i~=numRules
                fprintf(fileID,[log_grRule{i},';...\n\t']);
            else
                fprintf(fileID,[log_grRule{i},'];\n']);
            end
        end
        fprintf(fileID,['if all(targetRxnLog)\n\t',...
                            'targetRxnNum = [];\n\t',...
                            'targetRxnNum_i = [];\n\t',...
                            'targetRxnBounds = [0,0];\n\t',...
                            'targetRxnBounds_i = [0,0];\n\t',...
                            'return\n',...
                            'end\n']);
        fprintf(fileID,'targetRxnLog = targetRxnLog==0;\n');
        fprintf(fileID,'%% Calculate upper reaction bounds\n');
        fprintf(fileID,'ub_rxns_target = [');
        for i=1:numRules
            if i~=numRules
                fprintf(fileID,[bound_rules{i},';...\n\t']);
            else
                fprintf(fileID,[bound_rules{i},'];\n']);
            end
        end
        fprintf(fileID,'%% Determine numbers of reaction targets\n');
        fprintf(fileID,['targetRxnNum = targets.rxnNum(targetRxnLog);\n',...
                           'targetRxnNum_i = [targets.rev2irrev{targetRxnNum}]'';\n']);

        fprintf(fileID,'%% Determine target reaction bounds\n');  
        fprintf(fileID,'targetRxnBounds = [zeros(sum(targetRxnLog),1),targets.fd_ref(targetRxnNum).*ub_rxns_target(targetRxnLog)];\n'); 
        fprintf(fileID,'ub_i = targets.B_t_bound(:,targetRxnLog)*targetRxnBounds(:,2);\n');
        fprintf(fileID,'targetRxnBounds_i = [zeros(length(targetRxnNum_i),1),ub_i(targetRxnNum_i)];\n');  

        % fprintf(fileID,['gene = gene==0;\n',...
        %                     ]);
        % 
        % fprintf(fileID,'targetRxnNum = targets.rxnNum(targetRxnLog);\n');
        % fprintf(fileID,'[min_ub_bound,min_ub_bound_pos] = min(abs(targets.ub_gene_rxn_mat(gene,targetRxnLog)),[],1);\n');
        % fprintf(fileID,'lb_gene_rxn_mat_red = targets.lb_gene_rxn_mat(gene,:);\n');
        % fprintf(fileID,'targetRxnBounds = lb_gene_rxn_mat_red(min_ub_bound_pos,targetRxnLog);\n');
        % % fprintf(fileID,'targetRxnNum_i = targets.rxnNum_i(targetRxnLog);\n');
        % fprintf(fileID,'targetNum = find(targetRxnLog);\n');
        % fprintf(fileID,'numTarget = length(targetNum);\n');
        % fprintf(fileID,'targetRxnBounds = zeros(numTarget,2);\n');
        % fprintf(fileID,['for i=1:numTarget\n\t',...
        %                 'genes_rule = targets.genes_rule{targetNum(i)};\n\t',...
        %                 'genes_rule(~ismember(genes_rule,geneTargets)) = [];\n\t',...
        %                 'targetRxnBounds(i,:) = targets.bounds(genes_rule(1),:);\n',...
        %                 'end\n']);
        fprintf(fileID,'end');

        fclose(fileID);
end

end

