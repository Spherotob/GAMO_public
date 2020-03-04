function [pop,pop_Tbin, pop_bin,enc,targets] = initializePopulation(model,targets,popSize,b,K,K_hri,optFocus,threads,initPopType,opt_ga)

% Create initial population (generation) 

% INPUT
% model:        Metabolic model in COBRA format

% targets:      Struct containing all possible deletion or knockdown targets
%               including reaction numbers, scores and bounds (rMMA) as well
%               as target map.

% threads:      Number of active parallel workers

% b:            Selection pressure for target ranking

% K:            Maximal number of interventions

% initPopType:  Type of initial population (0): randomly select interventions; (1) only single intervention mutants

% OUTPUT
% pop:        initial Population

% pop_Tbin:   Initial population in target binary format 

% pop_bin:    Initial population in binary format

% enc:        struct containing encoding information

%% General parameter
elite   = opt_ga.elite;
K_tot   = K+K_hri;                  % total number of targets
Nt      = length(targets.rxnNum);   % Number of targets (reactions)
Nt_hri  = targets.Nt_hri;   % number of insertion targets
Np      = popSize*threads;          % Number of chromosomes 
RN      = rand(Np,K_tot);         % Array of uniformily distributed random numbers
a       = 2-b;                      % Selective pressure for the worst individual
critSampleRate  = 100;              % Critical maximal sample rate
discrFac = 50;                      % Discretization factor for encoding targets as binaries

%% Use gene or reaction deletions
if strcmp(optFocus,'gene')
    % search for gene interventions
    rxnGeneMat  = model.rxnGeneMat;
    % find target genes according to target reactions
    g_targets_KD            = [];
    g_targets_KO            = [];
    g_targets_KD_relBound      = [];    % bound relative to wildtype flux
    numTargets              = 0;    % number of gene targets
    for i=1:Nt
        if any(targets.bound(i,:))
            % Knockdown target, every target is unique due to differing
            % bounds
            targetCluster           = find(rxnGeneMat(targets.rxnNum(i),:));
            g_targets_KD            = [g_targets_KD;targetCluster(:)];
            for j=1:length(targetCluster)
                g_targets_KD_relBound   = [g_targets_KD_relBound;targets.bound(i,:)./model.fd_ref(targets.rxnNum(i))];
            end
            numTargets              = numTargets+length(targetCluster);
        else
            % KO target
            g_targets_KO           = [g_targets_KO,find(rxnGeneMat(targets.rxnNum(i),:))];  % necessary gene targets for reaction deletion
        end 
    end
    g_targets_KO    = unique(g_targets_KO)';     % erase duplicates in gene targets
    
    % number of gene targets
    numTargets      = numTargets+length(g_targets_KO);
    
    % cover all possible interventions for each target gene and construct
    % upper bound matrix
    r_targets_KO    = [];
    r_targets_KD    = [];

    for i=1:length(g_targets_KO)
        r_targets_KO    = [r_targets_KO;find(rxnGeneMat(:,g_targets_KO(i)))];
    end
    for i=1:length(g_targets_KD)
        r_targets_KD    = [r_targets_KD;find(rxnGeneMat(:,g_targets_KD(i)))];
    end
    % merge target vectors
    r_targets_KO    = unique(r_targets_KO);
    r_targets_KD    = unique(r_targets_KD);
    targets.rxnNum  = [r_targets_KO;r_targets_KD];
    targets.KDID    = [zeros(length(r_targets_KO),1);ones(length(r_targets_KD),1)];     % marks KD targets with '1'
    targets.geneNum  = [g_targets_KO;g_targets_KD];
    % create new reaction target bound vector ONLY KNOCKOUTS CONSIDERED!!!!!!!!!!
    targets.bound   = zeros(length(targets.rxnNum),2);
    
    % number of reaction targets
    numRxnTargets   = length(targets.rxnNum);
    
     
    % create vector with gene IDs and bounds
    genes    = cell(numTargets,1);
    bounds   = zeros(numTargets,2);
    c       = 1;
    % KO targets
    for i=1:length(g_targets_KO)
        genes{c}     = model.genes{g_targets_KO(i)};
        c   = c+1;
    end
    % KD targets
    for i=1:length(g_targets_KD)
        bounds(c,:)     = g_targets_KD_relBound(i,:);
        genes{c}        = [model.genes{g_targets_KD(i)},'KD'];     % mark gene as KD target
        c   = c+1;
    end
    
    % construct upper and lower bound gene matrix
    ub_gene_rxn_mat = NaN(numTargets,numRxnTargets);
    lb_gene_rxn_mat = ub_gene_rxn_mat;
    for i=1:numTargets
        for j=1:numRxnTargets
            if rxnGeneMat(targets.rxnNum(j),targets.geneNum(i))
                ub_gene_rxn_mat(i,j)    = rxnGeneMat(targets.rxnNum(j),targets.geneNum(i))*bounds(i,2);
                lb_gene_rxn_mat(i,j)    = rxnGeneMat(targets.rxnNum(j),targets.geneNum(i))*bounds(i,1);
            end
        end
    end
    
    targets.lb_gene_rxn_mat     = lb_gene_rxn_mat;
    targets.ub_gene_rxn_mat     = ub_gene_rxn_mat;
    
    targets.genes           = genes;
    targets.bound_gene      = bounds;
    targets.Nt              = length(genes);   

    
    % write m-File for gene logic evaluation in additional folder
    targets     = writeLogGeneRules(model,targets,opt_ga);
    
    % recalculate number of targets
    Nt      = length(targets.genes);  
    
    % target map (überflüssig)
    targets.map     = eye(Nt);
    targetMap       = targets.map;
    
    % assign uniform scores to gene targets
    score               = ones(Nt,1)./Nt;
    targets.score       = score;
    
    
    c = fix(clock);
    disp([num2str(c(4:end)),': Number target genes: ',num2str(length(targets.genes))])
else
    % search for reaction interventions
    % load scores
    score       = targets.score;
    % load target map
    targetMap   = targets.map;
    % Allocate bound vaector to cell
%     targets.bound   = num2cell(targets.bound);

    % write m-File to get target numbers and bounds
    writeRxnRules(opt_ga);
        
    targets.Nt      = Nt;
end

targets.rxnNum_KO  = targets.rxnNum;
% declare a shift operator for mapping genes of an individual) to actual targets
targets.shift   = length(targets.rxnNum)-Nt;

% Add heterologous reaction insertion targets
targets.rxnNum  = [targets.rxnNum;targets.rxnNum_hri];
targets.bound   = [targets.bound;targets.bound_hri];
targets.score   = [targets.score;targets.score_hri];


% calculate number of targets
targets.Nt      = Nt;
targets.Nt_tot  = Nt+targets.Nt_hri;


%% Checkings
% enough combinations for population size?
maxComb     = 0;
for i=1:K
    maxComb     = maxComb + nchoosek(Nt,i);
end
if maxComb < Np
    warning(['Population size exceeds maximal number of target combinations. Maximal combinations: ',num2str(maxComb)]);
    % reduce population size
    popSize     = floor(maxComb/threads);
    Np          = popSize*threads;
    RN          = rand(Np,K_tot);
    % check new population size
    if popSize<=elite
        error('Population size could not be reduced appropriately. Reduce number of threads if possible.')
    else
        disp(['New populations size: ',num2str(popSize)])
    end        
end

%% set up first population
% sort targets according to the score
[~,P]   = sort(score,'ascend');
% assign probability to each target (rank weighting or equal probability)
if any(abs(score-score(1))>1e-4)
    % rank weighting
    W       = (a+(([1:Nt]'./Nt).*(b-a)))./(Nt+1);
    sumW    = cumsum(W);
    sumW2   = cumsum(W-W(1));
else
    % uniform probability
    W       = ones(Nt,1)./Nt;
    sumW    = cumsum(W);
    sumW2   = sumW-W(1);
end

% assign targets to population
pop     = zeros(Np,K_tot);
pop_Tbin = zeros(Np,Nt+Nt_hri);     % binary representation of the population

blank               = zeros(1,Nt);  % blank binary chromosome
sampleRate_save     = zeros(Np,1);

% consider user-defined default targets
numDefaultSolutions     = 0;    % total number of default solutions
if ~isempty(opt_ga.defaultTargets)
    % check format of default targets field
    if iscell(opt_ga.defaultTargets)
        % valid cell array containing gene or reaction identifier
        for i=1:length(opt_ga.defaultTargets)
            defaultSolution = [];  % default encoded solution 
            % each entry is a unique solution
            defaultSolutionIDs     = opt_ga.defaultTargets{i};
            for j=1:length(defaultSolutionIDs)
                % identify target type
                if any(strcmp(model.genes,defaultSolutionIDs{j}))
                    % gene target
                    % get position in target vector
                    if strcmp(optFocus,'gene')            
                        defaultSolution   = [defaultSolution;find(ismember(targets.geneNum,...
                            find(strcmp(model.genes,defaultSolutionIDs{j}))))];
                    else
                        % match reaction to gene(s)
                        rxnPos  = find(model.rxnGeneMat(:,...
                            find(strcmp(model.genes,defaultSolutionIDs{j}))));
                        defaultSolution     = [defaultSolution;...
                            find(ismember(targets.rxnNum,rxnPos))];
                    end
                elseif any(strcmp(model.rxns,defaultSolutionIDs{j}))
                    % reaction target
                    % get position in target vector
                    if strcmp(optFocus,'rxns')            
                        defaultSolution   = [defaultSolution;find(ismember(targets.rxnNum,...
                            find(strcmp(model.rxns,defaultSolutionIDs{j}))))];
                    else
                        % match gene to reaction(s)
                        genePos  = find(model.rxnGeneMat(...
                            find(strcmp(model.rxns,defaultSolutionIDs{j})),:));
                        defaultSolution     = [defaultSolution;...
                            find(ismember(targets.geneNum,genePos))];
                    end                  
                end
            end
            defaultSolution     = unique(defaultSolution);
            % check length of default solution
            lenDefaultSolution  = length(defaultSolution);
            if lenDefaultSolution<=K && ~isempty(defaultSolution)
                numDefaultSolutions     = numDefaultSolutions+1;
                pop(numDefaultSolutions,1:lenDefaultSolution)    = defaultSolution';
                % fill remaining empty target spaces
                for k=(lenDefaultSolution+1):K
                    % pick a random target of the default solution
                    pop(numDefaultSolutions,k)    = pop(i,randperm(lenDefaultSolution,1));
                end
                
            else
                if strcmp(optFocus,'gene')  
                    warning(['Default solution ',num2str(i),...
                        ': Translates to more gene targets than allowed (or none)']) 
                else
                    warning(['Default solution ',num2str(i),...
                        ': Translates to more reaction targets than allowed (or none)'])   
                end
            end
        end
        
        
    else
        warning('Default targets field is not in the correct format. Provide a cell array containing gene or reaction identifiers!');
    end  
end

switch initPopType
    case 0
        % randomly select interventions
        for i=(numDefaultSolutions+1):Np
            % sample until a valid chromosome is found
            sample      = 1;
            sampleRate  = 1;
            while sample
                % Choose target for a chromosome
                for j=1:K
                    pop(i,j)    = P(sumW>=RN(i,j) & sumW2<RN(i,j));
                end

                % check if each intervention tackles a unique reaction
                rxnFlag             = blank;
                rxnFlag(pop(i,1:K))   = 1;
                if any(rxnFlag*targetMap > 1)
                    % resample chromosome
                    RN(i,:)     = rand(1,K_tot);
                    sampleRate  = sampleRate+1;
                    continue;
                else
                    % chromosome valid, continue
                end

                
                % add target insertions (choose randomly)
                for j=(K+1):(K_tot)
                    pop(i,j)    = Nt+(ceil(RN(i,j)*Nt_hri)+(RN(i,j)==0));
                end
                
                pop_Tbin(i,pop(i,:))     = 1;    % save chromosome in binary code
                
                if i>1
                    % check for chromosome duplicates
        %             A   = pop_Tbin(1:(i-1),:);
        %             for j=1:(i-1)
        %                 A(j,:)  = A(j,:)==pop_Tbin(i,:);
        %             end
        %             if any(all(A,2))    
                    if any(all(pop_Tbin(1:(i-1),:)==pop_Tbin(i,:),2))  
                        % resample chromosome
                        RN(i,:)     = rand(1,K_tot);
                        sampleRate  = sampleRate+1;
                    else
                        % chromosome valid, sample next chromosome
                        sample  = 0;
                    end
                else
                    sample = 0;
                end

                if critSampleRate < sampleRate
                    % critical sample rate exceeded
                    warning('Critical sample rate for population initialization exceeded!')
                    % Now only print a warning, take chromosome
                    sample  = 0;
                end
                if sample==0
                    sampleRate_save(i)  = sampleRate;
                end
            end
        end
        
    case 1
        % generate only single knockout mutants for the initial population
%         t1   = 1;    % target number counter
%         t2   = Nt+1;
%         for i=1:Np
%             for j=1:K
%                 pop(i,j)        = t1;    % write population file
%                 pop_Tbin(i,t1)   = 1;    % code population in binary form
%             end
%             t1   = t1+1;  % advance counter
%             % add insertion targets
%             for j=(K+1):K_tot
%                 pop(i,j)        = t2;
%                 pop_Tbin(i,t2)   = 1; 
%             end
%             t2  = t2+1;
%             
%             % check if target counter exceeds number of targets
%             if t1>Nt
%                 t1   = 1;
%             end
%             if t2>(Nt+Nt_hri)
%                 t2  = Nt+1;
%             end
%         end
        
        % randomized routine
        randTargetNum   = (rand(Np,1).*(Nt-1))+1;
        randInsertNum   = (rand(Np,1).*(Nt_hri-1))+1+Nt;
        for i=(numDefaultSolutions+1):Np
            % add deletions
            for j=1:K
                pop(i,j)                        = randTargetNum(i);    % write population file
                pop_Tbin(i,randTargetNum(i))    = 1;    % code population in binary form
            end
            % add insertions
            for j=(K+1):K_tot
                pop(i,j)                        = randInsertNum(i);
                pop_Tbin(i,randInsertNum(i))    = 1; 
            end
        end
        
    otherwise
        error('Choose how initial population should be constructed')
end

%% Translate population to binary format
% choose binary factor
% chosse sufficient number of bits to minimize bias during mutations
numBits     = round(log(discrFac*Nt)/log(2));           % number of bits
encodeVec   = ceil([1:(2^numBits)]./((2^numBits)/Nt));  % (d)encoding vector matching decimal to an intervention
% save start position of genes in binary population format
posGene     = 1:numBits:((numBits*K)-numBits+1);


% choose sufficient number of bits for insertion targets
if Nt_hri>0
    numBits_hri     = round(log(discrFac*Nt_hri)/log(2));
    encodeVec_hri   = ceil([1:(2^numBits_hri)]./((2^numBits_hri)/Nt_hri));
    % combine encoding vectors
    encodeVec       = [encodeVec,encodeVec_hri+Nt];
    % add start positions of gene
    posGene     = [posGene,[1:numBits_hri:((numBits_hri*K_hri)-numBits_hri+1)]+(numBits*K)];
else
    numBits_hri     = 0;
end



% for encoding choose middle discrete value for an intervention
encodeVec_midPos        = zeros(Nt+Nt_hri,1);
for i=1:(Nt+Nt_hri)
    discr                   = find(encodeVec==i);
    encodeVec_midPos(i)     = round((discr(1)+discr(end))/2);
end
% % do the same for insertion targets
% encodeVec_midPos_hri    = zeros(Nt_hri,1);
% for i=1:Nt_hri
%     discr                       = (2^numBits)+find(encodeVec(((2^numBits)+1):end)==(i+Nt));
%     encodeVec_midPos_hri(i)     = round((discr(1)+discr(end))/2);
% end


% encode population
pop_bin     = encode(pop,numBits,numBits_hri,encodeVec_midPos,K,K_hri);

% save and parse encoding information
enc.encodeVec           = encodeVec;
enc.encodeVec_midPos    = encodeVec_midPos;
% enc.encodeVec_midPos_hri = encodeVec_midPos_hri;
enc.numBits             = numBits;
enc.numBits_hri         = numBits_hri;
enc.posGene             = posGene;
enc.K                   = K;
enc.K_hri               = K_hri;
enc.K_tot               = K_tot;
enc.numBits_tot         = (numBits*K)+(numBits_hri*K_hri);
% enc.numRxn_KO           = numRxn_KO;

end