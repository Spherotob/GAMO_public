function [results,prob] = GAMO(model,opt,opt_fitFun,prob)

% starting routine for the Genetic Algorithm for Metabolic Optimization (GAMO)
% GAMO includes reverse Minimal Metabolite Balance Analysis (rMMA)

% INPUT
% model:        Metabolic model in COBRA format additionally specifying the
%               substrate uptake reaction (subsRxn), the target compound exchange
%               reaction (targetRxn), biomass reaction (bmRxn) and a reference flux distribution
%               (refFluxDist).
% opt:          Struct including general options for GAMO
% opt_fitFun:   Struct including fitness function options for GAMO
% prob:         Struct containing information about the same problem
%               tackled before

% OUTPUT
% results:      compact results including the best identified strain designs
% prob:         contains all information to restart optimization from the latest population

% .. AUTHOR
%     Tobias Alter    26/03/2018
%       - Institute of Applied Microbiology, RWTH Aachen University

c = fix(clock);
disp([num2str(c(4:end)),': Starting GAMO ...'])


%% Check opt


% Check if particular options have changed
fitFun_chngd    = 1;
if isfield(prob,'opt')
    % fitness function type
    if isfield(opt,'fitFun')
        if opt.fitFun==prob.opt.fitFun
            fitFun_chngd    = 0;
        end    
    else
        fitFun_chngd    = 0;
    end
end
% load previous problem
if isfield(prob,'model')
    model   = prob.model;
end
if isfield(prob,'model_i')
    model_i     = prob.model_i;
end
    
% save for rerunning problem
prob.opt    = opt;

% load previous total data struct
if isfield(prob,'totalData')
    if isempty(prob.totalData)
        pop_gen_tot         = [];
        popFit_gen_tot      = [];
        popFitObj_gen_tot   = [];
        timing_gen_tot      = [];
        timing_drifts_tot   = [];
    else
        pop_gen_tot         = prob.totalData.pop;
        popFit_gen_tot      = prob.totalData.popFit;
        popFitObj_gen_tot   = prob.totalData.popFitObj;
        if isfield(prob.totalData,'timing_gen')
            timing_gen_tot      = prob.totalData.timing_gen;
            timing_drifts_tot   = prob.totalData.timing_drifts;
        else
            timing_gen_tot      = cell(size(pop_gen_tot,1),1);
            timing_drifts_tot   = zeros(size(pop_gen_tot,1),1);
        end
    end
else
    pop_gen_tot         = [];
    popFit_gen_tot      = [];
    popFitObj_gen_tot   = [];
    timing_gen_tot      = [];
    timing_drifts_tot   = [];
end

% general options
if isfield(opt,'saveFile')
    if ischar(opt.saveFile)
        saveFile   = opt.saveFile;
    else
        saveFile    = [];
    end
else
    saveFile    = [];
end

if isfield(opt,'saveFolder')
    saveFolder  = opt.saveFolder;
else
    saveFolder  = [];
end

% fitness function
if isfield(opt_fitFun,'excl_rxns')
    % include MiMBl options 
    opt_rMMA.MiMBl.excl_rxns    = opt_fitFun.excl_rxns;
else
    opt_fitFun.excl_rxns     = [];
    opt_rMMA.MiMBl.excl_rxns = [];
end

if ~isfield(opt,'fitFun')
    opt.fitFun  = 0;
end

% special fitness function options
if opt.fitFun==1 % multiobjective function
    if ~isfield(opt_fitFun,'isObj')
        opt_fitFun.isObj    = 'GRO';
    end
    if ~isfield(opt_fitFun,'leadObj')
        opt_fitFun.leadObj    = 'M';
    end
    if ~isfield(opt_fitFun,'weighting')
        opt_fitFun.weighting    = [1,1,1,1];    % uniform weightings
    else
        if size(opt_fitFun.weighting,2)~=4
            error('Pass weighting factors for each objective function')
        end  
    end
end 
% minimal allowable growth rate
if ~isfield(opt_fitFun,'minGrowth')
    opt_fitFun.minGrowth    = 0;
end
% minimize number of interentions
if isfield(opt_fitFun,'minInt')
    if (opt_fitFun.minInt~=0) && (opt_fitFun.minInt~=1)
        warning('Minimization of interventions not initialized correctly! Not considered!');
        opt_fitFun.minInt   = 0;
    end
else
    opt_fitFun.minInt  = 0;
end
    

% check opt struct
if isfield(opt,'modelType')
    modelType   = opt.modelType;
else
    modelType   = 0;    % Default: Bigg Model
end

% optimization focus
if isfield(opt,'optFocus')
    if  strcmp(opt.optFocus,'gene')
        if isfield(model,'rxnGeneMat') && isfield(model,'grRules')
            optFocus    = opt.optFocus;
        else
            error('No gene rules associated with the provided model')
        end
    elseif strcmp(opt.optFocus,'rxns')
        optFocus    = opt.optFocus;
    else
        optFocus    = 'rxns';
    end
else
    optFocus    = 'rxns';
end
if strcmp(optFocus,'gene')
    opt_fitFun.optFocus     = 1;
elseif strcmp(optFocus,'rxns')
    opt_fitFun.optFocus     = 0;
end


% model related
if isfield(opt,'rMMA_flag')
    rMMA_flag     = opt.rMMA_flag;
else
    rMMA_flag     = 0; 
end
if isfield(opt,'redFlag')
    redFlag     = opt.redFlag;
else
    redFlag     = 0; 
end
if isfield(opt,'compress')
    compress     = opt.compress;
else
    compress     = 0; 
end
if isfield(opt,'numInt')
    K       = opt.numInt;
else
    K       = 3; 
end
if isfield(opt,'typeInt')
    typeInt     = opt.typeInt;
else
    typeInt       = 0; 
end

% genetic algorithm properties
if isfield(opt,'memPop')
    opt_ga.memPop   = opt.memPop;
else
    opt_ga.memPop   = 1;
end

if isfield(opt,'popSize')
    popSize    = opt.popSize;
else
    popSize    = 20; 
end
if isfield(opt,'slctPressure')
    b    = opt.slctPressure;
else
    b    = 2; 
end
if isfield(opt,'fitFun')
    fitFun_type     = opt.fitFun;
else
    error('No fitness function specified') 
end
if isfield(opt,'genSize')
    opt_ga.genSize      = opt.genSize;
else
    opt_ga.genSize      = 5;
end
if isfield(opt,'maxGen')
    opt_ga.maxGen      = opt.maxGen;
else
    opt_ga.maxGen      = 100;
end
if isfield(opt,'slctRate')
    opt_ga.slctRate    = opt.slctRate;
else
    opt_ga.slctRate    = 0.5;
end
if isfield(opt,'noMatingIdent')
    opt_ga.noMatingIdent    = opt.noMatingIdent;
else
    opt_ga.noMatingIdent    = 1;
end
if isfield(opt,'numKntChr')
    opt_ga.numKntChr    = opt.numKntChr;
else
    opt_ga.numKntChr    = 1;
end
if isfield(opt,'crossType')
    opt_ga.crossType    = opt.crossType;
else
    opt_ga.crossType    = 0;
end

% mutation options
if isfield(opt,'mutRate')
    opt_ga.mutRate    = opt.mutRate;
else
    opt_ga.mutRate    = 0.2;
end

% advanced settings
if isfield(opt,'elite')
    opt_ga.elite    = opt.elite;
else
    opt_ga.elite    = 2;
end


if isfield(opt,'initPopType')
    initPopType     = opt.initPopType;
else
    initPopType     = 1;
end

if isfield(opt,'hetRxnNum') && ~isfield(model,'hetRxnNum')
    model.hetRxnNum   = opt.hetRxnNum;
    % reference flux distribution must be provided
    if ~isfield(model,'fd_ref')
        error('Reference flux distribution must be provided. Add field "fd_ref" to model structure')
    end
else
    if ~isfield(model,'hetRxnNum')
        model.hetRxnNum   = [];
    end
end
if isempty(model.hetRxnNum)
    K_hri   = 0;
else
    if isfield(opt,'numInsertions')
        K_hri   = opt.numInsertions;
    else
        K_hri   = 1; 
    end
end


% non target reactions
if isfield(opt,'nonTarget')
    nonTarget     = opt.nonTarget;
else
    nonTarget     = [];
end


% technical
if isfield(opt,'threads')
    threads         = opt.threads;
    opt_ga.threads  = threads;
else
    threads         = 1;
    opt_ga.threads  = threads;
end


%% Additional preprocessing

% establish home directory and fodlers
if ispc
    slash   = '\';
else
    slash   = '/';
end
homeDir             = pwd;
opt_ga.homeDir      = homeDir;
opt_ga.AddFilesDir  = [homeDir,slash,'AddFiles'];
opt_ga.slash    = slash;

if exist(opt_ga.AddFilesDir,'dir')~=7
    mkdir(opt_ga.AddFilesDir)
end
addpath(opt_ga.AddFilesDir)

% Rescue Folder for saving files
if exist([homeDir,slash,'Rescue'],'dir')~=7
    mkdir([homeDir,slash,'Rescue'])
end


%% Re-load optimization problem




%% General parameter
numRxn  = length(model.rxns);   % number of reactions
numMet  = length(model.mets);   % number of metabolites


%% process model 
model.bmRxnNum      = find(ismember(model.rxns,model.bmRxn));
model.subsRxnNum    = find(ismember(model.rxns,model.subsRxn));
model.targetRxnNum  = find(ismember(model.rxns,model.targetRxn));



%% start parallel cluster
p   = gcp('nocreate');
if isempty(p)
    if threads <= 0
        % start default parallel cluster
        parpool;
    else
        % start paralel cluster
        parpool(threads);
    end
else
    if p.NumWorkers ~= threads
        delete(p)   % close parallel cluster
        parpool(threads)    % start new parallel cluster
    end
end
p           = gcp('nocreate');
threads     = p.NumWorkers;
disp(['Number of parallel workers: ',num2str(p.NumWorkers)])





%% reference flux distribution
c = fix(clock);
disp([num2str(c(4:end)),': Calculating reference flux distribution ...'])

% Check provided prob struct
refFlux_flag    = 0;
if isfield(prob,'model')
    % model of previous problem exists
    model   = prob.model;
    if ~isfield(model,'fd_ref')
        refFlux_flag    = 1;
    end
else
    refFlux_flag    = 1;
end
        
if refFlux_flag
    rxnNumBM    = find(ismember(model.rxns,model.bmRxn));
    if isfield(model,'fd_ref')
        if isempty(model.fd_ref)
            [refFD, ~]  = predictRefFluxDist(model,rxnNumBM,0);
        else
            refFD   = model.fd_ref;
        end
    else
        [refFD, ~]  = predictRefFluxDist(model,rxnNumBM,0);    
    end
    model.fd_ref   = refFD;
end

%% translate non target reaction identifiers to positions
if ~isempty(nonTarget)
    nonTargetNum    = zeros(length(nonTarget),1);
    for i=1:length(nonTarget)
        nonTargetNum(i)     = find(strcmp(model.rxns,nonTarget{i}));
    end
    % swap identifiers witch number
    nonTarget   = nonTargetNum;
else
    nonTarget   = [];
end


%% compress model
% save model
model_save          = model;
if ~isfield(prob,'model_compress') && compress
        
        % save uncompressed model options 
        opt_fitFun_save         = opt_fitFun;
        prob.opt_fitFun_save    = opt_fitFun_save;
        model                   = compressModel(model);     % compress model
        
        % translate reference flux distribution, excluded
        % reactions, heterologous insertion reactions
        model.fd_ref                = model.comprMapMat'*model.fd_ref;
        excl_rxns_help              = zeros(length(model_save.rxns),1);
        excl_rxns_help(opt_fitFun.excl_rxns)   = 1;
        opt_fitFun.excl_rxns        = find(model.comprMapMat'*excl_rxns_help);
        opt_rMMA.MiMBl.excl_rxns    = opt_fitFun.excl_rxns;
        model.hetRxnNum             = find(ismember(model.comprMapVec,model.hetRxnNum));
        nonTarget                   = find(ismember(model.comprMapVec,nonTarget));
        
        model.bmRxnNum      = find(model.comprMapVec==model_save.bmRxnNum);
        model.subsRxnNum    = find(model.comprMapVec==model_save.subsRxnNum);
        model.targetRxnNum  = find(model.comprMapVec==model_save.targetRxnNum);
        
        c = fix(clock);
        fprintf([num2str(c(4:end)),': Size of compressed model:\n\t\t\t',num2str(length(model.rxns)),...
                    ' reactions\n\t\t\t',num2str(length(model.mets)),' metabolites\n\t\t\t',...
                    'CF reactions: ',num2str(round(length(model.rxns)/length(model_save.rxns)*100)),'%%\n\t\t\t',...
                    'CF metabolites: ',num2str(round(length(model.mets)/length(model_save.mets)*100)),'%%\n'])

        prob.opt_fitFun     = opt_fitFun;
        prob.opt_rMMA       = opt_rMMA;
                
elseif compress
    % load compressed model
    model               = prob.model_compress;
    model_i             = prob.model_i_compress;
    opt_fitFun_save     = prob.opt_fitFun_save;
    opt_rMMA            = prob.opt_rMMA;
    opt_fitFun          = prob.opt_fitFun;
else
    % no compression, add dummy transformation vectors and matrices
    model.comprMapVec   = [1:length(model.rxns)]';
    model.comprMapMat   = speye(length(model.rxns));
end


    
%% reduce target space
if redFlag
    c = fix(clock);
    disp([num2str(c(4:end)),': Basic reduction of target space ...'])

    % exclude by reaction identifier
    [nonTargetRed]         = redTargetSpace(model, [], modelType);

    % exclude reaction that cannot carry any flux
    % determine blocked reactions
    [minFlux,maxFlux,~,~] = manualFVA(model);    % flux variability analysis
    blckd_rxns      = find((minFlux==0) & (maxFlux==0));
    % add to non targets
    nonTargetRed    = unique([nonTargetRed;blckd_rxns]);       
end
nonTarget   = [nonTarget;nonTargetRed];



%% Conduct rMMA and specify targets

% exclude heterologous insertion reaction 
nonTarget   = [nonTarget;model.hetRxnNum];

% setup target struct
if isfield(prob,'targets')
    targets     = prob.targets;
else
    tar_rxnNum              = 1:length(model.rxns);
    tar_rxnNum              = tar_rxnNum';
    tar_rxnNum(nonTarget)   = [];
    targets.rxnNum          = tar_rxnNum;
    targets.score           = ones(length(tar_rxnNum),1)./length(tar_rxnNum);   % assign equal probability to all targets
    targets.bound           = zeros(length(tar_rxnNum),2);                      % consider KO targets only
    targets.map             = eye(length(tar_rxnNum));
end

% add heterologous reaction insertion targets WRITE A "LOAD MODEL" PART
numHetRxns          = length(model.hetRxnNum);    % number of heterologous insertion reactions
targets.rxnNum_hri  = model.hetRxnNum;
targets.score_hri   = ones(numHetRxns,1);
targets.bound_hri   = [model.lb(model.hetRxnNum),model.ub(model.hetRxnNum)];
targets.map_hri     = eye(numHetRxns);
targets.Nt_hri      = numHetRxns;

% update target structure
targets.K       = K;
targets.K_hri   = K_hri;

% update reference flux distribution in targets struct
targets.fd_ref  = model.fd_ref;
% save for future (re-)runs
prob.targets    = targets;

c = fix(clock);
disp([num2str(c(4:end)),': Number target reactions: ',num2str(length(targets.rxnNum))])



%% If intervention size equals one do a brute force search
if K==1
    
    
    
end


%% Build initial population
c = fix(clock);
disp([num2str(c(4:end)),': Initialize genetic optimization algorithm ...'])
if ~isfield(prob,'pop') && ~isfield(prob,'enc')
    % initializePopulation
    [pop,pop_Tbin, pop_bin,enc,targets]     = initializePopulation(model,targets,popSize,...
                                                b,K,K_hri,optFocus,threads,initPopType,opt_ga);
    
    prob.enc    = enc;
    % save initial population
    prob.pop_init       = pop;  
else
    % load population
    pop     = prob.pop;
    enc     = prob.enc;
    if ~isfield(prob,'pop_bin')
        pop_bin = encode(pop,enc.numBits,enc.encodeVec_midPos,enc.encodeVec_midPos_hri,enc.K,enc.K_hri,enc.numRxn_KO);
    else
        pop_bin     = prob.pop_bin;
    end
    if ~isfield(prob,'pop_bin')
        pop_Tbin    = zeros(size(pop,1),length(targets.rxnNum));
        for i=1:size(pop,1)
            pop_Tbin(i,pop(i,:))    = 1;
        end
    else
        pop_Tbin    = prob.pop_Tbin;
    end
    % initialize additional m-files
    if strcmp(optFocus,'gene')
        % write m-File for gene logic evaluation in additional folder
        targets     = writeLogGeneRules(model,targets,opt_ga);
    elseif strcmp(optFocus,'rxns')
        % write m-File to get target numbers and bounds
        writeRxnRules(opt_ga);
    end
end
  



%% Initialize fitness function and storage
if fitFun_chngd
    % change in fitness function, recalculate fitness for the initial
    % population
    [popFit,popFitObj,chr_map,targets,model,model_i,opt_fitFun]   ...
        = initializeFitFun(model,targets,pop,pop_Tbin,fitFun_type,opt_fitFun,optFocus);
    % initializeFitFun
else
    chr_map     = prob.chr_map;
    popFit      = prob.popFit;
    popFitObj   = prob.popObjFit;
    opt_fitFun  = prob.opt_fitFun;
    opt_fitFun.excl_rxns_i  = prob.opt_fitFun.excl_rxns_i;
end
 


% save for future runs
prob.model          = model;
prob.model_i        = model_i;
prob.targets        = targets;
prob.opt_fitFun     = opt_fitFun;
prob.popFit_init    = popFit;
prob.popObjFit_init = popFitObj;

% save progress
save('Rescue\GAMO_PreGA_SaveFile')

%% Genetic algorithm
c = fix(clock);
disp([num2str(c(4:end)),': Execute genetic algorithm (best initial fitness: ',num2str(max(popFit)),') ...'])


[finalPop,chr_map,totalData] = geneticAlgorithm(model,model_i,targets,pop,pop_Tbin,pop_bin,popFit,enc,...
                              chr_map,fitFun_type,opt_ga,opt_fitFun);
% geneticAlgorithm

prob.pop            = finalPop.pop;
prob.pop_bin        = finalPop.pop_bin;
prob.pop_Tbin       = finalPop.pop_Tbin;
prob.popFit         = finalPop.popFit;
prob.popObjFit      = finalPop.popObjFit;
prob.popObjVal      = finalPop.popObjVal;
prob.chr_map        = chr_map;

% add total data to previous data
if ~isempty(pop_gen_tot) && ~isempty(popFit_gen_tot)
    [gen_prev,threads_prev]     = size(pop_gen_tot);
    if threads_prev<threads
        totalData.pop       = [pop_gen_tot,cell(gen_prev,threads-threads_prev);totalData.pop];
        totalData.popFit    = [popFit_gen_tot,cell(gen_prev,threads-threads_prev);totalData.popFit];
        totalData.popFitObj = [popFitObj_gen_tot,cell(gen_prev,threads-threads_prev);totalData.popFitObj];
    else
        totalData.pop       = [pop_gen_tot;totalData.pop,cell(opt_ga.maxGen,threads_prev-threads)];
        totalData.popFit    = [popFit_gen_tot;totalData.popFit,cell(opt_ga.maxGen,threads_prev-threads)];
        totalData.popFitObj = [popFitObj_gen_tot;totalData.popFitObj,cell(opt_ga.maxGen,threads_prev-threads)];
    end
    totalData.timing_drifts     = [timing_drifts_tot;totalData.timing_drifts];
    totalData.timing_gen        = [timing_gen_tot;totalData.timing_gen];

end
        
prob.totalData      = totalData;

% save progress
save('Rescue\GAMO_Final_SaveFile')

%% Postprocess results
c = fix(clock);
disp([num2str(c(4:end)),': Postprocess results ...'])


% decompress results
if compress
    % save compressed model
    prob.model_compress     = model;
    prob.model_i_compress   = model_i;
    prob.model              = model_save;
    % deactivate insertion target reactions
    model_save.lb(model_save.hetRxnNum)   = 0;
    model_save.ub(model_save.hetRxnNum)   = 0;
        
end


results     = [];
% sort population
[popFit_sort,popFit_sort_I]     = sort(finalPop.popFit,'descend');
popFit_sort(popFit_sort<0)      = 0;
pop_sort                        = finalPop.pop(popFit_sort_I,:);

% assign target type
numIS       = length(popFit_sort);
res         = [];
parfor i=1:numIS   
    
    indv                            = pop_sort(i,:);
    [KO,KD,Gene_KO,Gene_KD,Ins]     = translatePop(model_save,model,indv,optFocus,targets,enc.K,enc.K_hri);
    
    res(i).KO   = KO;
    res(i).KD   = KD;
    res(i).Gene_KD  = unique(Gene_KD);
    res(i).Gene_KO  = unique(Gene_KO);
    res(i).Ins      = Ins;
    res(i).fitness  = popFit_sort(i);
    res(i).pop      = pop_sort(i,:);
    
end


% get best solutions from map
res_best    = [];
allKeys     = keys(chr_map);
numBest     = length(allKeys);
allFit_cell = values(chr_map,allKeys);
if ~isempty(allFit_cell)
    allFit      = zeros(length(allFit_cell),length(allFit_cell{1}));
else
    allFit      = [];
end
for i=1:length(allFit_cell)
    allFit(i,:)   = allFit_cell{i};
end

[allFit_sort_I,allFit_sort_P]   = sort(allFit(:,1),'descend');
if length(allFit_sort_P)<numBest
    numBest     = length(allFit_sort_P);
end

parfor i=1:numBest

    indv                            = str2num(allKeys{allFit_sort_P(i)});   % get individual  
    [KO,KD,Gene_KO,Gene_KD,Ins]     = translatePop(model_save,model,indv,optFocus,targets,enc.K,enc.K_hri);
    
    % save results
    res_best(i).KO          = KO;
    res_best(i).KD          = KD;
    res_best(i).Ins         = Ins;
    res_best(i).Gene_KD     = Gene_KD;
    res_best(i).Gene_KO     = Gene_KO;
    res_best(i).fitness     = allFit_sort_I(i);
    res_best(i).pop         = indv;
    if fitFun_type==1
        res_best(i).objVal      = allFit(allFit_sort_P(i),2:end);
    end
end
i           = [];
Gene_KO     = [];
Gene_KD     = [];
 

% calculate flux distribution for the best 50 solutions
opt_fitFun.minInt   = 0;    % do not account for intervention size
numBS   = length(res_best);
if numBS>50
    numBS   = 50;
end
fluxes      = cell(numBS,1);
popObjVal   = cell(numBS,i);
numRxns_orig    = size(model_save.rxns,1);
parfor i=1:numBS
    [~,~,popFD,popObjVal_out] = evalFitness(res_best(i).pop,[],model,model_i,1,targets,fitFun_type,opt_fitFun);
    % decompress results
    if ~isempty(popFD{1})
        fluxes{i}                       = zeros(numRxns_orig,1);
        fluxes{i}(model.comprMapVec)    = popFD{1}; 
    else
        fluxes{i}   = [];
    end
    popObjVal{i}                    = popObjVal_out{1};
end
i   = [];

for j=1:numBS
    res_best(j).fluxes      = fluxes{j};
    res_best(j).popObjVal   = popObjVal{j};
end

% calculate flux distribution for the best 50 solutions in the final
% population
if numIS>50
    numIS   = 50;
end
fluxes      = cell(numIS,1);
popObjVal   = cell(numIS,i);

parfor i=1:numIS
    [~,~,popFD,popObjVal_out] = evalFitness(pop_sort(i,:),[],model,model_i,1,targets,fitFun_type,opt_fitFun);
    % decompress results
    if ~isempty(popFD{1})
        fluxes{i}                       = zeros(numRxns_orig,1);
        fluxes{i}(model.comprMapVec)    = popFD{1}; 
    else
        fluxes{i}   = [];
    end
    popObjVal{i}                    = popObjVal_out{1};
end
for i=1:numIS
    res(i).fluxes      = fluxes{i};
    res(i).popObjVal   = popObjVal{i};
end

% Test combinations of target reactions of the best solutions for fitness
best_sol    = 10;   % size of best solutions
if length(res_best)<best_sol
    best_sol    = length(res_best);
end

K_tot       = enc.K_tot;
K           = enc.K;
K_hri       = enc.K_hri;
res_comb    = cell(best_sol,1);
parfor i=1:best_sol
    comb        = [];
    comb.KO     = [];
    comb.KD     = [];
    comb.fit    = [];
    c   = 1;    % counter
    for j=1:(K_tot-1)
        % determine insertion targets
        if K_hri>0
            targets_hri     = res_best(i).pop((K+1):K_tot);
        else
            targets_hri  = [];
        end
        % create combinations of targets
        pop_comb    = nchoosek(res_best(i).pop,j); 
        targets_h   = targets;  % copy targets structure          
        for l=1:size(pop_comb,1)

            indv    = pop_comb(l,:);
            % determine number of insertions targets
            K_hri_c             = sum(ismember(indv,targets_hri));   
            K_c                 = length(indv)-K_hri_c;
            targets_h.K         = K_c;
            targets_h.K_hri     = K_hri_c;
         
            
            [KO,KD,~,~,Ins]     = translatePop(model_save,model,indv,optFocus,targets_h,K_c,K_hri_c);
            
            % calculate fitness (use alternative function)
            [popFit,~,~]    = evalFitness_comb(indv,model,model_i,1,targets_h,fitFun_type,opt_fitFun);
            
            comb(c).KO  = KO;
            comb(c).KD  = KD;
            comb(c).Ins = Ins;
            comb(c).fit = popFit;
            comb(c).pop = pop_comb(l,:);
            c   = c+1;     
        end
    end 
    res_comb{i}     = comb;
end




% save results
results.res_comb    = res_comb;
results.res         = res;
results.res_best    = res_best;
% results.chr_map     = chr_map;
% results.chr_map_obj = chr_map_obj;
% results.pop         = pop_sort;
% results.fitness     = popFit_sort;
% results.targets     = targets;
% results.genFit      = finalPop.genFit;
% results.model       = model;
results.timing      = finalPop.timing;
results.totalData   = totalData;
% results.opt         = opt;
% results.opt_fitFun  = opt_fitFun;

if ~isempty(saveFile)
    if isempty(saveFolder)
        save(saveFile,'results','prob')
    else
        if ~isdir(saveFolder)
            mkdir(saveFolder)
        end
        save([saveFolder,'/',saveFile],'results','prob')
    end
end

%% clean up workspace/folder
rmpath(opt_ga.AddFilesDir);

c = fix(clock);
disp([num2str(c(4:end)),': ... Done!'])


end