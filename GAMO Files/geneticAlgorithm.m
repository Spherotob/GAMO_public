function [finalPop,chr_map,totalData] = geneticAlgorithm(model,model_i,targets,pop,pop_Tbin,pop_bin,popFit,enc,chr_map,fitFun_type,opt_ga,opt_fitFun)
% Genetic Algorithm for GAMO

%% INPUT
% model:            Stoichiometric metabolic model
% model_i:          Irreversible stoichiometric metabolic model
% targets:          Struct containing all possible targets
% pop:              Population containing intervention strategies as
%                   chromosomes
% pop_Tbin:         Population in binary intervention format
% pop_bin:          Population in binary format 
% popFit:           Fitness of each chromosome in the population
% enc:              Struct containing encoding information for the
%                   conversion from and to the binary format
% chr_map:          Container map holding fitness of all unique
%                   chromosomes ever created. Keys are strings of the
%                   binary intervention format
% fitFun_type:      Type of fitness function
% opt_ga:              Mandatory and optional options for the genetic
%                   algorithm
% opt_fitFun:       fitness function options

%% OUTPUT
% finalPop:         Struct containg the results of the genetic optimization
% chr_map:          Container map including fitness of all evaluated
%                   chromosomes

gaStart     = tic;

analyzePop  = 0;
genPeriod   = 1;

%% General parameter

Nt              = targets.Nt;           % number of possible targets

memPop          = opt_ga.memPop;        % Memorize chromosomes and their related fitness
genSize         = opt_ga.genSize;      % Number of generations which are iterated until parallel strands are merged
maxGen          = opt_ga.maxGen;       % Maximal number of sequential generations
threads         = opt_ga.threads;
% threads         = 2;                    % for development
X               = opt_ga.slctRate;      % selection rate
noMatingIdent   = opt_ga.noMatingIdent;
numKntChr       = opt_ga.numKntChr;     % Number of kinetochores (max. 2)
crossType       = opt_ga.crossType;     % Type of crossover routine
mutRate         = opt_ga.mutRate;       % mutation rate
elite           = opt_ga.elite;         % Number of elite solutions


%% start genetic algorithm
% prepare split complete population to several subpopulations (islands, demes)
[Np,K_tot]  = size(pop);
Np_s    = Np/threads;   % number of chromosomes in each subpopulations
popIdx  = 1:threads;
popIdx  = ((popIdx-1).*Np_s)+1;

% define general properties
Nslct       = ceil(X*Np_s);     % Number of chromosomes to be kept during selection
Ndel        = Np_s-Nslct;       % Number of chromosomes to be dismissed during selection
opt_fitFun.K_tot    = K_tot;    % maximal number of interventions    
% ensure even number of deleted chromosomes during selection
if mod(Ndel,2)~=0
    Nslct       = Nslct-1;
    Ndel        = Np_s-Nslct;
end
numPairs    = Ndel/2;       % Number of mating pairs (2 offsprings per pair)

% define factor for hamming distance calculations
hd_fac          = size(pop_Tbin,2)/((Np_s/2)*(Np_s-1)*K_tot*2);
hd_first        = sum(pdist(double(rand(Np_s,size(pop_Tbin,2))<(K_tot/size(pop_Tbin,2))),'hamming'))...
                    *hd_fac;
hd_fac          = hd_fac/hd_first;
opt_mut.hd_fac  = hd_fac;          % hamming distance factor


%
pop_tot         = cell(threads,1);  % struct to save final population of each workers
gen_fit         = cell(maxGen,1);   % store fitness of chromosomes for each generation
timing.drifts   = zeros(maxGen,1);  % Computation time of each drift
timing.gen      = cell(maxGen,1);  % Computation time of each generation
timing.fit      = cell(maxGen,1);   % Computation time of fitness evaluation of a whole population

pop_gen_tot         = cell(maxGen,threads);
popFit_gen_tot      = cell(maxGen,threads);
popFitObj_gen_tot   = cell(maxGen,threads);


genNum      = 1;
genFlag     = 1;
while (genNum<=maxGen) && genFlag
    genStart    = tic;  % start generation timer
    % display
    c = fix(clock);
    fprintf([num2str(c(4:end)),': Gendrift ',num2str(genNum),' of ',num2str(maxGen)]);
    
    % setup population arrays
%     pop_temp        = zeros(size(pop));
%     pop_Tbin_temp   = zeros(size(pop_Tbin));
%     pop_bin_temp    = zeros(size(pop_bin));
%     popFit_temp     = zeros(size(popFit));
        
    % distribute best strategies to different workers
    [~,popSort_idx]     = sort(popFit,'descend');
    shuffle             = randperm(Np);
    pop_temp      = pop(popSort_idx(shuffle),:);
    pop_Tbin_temp = pop_Tbin(popSort_idx(shuffle),:);
    pop_bin_temp  = pop_bin(popSort_idx(shuffle),:);
    popFit_temp   = popFit(popSort_idx(shuffle),:);
        
    % assign new population        
    pop         = pop_temp;
    pop_Tbin    = pop_Tbin_temp;
    pop_bin     = pop_bin_temp;
    popFit      = popFit_temp;
          
    % initialize parallel populations
    % Split population
    pop_dist        = cell(threads,1);
    pop_Tbin_dist   = pop_dist;
    pop_bin_dist    = pop_dist;
    popFit_dist     = pop_dist;
    for p=1:threads
        popStart        = popIdx(p);
        popStop         = popIdx(p)+Np_s-1;
        pop_dist{p}     = pop(popStart:popStop,:);
        pop_Tbin_dist{p}= pop_Tbin(popStart:popStop,:);
        pop_bin_dist{p} = pop_bin(popStart:popStop,:);
        popFit_dist{p}  = popFit(popStart:popStop,:);
    end
   
    
    parfor p=1:threads        
        % split population
        pop_s       = pop_dist{p};
        pop_Tbin_s  = pop_Tbin_dist{p};
        pop_bin_s   = pop_bin_dist{p};
        popFit_s    = popFit_dist{p};
        popFitObj_s = [];
        popObjVal_s = [];
        % assign container map
        if memPop
            % allocate container map
            chr_map_s       = chr_map;
        else
            % create empty container file
            chr_map_s   = containers.Map;
        end
        
%         % assign parameters
%         [~,hd_prev]     = adaptMutRate(pop_Tbin_s,0,0,opt_mut);
%         maxFit_prev     = max(popFit_s);    % maximal fitness of previous generation
%         staticGen_c_s   = round(mean(staticGen_c));    % static fitness generation counter
        
        % start sequential generation loop
        pop_gen_s       = cell(genSize,1);
        popFit_gen_s    = cell(genSize,1);
        popFitObj_gen_s = cell(genSize,1);   
        gen_timing      = zeros(genSize,1);
%         fit_timing      = zeros(genSize,1);
        for g=1:genSize
            t = tic; 
            
            % Selection
            [matingMat,pop_slct] = selection(pop_s,pop_Tbin_s,pop_bin_s,popFit_s,...
                                            Nslct,Ndel,numPairs,noMatingIdent);

            % Mating
            offspring_bin   = mating(pop_bin_s,K_tot,matingMat,numPairs,numKntChr,crossType,enc);
          
            % merge new population
            pop_bin_s_new   = [pop_bin_s(pop_slct,:);offspring_bin];
            
%             % adapt mutation rate
%             popFit_s_max            = max(popFit_s);
%             staticGen_c_s           = (staticGen_c_s+1)*(popFit_s_max<=maxFit_prev);
%             [mutRate_adapt,hd_prev] = adaptMutRate(pop_Tbin_s,hd_prev,staticGen_c_s,opt_mut);
%             
% %             disp([p,mutRate_adapt,hd_prev,staticGen_c_s])
% 
%             mutRate     = mutRate_adapt;
            
            % Mutation
            popFit_slct     = popFit_s(pop_slct);    % fitness of selected individuals
            [pop_s_new,pop_Tbin_s_new,pop_bin_s_new]   ...
                = mutation(pop_bin_s_new,popFit_slct,Np_s,Nt,K_tot,Nslct,enc,mutRate,elite,targets,chr_map_s,opt_mut);

            % calculate fitness of new population
%             tFit    = tic;
%             
           
            if memPop
                [popFit_s_new,popFitObj_s_new,~,popObjVal_s,chr_map_s]  ...
                            = evalFitness_mem(pop_s_new,pop_Tbin_s_new,chr_map_s,...
                                            model,model_i,Np_s,targets,fitFun_type,opt_fitFun);
            else
                [popFit_s_new,popFitObj_s_new,~,popObjVal_s]  ...
                            = evalFitness(pop_s_new,pop_Tbin_s_new,model,model_i,Np_s,targets,fitFun_type,opt_fitFun);
            end

%             fit_timing(g)  = toc(tFit);
            
            gen_timing(g)  = toc(t);
            
            
            % prepare for next generation
            pop_s       = pop_s_new;
            pop_Tbin_s  = pop_Tbin_s_new;
            pop_bin_s   = pop_bin_s_new;
            popFit_s    = popFit_s_new;
            popFitObj_s = popFitObj_s_new;
            % save generation
            pop_gen_s{g}        = pop_s_new;
            popFit_gen_s{g}     = popFit_s_new;
            popFitObj_gen_s{g}  = popFitObj_s_new;
           
            
        end
        % save last population for merging
        pop_struct              = [];
        pop_struct.pop          = pop_s;
        pop_struct.pop_Tbin     = pop_Tbin_s;
        pop_struct.pop_bin      = pop_bin_s;
        pop_struct.popFit       = popFitObj_s;      % save non-transformed fitnesses and transform with merged population
        pop_struct.chr_map      = chr_map_s;
        pop_struct.popObjVal    = popObjVal_s;
        pop_struct.gen_timing   = gen_timing;
%         pop_struct.fit_timing   = fit_timing;
%         pop_struct.staticGen_c   = staticGen_c_s;   % counter for subsequent generations without change in maximal fitness

        pop_tot{p}  = pop_struct;
        pop_gen_tot{genNum,p}       = pop_gen_s;
        popFit_gen_tot{genNum,p}    = popFit_gen_s;
        % save only if fitness was additionally transformed
        if opt_fitFun.minInt==1
            popFitObj_gen_tot{genNum,p} = popFitObj_gen_s;
        end
       
    end
    
    % merge populations
    pop_merge           = [];
    pop_Tbin_merge      = [];
    pop_bin_merge       = [];
    popFit_merge        = [];
    chr_map_merge       = pop_tot{1}.chr_map;
    popObjVal_merge     = [];
    gen_timing          = zeros(genSize,threads);
%     fit_timing          = zeros(genSize,threads);
    for i=1:threads
        pop_merge           = [pop_merge;pop_tot{i}.pop];
        pop_Tbin_merge      = [pop_Tbin_merge;pop_tot{i}.pop_Tbin];
        pop_bin_merge       = [pop_bin_merge;pop_tot{i}.pop_bin];
        popFit_merge        = [popFit_merge;pop_tot{i}.popFit];     % contains non-transformed fitnesses
        popObjVal_merge     = [popObjVal_merge;pop_tot{i}.popObjVal];
        chr_map_merge       = [chr_map_merge;pop_tot{i}.chr_map];
        gen_timing(:,i)     = pop_tot{i}.gen_timing;
%         staticGen_c(i)      = pop_tot{i}.staticGen_c;
%         fit_timing(:,i)     = pop_tot{i}.fit_timing;
    end
    
    % transform merged population
    popObjFit_merge     = popFit_merge;
    if opt_fitFun.minInt
        popFit_merge    = transformFitness(popFit_merge,pop_Tbin_merge,opt_fitFun);
    end
    
    
    
    
    % randomly shuffle complete population
    shuffle     = randperm(Np);     % get random permutations of vector indices
    pop         = pop_merge(shuffle,:);
    pop_Tbin    = pop_Tbin_merge(shuffle,:);
    pop_bin     = pop_bin_merge(shuffle,:);
    popFit      = popFit_merge(shuffle,:);
    popObjVal   = popObjVal_merge(shuffle,:);
    popObjFit   = popObjFit_merge(shuffle,:);

    % merge fitness container map
    chr_map     = chr_map_merge; 

%     a   = struct(chr_map);
%     b   = struct(chr_map_obj);
%     whos('a','b')
    
    % store fitness of current generation
    gen_fit{genNum}     = sort(popFit_merge,'descend');
   
    
    % save compuation time 
    timing.drifts(genNum,1)     = toc(genStart);
    timing.gen{genNum}          = gen_timing;
%     timing.fit{genNum}          = fit_timing;
    
    % display
    fprintf([' completed (best fitness: ',num2str(max(popFit)),')\n']);
    
    genNum  = genNum+1;
    
    % Dummy save select
%     save('Rescue\GAMO_SaveFile','pop_gen_tot','popFit_gen_tot','pop','pop_Tbin','pop_bin','popFit')
    % Dummy save all
    save('Rescue\GAMO_SaveFile_GA') 
    
    % clear space and reallocate variables
    clearvars chr_map_merge pop_tot
    pop_tot     = cell(threads,1);
    

%     whos
%     clearvars -except genNum pop pop_Tbin pop_bin popFit gen_fit timing gaStart ...
%                         pop_gen_tot popFit_gen_tot chr_map 
    
    
end

% merge results
finalPop.pop        = pop;
finalPop.pop_Tbin   = pop_Tbin;
finalPop.pop_bin    = pop_bin;
finalPop.popFit     = popFit;
finalPop.popObjFit  = popObjFit;
finalPop.popObjVal  = popObjVal;
finalPop.genFit     = gen_fit;
finalPop.timing     = timing;
finalPop.timing.tot = toc(gaStart);
% save all results
totalData.pop       = pop_gen_tot;
totalData.popFit    = popFit_gen_tot;
totalData.timing_drifts = timing.drifts;
totalData.timing_gen    = timing.gen;
if opt_fitFun.minInt
    totalData.popFitObj     = popFitObj_gen_tot;
else
    totalData.popFitObj     = [];
end