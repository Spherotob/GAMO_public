%% Genetic Algorithm for Metabolic Optimization (GAMO)
% Example for GAMO employing the E. coli iAF1260 core and GEM iJO1366 model
% (http://bigg.ucsd.edu/)

clearvars -except prob 

opt             = [];


%%
% set path for additional files
addpath('GAMO Files')
addpath('Data Files')            % flux data for reference flux distribution
addpath('Additional Files')  


%% load model
% iAF1260 core model
% load('iAF1260Core.mat')
load('iJO1366.mat')
model.csense    = model.csense';

model           = changeRxnBounds(model,'EX_glc__D_e',-13.34,'l');   % activate glucose uptake
model           = changeRxnBounds(model,'EX_o2_e',-20,'l');   % activate oxygen uptake

% define target (exchange) reaction, biomass formation reaction and
% substrate uptake reaction
model.targetRxn     = 'EX_succ_e';
% model.bmRxn         = 'BIOMASS_Ecoli_core_w_GAM';   % E. coli core model
model.bmRxn         = 'BIOMASS_Ec_iJO1366_WT_53p95M';   % E. coli GEM iJO1366
model.subsRxn       = 'EX_glc__D_e';

% Reference flux distribution (Ishii 2007, growth on glucose, mu=0.7 1/h)
opt.filename    = 'fluxData_Ishii_2007';
[model,~,res]   = createRefFD(model,[],1,opt);  % creates reference flux distribution (see field "fd_ref" in model struct)


% GAMO OPTIONS
%% general options
opt.saveFile    = 'iAF1260_Core_TEST';  % Filename for saving data
opt.saveFolder  = 'Results';            % folder for saving results
opt.threads     = 7;                    % Number of parallel workers


%% model preprocessing
opt.redFlag     = 1;    % (0): no target space reduction; (1): target space reduction routine
opt.compress    = 1;    % model compression (0): No; (1): active

%% General (crucial) optimization options
opt.numInt          = 5;        % Maximal number of interventions
opt.optFocus        = 'gene';   % 'rxns': reaction deletions; 'gene': gene deletions

%% genetic algorithm parameter (crucial)
opt.memPop          = 1;    % (1): Memorize fitness of each generated chromosome
                            % (0): only final population fitness is passed
                            
opt.popSize         = 20;   % Population size/Number of chromosomes per generation     
opt.maxGen          = 2;   % Numberof Gene-Drift-Events
opt.genSize         = 2;   % Number of generations between two Gene-Flow-Events 
opt.slctRate        = 0.25;  % Selection rate
opt.mutRate         = 0.05;  % Mutation rate (0-1) related to the whole population and its number of bits
opt.elite           = 1;    % Number of elite chromosomes which are not to be mutated

%% genetic algorithm parameter (additional)
opt.initPopType     = 0;    % Type of initial population
                            % (0): randomly select interventions
                            % (1): only single intervention mutants 

opt.slctPressure    = 2;    % Selection Pressure for targets at population initilization (1-2)


%% fitness function options
opt.fitFun                 = 0;    % Choose fitness function (0): MiMBL; (1): Multiobjective
% (G):  gcOpt (4)
% (R):  RobustKnock (3)
% (O):  Optknock (2)   
% (M):  MiMBl (1)
opt_fitFun.weighting       = [1,1,1,1];     % weighting factor for each objective function (see numbering above)
opt_fitFun.isObj           = 'OG';          % choose objectives to be evaluated (Include characters from above list in the string)
opt_fitFun.leadObj         = 'M';           % Choose one master objective, if this objective value 
                                            % is zero the fitness of the individual is zero
                                            % disregarding the other objectives
                                            
opt_fitFun.minGrowth       = 0.1;   % Assign minimal allowable growth rate for LP objectives (OptKnock etc)

%% exclude reactions from metabolite balance (basically all non-gene-related
% reactions). These will not be considered for the MiMBl optimization
excl_rxns   = find(any(model.rxnGeneMat,2)==0);
excl_rxns   = [excl_rxns; find(not(cellfun('isempty',(strfind(model.rxnNames,'diffusion')))))];
excl_rxns   = [excl_rxns; find(not(cellfun('isempty',(strfind(model.rxnNames,'spontaneous')))))];
opt_fitFun.excl_rxns    = unique(excl_rxns);


%% exclude reactions from target space
opt.nonTarget   = {};   % manually specify reaction identifiers of reactions not to be targeted


% ADDITIONAL FEATURES
%% heterologous reaction insertions
opt.numInsertions   = 0;    % Maximal number of heterologous reaction insertions
if opt.numInsertions > 0
    % if heterologous insertions are considered include novel reactions in
    % the model. Heterologous reactions were derived from MetaNetX databank
    % search.
    % load database model
    load('heterogeneous_rxn_list_iAF1260Core.mat')  % load list including novel network edges (list for genome-scale model iJO1366 also provided)
    model           = addNetworkBranches(model,rxn_list_trunc);     % add novel reaction to the model
    opt.hetRxnNum   = model.hetRxnNum;
    excl_rxns       = [opt_fitFun.excl_rxns;model.hetRxnNum];    % add reactions exclusion list
    opt_fitFun.excl_rxns    = unique(excl_rxns);
end

%% Minimization of intervention set sizes
opt_fitFun.minInt          = 0;     % (1): Minimize intervention set size during optmization
                                    % (0): Disable intervention set size minimization
opt_fitFun.FIRF            = 0.1;   % Fitness-Intervention Relation Factor


%% execute GAMO
start   = tic;

% execute GAMO
prob    = [];
[results,prob] = GAMO(model, opt, opt_fitFun, prob);

toc(start)
%% INFORMATION:

% The results struct contains the following relevant fields
% res:        contains positions of the knockouts (KO) and insertions (Ins) (see model struct),
%             fitness values and flux distributions of each individual of the final population
% res_best:   same as field "res" but including all individuals which emerged 
%             during the genetic algorithm run (only if opt.memPop=1)
% res_comb:   fitness values of all combinations of targets for the best 10 final individuals
% timing:     timing data of all the conducted Gene-Flow-Events and generations
% totalData:  Fitness values and targets of each individual of every population

% To restart GAMO using the final population as the initial population load
% the prob struct and parse it to the GAMO function. Via the opt struct,
% changes in the number of generations or Gene-Flow-Events can be made.


