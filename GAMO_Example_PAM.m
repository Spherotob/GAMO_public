%% Genetic Algorithm for Metabolic Optimization (GAMO)
% Example for GAMO employing the E. coli iML1515 (http://bigg.ucsd.edu/) protein allocation model

clearvars -except prob 

opt             = [];


%%
% set path for additional files
addpath('GAMO Files')
addpath('Data Files')            % flux data for reference flux distribution
addpath('Additional Files')  


%% load model
load('C:\Users\altert\Promo\Matlab Code\Protein Allocation\AA independent model\iML1515_PAM_BRENDA_EnzymeMassBalance_MeasuredProtConc_kcatMean20190205_UnUnderUitilizedProtein_20190315_SubsDepend.mat')
model_i_pa.csense   = model_i_pa.csense';
model   = model_i_pa;
% change pam parameter
params.ms_max   = 11.1;
model   = changePAAModelParameter(model,params);

model           = changeRxnBounds(model,'EX_glc__D_e_b',11.1,'u');   % activate glucose uptake
model           = changeRxnBounds(model,'EX_o2_e_b',100,'u');   % activate oxygen uptake

% define target (exchange) reaction, biomass formation reaction and
% substrate uptake reaction
model.targetRxn     = 'EX_succ_e_f';
model.bmRxn         = 'BIOMASS_Ec_iML1515_WT_75p37M';   % E. coli iML1515
model.subsRxn       = 'EX_glc__D_e_b';


% GAMO OPTIONS
%% general options
opt.saveFile    = 'iML1515_pa_TEST';  % Filename for saving data
opt.saveFolder  = 'Results';            % folder for saving results
opt.threads     = 7;                    % Number of parallel workers


%% model preprocessing
opt.redFlag     = 1;    % (0): no target space reduction; (1): target space reduction routine
opt.compress    = 0;    % model compression (0): No; (1): active

%% General (crucial) optimization options
opt.numInt          = 3;        % Maximal number of interventions
opt.optFocus        = 'rxns';   % 'rxns': reaction deletions; 'gene': gene deletions

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
opt.fitFun                 = 2;    % Choose fitness function (0): MiMBL; (1): Multiobjective; (2): Protein allocation model optimization

% Choose fitness parameter to be calculated from flux distribution                                           
opt_fitFun.fitParam        = 0;             % 0: Biomass-Product coupled yield (BPCY) [mol/mol/h]
                                            % 1: Yield [mol/mol]
                                            % 2: Production rate [mmol/g/h]
                                            
opt_fitFun.minGrowth       = 0.1;   % Assign minimal allowable growth rate for LP objectives (OptKnock etc)

% specify protein allocation options
opt_fitFun.subsGrid     = 3:0.01:11.1;  % substrate uptake rate grid
opt_fitFun.initGridNum  = 30;   % number of initial grid size
opt_fitFun.maxOEELimit  = 3.6;  % maximally allowable maxmimal enzyme overexpression capacity [mg/g/h]
% determine wild-type solutions for substrate uptake grid
wildtypeGridSolutions   = cell(length(opt_fitFun.subsGrid),1);
parfor i=1:length(opt_fitFun.subsGrid)
    globalSolverVariable('gurobi');
    model_m     = model;
    model_m     = changeRxnBounds(model_m,model_m.subsRxn,opt_fitFun.subsGrid(i),'u');
    sol         = optimizeCbModel(model_m,'max');
    wildtypeGridSolutions{i}     = sol.x;  
end
opt_fitFun.wildtypeGridSolutions    = wildtypeGridSolutions;

%% exclude reactions from metabolite balance (basically all non-gene-related
% reactions). These will not be considered for the MiMBl optimization
excl_rxns   = find(any(model.rxnGeneMat,2)==0);
excl_rxns   = [excl_rxns; find(not(cellfun('isempty',(strfind(model.rxnNames,'diffusion')))))];
excl_rxns   = [excl_rxns; find(not(cellfun('isempty',(strfind(model.rxnNames,'spontaneous')))))];
opt_fitFun.excl_rxns    = unique(excl_rxns);


%% exclude reactions from target space
opt.nonTarget   = {};   % manually specify reaction identifiers of reactions not to be targeted


% ADDITIONAL FEATURES
%% heterologous reaction insertions (not applicable for protein allocation models
opt.numInsertions   = 0;    % Maximal number of heterologous reaction insertions

%% Minimization of intervention set sizes (not applicable for protein allocation models)
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


