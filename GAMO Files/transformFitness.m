function popFitTrans = transformFitness(popFit,pop_Tbin,opt_fitFun)
% transforms fitness values according to the effective number of
% interventions of individuals

% INPUT
% popFit:             Fitnesses of the population
% pop_Tbin:           Population in target-specific binary format
% opt_fitFun:         Options for fitness function

% OUTPUT
% popFitTrans:        Transformed fitness according to the number of interventions of the individuals


%%
% calculate effetive number of interventions for each individual
numIntEff   = sum(pop_Tbin,2);
% % determine best individual
% [~,posMaxFit]   = max(popFit);
% % determine effective intervention size of best individual
% numIntMaxFit    = numIntEff(posMaxFit);
% transform fitnesses
popFitTrans     = popFit+(popFit.*opt_fitFun.FIRF.*(opt_fitFun.numIntMaxFit-numIntEff));


end