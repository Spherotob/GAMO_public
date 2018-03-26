function model_d = addNetworkBranches(model,reactionList)
%% Add network branches to wild-type model
% novel branches/reactions are selected from a survey of all reactions
% found at MetaNetX

% reaction list has to have the following structure
% Col 1: reaction formula in COBRA format
% Col 2: reaction identifier
% Col 3: reaction (full) name
% Col 4: specifier for main database (MetaNetX)
% Col 5: specifier for secondary database
% Col 6: reversibility flag
% Col 7: maximal Gibbs Free Energy of Reaction
% Col 8: predicted Gibbs Free Energy of Reaction for standard metabolite concentrations
% Col 9: standard Gibss Free Energy of Reaction

% function also predicts a reference flux for heterologous reactions

%%
% define parameters
lList       = length(reactionList);     % number of reactions to be added
model_d     = model;    % copy model

muFac           = [0.98,1];     % factor for maximal growth rate at which reference flux is calculated
fluxRate_std    = 2;    % standard reference flux rates for unknown flux directions
maxFluxRate     = 20;   % maximal flux rate
limGibbs        = 15;   % threshold for Gibbs Free Energy of Reaction to allow for maximal reaction velocity
refMu           = model.fd_ref(ismember(model.rxns,model.bmRxn));
subsUpRate      = model.fd_ref(ismember(model.rxns,model.subsRxn));

% possible warnings during reaction addition
warn    = {'not in model - added to the model',...
            'Model already has the same reaction you tried to add:',...
            'Reaction with the same name already exists in the model'};

% add reaction to model

rxn_valid   = zeros(lList,1);   % valid vector
refFlux     = zeros(lList,1);   % reference fluxes
for r=1:lList  
    lastwarn('') % clear last warning message
    model_p     = addReaction(model_d,reactionList{r,2},...
                                'reactionName',reactionList{r,3},...
                                'reactionFormula',reactionList{r,1},...
                                'subSystem','Heterologous reaction');    
    % check if warning has occured
    msg         = lastwarn;
    if isempty(msg)
        % reaction insertion was valid, continue
        rxn_valid(r)    = 1;
    else
        % some warnings occured
        for i=1:length(warn)
            % check if warning is critical
            if isempty(strfind(msg,warn{i}))
                rxn_valid(r)    = 1;
            else
                % discard reaction
                rxn_valid(r)    = 0;
                break
            end
        end
    end
    % check reaction consistency
    if ~any(model_p.S(:,end))
        % no metabolites participate in reaction!
        rxn_valid(r)    = 0;
    end
    % continue if reaction is not valid
    if rxn_valid(r)==0
        continue
    end
    
    % determine reference flux for new reaction
    model_t     = addReaction(model,reactionList{r,2},...
                                'reactionName',reactionList{r,3},...
                                'reactionFormula',reactionList{r,1},...
                                'subSystem','Heterologous reaction');   % set up testing model
    
    model_t     = changeRxnBounds(model_t,model_t.subsRxn,subsUpRate,'l');  % constrain substrate uptake rate according to reference                     
    % calculate maximal biomass formation                        
    model_t     = changeObjective(model_t,model_t.bmRxn);
    sol         = optimizeCbModel(model_t,'max','one');
    if ~strcmp(sol.origStat,'OPTIMAL')
        rxn_valid(r)    = 0;
        warning(['FBA failed considering reaction: ',reactionList{r,2}])
        continue
    end
    maxMu       = sol.f;
    flux_maxMu  = sol.x(end);
    
    % calculate lower and upper bound of heterologous reaction at maximal
    % and sub-maximal growth
    fluxRates_ub    = zeros(2,1);
    fluxRates_lb    = zeros(2,1);
    model_t     = changeObjective(model_t,model_t.rxns(end));
    for t=1:2
        model_t     = changeRxnBounds(model_t,model_t.bmRxn,maxMu*muFac(t),'b');
        sol         = optimizeCbModel(model_t,'max');
        fluxRates_ub(t)   = sol.f;
        sol         = optimizeCbModel(model_t,'min');
        fluxRates_lb(t)   = sol.f;
    end   
    fluxRates_ranges    = fluxRates_lb.*fluxRates_ub;
    
    % choose reference flux rate
    if flux_maxMu==0
        % consult flux rates bounds
        
        if ~any([fluxRates_ub;fluxRates_lb])
            % reaction cannot be active at (sub) maximal growth
            % exclude reaction since it may be too severe for the
            % metabolism
            warning(['Reaction: ',reactionList{r,2},' cannot carry any flux at (sub-)optimal growth'])
            rxn_valid(r)     = 0;
            
        elseif all(fluxRates_ranges<0)
            % reaction direction not determinable -> consult Gibbs
            % Free Energy
            if abs(reactionList{r,8})>limGibbs
                refFlux(r)  = fluxRate_std-((fluxRate_std*2)*(reactionList{r,8}>0));
            else
                refFlux(r)  = (fluxRate_std)*-reactionList{r,8}/limGibbs;
            end
            
        else
            for i=1:2
                if fluxRates_ub(i)>0 && abs(fluxRates_ub(i))<maxFluxRate
                    refFlux(r)  = (fluxRates_ub(i)/(maxMu*muFac(i)))*refMu;
                    break;

                elseif fluxRates_lb(i)<0 && abs(fluxRates_lb(i))<maxFluxRate
                    refFlux(r)  = (fluxRates_lb(2)/(maxMu*muFac(i)))*refMu;
                    break;
                end
            end
               
        end
        
        
    else
        if abs(flux_maxMu)<maxFluxRate
            refFlux(r)  = (flux_maxMu/maxMu)*refMu;
        else
            refFlux(r)  = fluxRate_std-((fluxRate_std*2)*(flux_maxMu>0));
        end
        
    end
    
    % check validity and continue
    if rxn_valid(r)
        % reaction and reference flux are valid
        model_d         = model_p;
    end
end

hetRxnNum   = zeros(lList,1);   % position of inserted reactions
hetRxns     = reactionList;     % information about the inserted reactions
fd_ref      = model_d.fd_ref;   % load reference flux distribution
for r=1:lList
    if rxn_valid(r)
        hetRxnNum(r)            = find(ismember(model_d.rxns,reactionList{r,2}));  
        fd_ref(hetRxnNum(r))    = refFlux(r);
    end
end
hetRxns(rxn_valid==0,:)       = [];
hetRxnNum(hetRxnNum==0)     = [];
model_d.hetRxns             = hetRxns;
model_d.hetRxnNum           = hetRxnNum;
model_d.fd_ref              = fd_ref;

% dismiss problems with subsystem annotation
for i=1:length(model_d.rxns)
    if isempty(model_d.subSystems{i})
        model_d.subSystems{i}  = 'NONE';
    end
end
model_d.subSystems  = [model_d.subSystems{:}]';

end