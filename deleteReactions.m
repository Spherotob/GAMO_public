function model = deleteReactions(model,del_rxns)
% deleting reactions from a stoichiometric model



%% delete reactions
% define reaction to be kept
keepRxns                = [1:length(model.rxns)]';
keepRxns(del_rxns)      = []; 

model.S   = model.S(:,keepRxns);    % stoichiometric matrix
model.rxnGeneMat  = model.rxnGeneMat(keepRxns,:);   % reaction gene matrix
model.grRules     = model.grRules(keepRxns);    % gene rules
model.rxns        = model.rxns(keepRxns);   % reaction identifiers
model.rxnNames    = model.rxnNames(keepRxns);   % reaction names
if isfield(model,'subSystems')
    model.subSystems  = model.subSystems(keepRxns);     % reaction compartments
end
model.lb  = model.lb(keepRxns);     % lower reaction bounds
model.ub  = model.ub(keepRxns);     % upper reaction bounds

model.c   = zeros(length(keepRxns),1);  % optimization objective
if isfield(model,'rev')
    model.rev     = model.rev(keepRxns);    % reaction reversibilities
end


%% delete metabolites
% define metabolites to be kept
delMet              = find(~any(model.S,2));
keepMet             = [1:length(model.mets)]';
keepMet(delMet)     = [];
model.S             = model.S(keepMet,:);   % stoichiometric matrix
model.mets          = model.mets(keepMet);  % metabolite identifiers
model.metNames      = model.metNames(keepMet);  % metabolite names
if isfield(model,'metFormulas')
    model.metFormulas     = model.metFormulas(keepMet);
end
if isfield(model,'csense')
    model.csense  = model.csense(keepMet);
end
if isfield(model,'b')
    model.b   = model.b(keepMet);
end

end