function model_i = rev2irr(model)
% convert reversible COBRA model to irreversible model

[model_i,matchRev,rev2irrev_map,irrev2rev_map]  = convertToIrreversible(model);

N   = length(irrev2rev_map);    % number reactions in irreversible model
R   = length(rev2irrev_map);    % number reactions in reversible model

% save mapping vectors in model struct
model_i.matchRev        = matchRev;
model_i.rev2irrev       = rev2irrev_map;
model_i.irrev2rev       = irrev2rev_map;

% Create mapping matrix B
B   = zeros(N,R);
for i=1:N
    B(i,irrev2rev_map(i))   = 1;
end
model_i.B   = B;

% create mapping matrix for transforming irreversible flux distreibution to
% reversible
mapIrr2Rev  = zeros(R,N);
for i=1:R
    irrRxnNum   = rev2irrev_map{i};
    if length(irrRxnNum)>1
        mapIrr2Rev(i,ismember(model_i.rxns,[model.rxns{i},'_f']))  = 1;
        mapIrr2Rev(i,ismember(model_i.rxns,[model.rxns{i},'_b']))  = -1;
    else
        mapIrr2Rev(i,ismember(model_i.rxns,model.rxns{i}))  = 1;
    end
end
model_i.mapIrr2Rev  = mapIrr2Rev;
        

% append information from reversible model
if isfield(model,'bmRxn')
    model_i.bmRxn   = model.bmRxn;
end
if isfield(model,'subsRxn')
    % substrate uptake reaction points to extracellular space
    model_i.subsRxn   = [model.subsRxn,'_b'];
    
end
if isfield(model,'targetRxn')
    model_i.targetRxn   = model.targetRxn;
end

end