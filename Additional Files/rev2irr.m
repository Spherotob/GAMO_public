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
    ub_subs     = model.ub(strcmp(model.rxns,model.subsRxn));
    lb_subs     = model.lb(strcmp(model.rxns,model.subsRxn));
    if (lb_subs<0) && (ub_subs<=0)
        % pure backward reaction
        warning('Consider to make substrate uptake reaction reversible to avoid problems with the calculation of the fitness!')
        model_i.subsRxn   = [model.subsRxn,'_r'];
    elseif (lb_subs<0) && (ub_subs>0)
        % reversible reaction, uptake is represented by a negative reaction rate
        model_i.subsRxn   = [model.subsRxn,'_b'];
    elseif (lb_subs>=0) && (ub_subs>0)
        % RARE CASE, uptake reaction is irreversible and is represented by
        % a forward flux
        warning('Substrate uptake reaction bounds are positive! Generally, uptake fluxes should be represented by a backward reaction')
        model_i.subsRxn   = model.subsRxn;
    end
    
end
if isfield(model,'targetRxn')
    ub_target     = model.ub(strcmp(model.rxns,model.targetRxn));
    lb_target     = model.lb(strcmp(model.rxns,model.targetRxn));
    if (lb_target<0) && (ub_target>0)
        % target reaction is reversible, treat forwards reaction as the
        % main objective
        model_i.targetRxn   = [model.targetRxn,'_f'];
    elseif (lb_target>=0) && (ub_target>0)
        % standard irreversible reaction
        model_i.targetRxn   = model.targetRxn;
    elseif (lb_target<0) && (ub_target<=0)
        % RARE CASE, target reaction backward flux is the objective
        warning('Target reaction bounds are negative! Generally, target fluxes should be represented by a forward reaction')

        model_i.targetRxn   = [model.targetRxn,'_r'];
    end
end

end