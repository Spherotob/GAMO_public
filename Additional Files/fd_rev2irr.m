function fd_i = fd_rev2irr(model,model_i,fd)
% transform reversible flux distribution to irreversible fd

fd_i    = zeros(length(model_i.rxns),1);
for i=1:length(fd)
    if length(model_i.rev2irrev{i}) > 1
        % reversible reaction
        if fd(i) > 0
            % positive forward flux of reversible reaction
            fd_i(ismember(model_i.rxns,[model.rxns{i},'_f']))   = fd(i);
        elseif fd(i) < 0
            % positive backward flux of reversible reaction
            fd_i(ismember(model_i.rxns,[model.rxns{i},'_b']))   = -fd(i);
        end
    else
        % irreversible reaction
        fd_i(model_i.rev2irrev{i})  = fd(i);
    end
end

end