function [model,model_i,res] = createRefFD(model,model_i,src_flag,opt)
% calculate reference flux distribution using a stoichiometric and, if
% provided, experimental flux data

% INPUT
% model:    Stoichiometric model (COBRA format)
% model_i:  Irreversible stoichiometric model (optional)
% src_flag: Identifier for choosing reference flux calculation method
%             (0): Maximization of biomass reaction flux and minimization
%             of fluxes
%             (1): Applying experimental flux data as flux bounds and
%             minimizing fluxes (pFBA)
% opt:      additional options

% OUTPUT
% model:    Adapted stoichiometric model
% model_i:  Adapted ireversible stoichiometric model
% res:      Results structure

if nargin<4
    opt    = [];
end

%% Check options
if isfield(opt,'filename')
    filename    = opt.filename;
end
if isfield(opt,'fluxFac')
    fluxFac     = opt.fluxFac;
else
    fluxFac     = 1;
end

%% decides on which basis reference flux distribution is calculated
switch src_flag
    % Maximization of biomass and minimization of fluxes
    case 0
        model           = changeObjective(model,model.bmRxn);
        sol             = optimizeCbModel(model,'max','one');
        fd_ref          = sol.x;
        model.fd_ref    = fd_ref;
        res             = {};
    % Use pFBA (Lewis 2010, Long 2016) and experimental extracellular/intracellular fluxes    
    case 1  
        if isempty(filename)
            error('No experimental flux data provided for flux distribution calculation')
        end
        
        % load experimental data
        % Column 1: reaction identifiers
        % Column 2+3: Lower and upper bounds of reaction
        [num,txt,raw]   = xlsread(filename);
        % Check consistency
        if ~strcmp(raw{1,1},'Rxn Identifier') || ~strcmp(raw{1,2},'Mean') || ~strcmp(raw{1,3},'LB') || ~strcmp(raw{1,4},'UB') 
            error('Excel file with measured fluxes is not consistent')
        end
        % load data and apply factor for fluxes
        model_m     = model;
        rxns_exp    = txt(2:end,1);
        rxns_exp_lb     = zeros(length(rxns_exp),1);
        rxns_exp_ub     = rxns_exp_lb;
        rxns_exp_f      = rxns_exp_lb;
        for i=1:length(rxns_exp)
            if ~strcmp(rxns_exp{i},model.bmRxn)
                rxns_exp_lb(i) = num(i,2).*fluxFac;
                rxns_exp_ub(i) = num(i,3).*fluxFac;
                rxns_exp_f(i)  = num(i,1).*fluxFac;
            else
                rxns_exp_lb(i) = num(i,2);
                rxns_exp_ub(i) = num(i,3);
                rxns_exp_f(i)  = num(i,1);
            end
        end

        % General parameter
        numRxns     = length(model.rxns);
        numMets     = length(model.mets);
        numExpRxns  = length(rxns_exp);     % number of experimentally determined fluxes
        rxnNum_exp  = zeros(numExpRxns,1);
        for i=1:numExpRxns
            rxnNum_pre      = find(ismember(model_m.rxns,rxns_exp{i}));
            if ~isempty(rxnNum_pre)
                rxnNum_exp(i)  = rxnNum_pre;
            else
                % Check for biomass reaction
                if strcmp('BIOMASS',rxns_exp{i})
                    rxnNum_exp(i)   = find(ismember(model_m.rxns,model_m.bmRxn));
                    rxns_exp{i}     = model_m.bmRxn;
                end
            end
            
        end
        % erase reactions that are not in the model
        eraseRxn                = rxnNum_exp==0;
        disp(rxns_exp(eraseRxn))
        if any(eraseRxn)  
            rxnNum_exp(eraseRxn)    = [];
            rxns_exp                = rxns_exp(~eraseRxn);
            rxns_exp_lb(eraseRxn)   = [];
            rxns_exp_ub(eraseRxn)   = [];
            rxns_exp_f(eraseRxn)    = [];
            numExpRxns              = length(rxns_exp);
        end
        
        
        
        bm_constr_flag  = 0;    % flag for biomass formation constraint
        for i=1:length(rxns_exp)
            % constrain model (only exchange reactions)
            is_exch     = ~isempty(strfind(model.rxnNames{rxnNum_exp(i)},'exchange'));
            if is_exch
                model_m     = changeRxnBounds(model_m,rxns_exp{i},...
                                rxns_exp_lb(i),'l');
                model_m     = changeRxnBounds(model_m,rxns_exp{i},...
                                rxns_exp_ub(i),'u');
            end
            % fix biomass formation rate
            if strcmp(rxns_exp{i},model.bmRxn)
                bm_constr_flag  = 1;
                model_m     = changeRxnBounds(model_m,rxns_exp{i},...
                                rxns_exp_lb(i),'l');
                model_m     = changeRxnBounds(model_m,rxns_exp{i},...
                                rxns_exp_ub(i),'u');
            end
        end
        
        % setup flux difference LP
        numVar                  = numExpRxns+numRxns;
        % constraints for auxiliary variables (quadratic varaibles)
        add_constr              = zeros(numExpRxns,numVar);
        for i=1:numExpRxns
            add_constr(i,rxnNum_exp(i))     = 1;
            add_constr(i,numRxns+i)         = -1;
        end
    
        gur.A                   = [model_m.S,sparse(numMets,numExpRxns);...
                                    sparse(add_constr)];    % mass balance equations
        gur.obj                 = zeros(numVar,1);  % only quadratic objective function
        gur.Q                   = [sparse(numRxns,numVar);...
                                    sparse(numExpRxns,numRxns),speye(numExpRxns)];
        gur.sense(1:(numMets+numExpRxns),1)   = '=';
        gur.rhs                 = [zeros(numMets,1);... 
                                    rxns_exp_f];
        gur.lb                  = [model_m.lb;ones(numExpRxns,1).*-Inf];
        gur.ub                  = [model_m.ub;ones(numExpRxns,1).*Inf];
        gur.vtype(1:numVar,1)  = 'C';
        gur.modelsense          = 'min';
        % gurobi Parameter
        gurParams.OutputFlag    = 0;
        gurParams.Presolve      = 1;
        gurParams.Threads       = 1;
        gurParams.Method        = 1;    % use dual simplex for better convergense
        sol                     = gurobi(gur,gurParams);
        
        % save results
        res.sol_minFluxDiff     = sol;
        res.objval_fac          = sol.objval/fluxFac;
       
        
        rxnFlux_exp     = sol.x(rxnNum_exp);
        
        % constrain measured fluxes
        model_m     = changeRxnBounds(model_m,rxns_exp,rxnFlux_exp,'b');
        % if growth rate has not been fixed before, determine (maximal)
        % growth rate 
        if ~bm_constr_flag
            model_m     = changeObjective(model_m,model_m.bmRxn);
            opt         = optimizeCbModel(model_m,'max');
            bmFlux      = opt.x(strcmp(model_m.rxns,model_m.bmRxn));
            model_m     = changeRxnBounds(model_m,model_m.bmRxn,bmFlux,'l');          
        end
        
        % only consider gene related reactions for flux minimizations
        rxns_gene_flag              = any(model_m.rxnGeneMat,2);        
        rxns_gene                   = zeros(numRxns,1);
        rxns_gene(rxns_gene_flag)   = 1;
        
        % Setup minimization of absolute fluxes while constraining measured fluxes
        % and growth rate according to previous optimizations
        gurParams.Method        = -1;   % automatic choice
        
        gur                     = [];
        gur.A                   = [model_m.S,sparse(numMets,numRxns);...
                                    speye(numRxns),-speye(numRxns);...
                                    -speye(numRxns),-speye(numRxns)];
        gur.rhs                 = zeros(numMets+(2*numRxns),1);
        gur.sense(1:numMets,1)  = '=';
        gur.sense((numMets+1):(numMets+(2*numRxns)),1)  = '<';
        gur.obj                 = [zeros(numRxns,1);...
%                                     ones(numRxns,1)];
                                    rxns_gene];
        gur.modelsense          = 'min';
        gur.vtype(1:(2*numRxns),1)  = 'C';
        gur.lb                  = [model_m.lb;zeros(numRxns,1)];
        gur.ub                  = [model_m.ub;ones(numRxns,1)*Inf];
        sol2                    = gurobi(gur,gurParams);
        
        res.sol_minFluxNorm     = sol2;
        
        fd_ref                  = sol2.x(1:numRxns);
        
        model.fd_ref            = fd_ref;
        
        res.expFlux_mean        = rxns_exp_f;
        res.simFlux             = fd_ref(rxnNum_exp);
        res.diffFlux            = (res.expFlux_mean-res.simFlux)./res.expFlux_mean;
        res.expRxn              = rxns_exp;
%         save TEST
end

%% translate to irreversible flux distribution
if ~isempty(model_i)
    fd_ref_i    = zeros(length(model_i.rxns),1);
    for i=1:length(fd_ref)
        if length(model_i.rev2irrev{i}) > 1
            % reversible reaction
            if fd_ref(i) > 0
                % positive forward flux of reversible reaction
                fd_ref_i(ismember(model_i.rxns,[model.rxns{i},'_f']))   = fd_ref(i);
            elseif fd_ref(i) < 0
                % positive backward flux of reversible reaction
                fd_ref_i(ismember(model_i.rxns,[model.rxns{i},'_b']))   = -fd_ref(i);
            end
        else
            % irreversible reaction
            fd_ref_i(model_i.rev2irrev{i})  = fd_ref(i);
        end
    end
    model_i.fd_ref_i    = fd_ref_i;
    model_i.fd_ref      = fd_ref;
else
    % no irreversible model provided
    model_i     = [];
end


end