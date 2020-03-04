%% Function to reduce the number of (KO) target reactions

function [noTarget] = redTargetSpace(model, noTarget, modelType)
% INPUTS
% noTarget:           Provided list of non-target reactions
% modelType:          Type of metabolic model concerning identifier standards
% coFacNum:           metabolite numbers of cofactors, energy equivalents etc

%% define subSystems irrelevant as knockout targets
nonTargetSubsystems = {'Cell Envelope Biosynthesis',...
                        'Transport',...
                        'transport',...
                        'Membrane Lipid Metabolism',...
                        'Murein Biosynthesis',...
                        'tRNA Charging',...
                        'Glycerophospholipid Metabolism'};
                        
%%
switch modelType
    case 0  % Bigg identifiers (like iAF1260)
        
        % Exclude certain types of reactions, namly
        %   spontaneous, transport, diffusion, exchange
        noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.rxnNames,'exchange')))))];
        noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.rxnNames,'diffusion')))))];
        noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.rxnNames,'transport')))))];  
        noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.rxnNames,'spontaneous')))))]; 
        

        % Exclude reactions belonging to a certain subsystem
        if isfield(model,'subSystems')
            for i=1:length(model.subSystems)
                % identify data type
                if ischar(model.subSystems{i})
                    for t=1:length(nonTargetSubsystems)
                        if ~isempty(strfind(model.subSystems{i},nonTargetSubsystems{t}))
                            % exclude reaction
                            noTarget    = [noTarget;i];
                            break;
                        end
                    end


                elseif iscell(model.subSystems{i})           
                    for t=1:length(nonTargetSubsystems)
                        if any(not(cellfun('isempty',(strfind(model.subSystems{i},nonTargetSubsystems{t})))))
                            % exclude reaction
                            noTarget    = [noTarget;i];
                            break;
                        end
                    end
                end
            end
        end 

        % only consider gene related targets
        % exclude reactions that are not related to a gene
        if isfield(model,'grRules')
            geneRel     = cellfun('isempty',model.grRules);
            if any(not(geneRel))
                noTarget    = [noTarget; find(geneRel)];
            end
        end
 
     
        
end
        
noTarget   = unique(noTarget);  
        
end