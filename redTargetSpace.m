%% Function to reduce the number of (KO) target reactions

function [noTarget] = redTargetSpace(model, noTarget, modelType)
% INPUTS
% noTarget:           Provided list of non-target reactions
% modelType:          Type of metabolic model concerning identifier standards
% coFacNum:           metabolite numbers of cofactors, energy equivalents etc


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
            noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.subSystems,'Cell Envelope Biosynthesis')))))];
            noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.subSystems,'Transport')))))];
            noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.subSystems,'transport')))))];
            noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.subSystems,'Membrane Lipid Metabolism')))))];
            noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.subSystems,'Murein Biosynthesis')))))];
            noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.subSystems,'tRNA Charging')))))];
            noTarget = [noTarget; find(not(cellfun('isempty',(strfind(model.subSystems,'Glycerophospholipid Metabolism')))))];
        end
        
        % exclude reactions that are not related to a gene
        if isfield(model,'grRules')
            geneRel     = cellfun('isempty',model.grRules);
            if any(not(geneRel))
                noTarget    = [noTarget; find(geneRel)];
            end
        end
        
%         % exclude reactions that act on metabolite with more then 7
%         % carbons (if sepc
% %         tar     = [];
%         if isfield(model,'metFormulas')
%             isMetFormulas   = not(cellfun('isempty',model.metNames));
%             if any(isMetFormulas)
%                 for i=1:length(isMetFormulas)
%                     if isMetFormulas(i)
%                         met     = model.metNames{i};
%                         posC     = strfind(met,'C');
%                         if ~isempty(posC)
%                             for k=1:length(posC)
%                                 % determine number after a 'C'
%                                 for j=(posC(k)+1):length(met)
%                                     if isempty(sscanf(met(j),'%g'))
%                                         % break if character is encountered
%                                         posNum  = j-1;
%                                         break;
%                                     end
%                                 end
%                                 if posNum ~= posC(k)
%                                     num     = str2double(met((posC(k)+1):posNum));
%                                     if num >= 7
%                                         % determine related reaction of the metabolite
% %                                         tar     = [tar; i, num];
%                                         noTarget    = [noTarget; find(model.S(i,:))'];
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
                
            
        
        
end
        
noTarget   = unique(noTarget);  
        
end