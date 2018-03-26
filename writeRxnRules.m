function writeRxnRules(opt_ga)
% function writes m-File to evaluate target reactions and get target bounds
% INPUT

% OUTPUT
% writes m- or mex-file

%%


fileID  = fopen([opt_ga.AddFilesDir,opt_ga.slash,'evalTargets.m'],'w');
fprintf(fileID,'function [targetRxnNum,targetRxnNum_i,targetRxnBounds,targetRxnBounds_i] = evalTargets(targetNum,targets)\n');
fprintf(fileID,['targetRxnNum_i = [targets.rxnNum_i{targetNum}];\n',...
                'targetRxnNum = targets.rxnNum(targetNum);\n',...
                'targetRxnBounds_i_cell = targets.bound_i(targetNum,:);\n',...
                'targetRxnBounds_i = [[targetRxnBounds_i_cell{:,1}]'',[targetRxnBounds_i_cell{:,2}]''];\n',...
                'targetRxnBounds = targets.bound(targetNum,:);\n',...
                ]);
fprintf(fileID,'end');      
fclose(fileID);

%% mex file compatible
% fileID  = fopen('evalTargets.m','w');
% fprintf(fileID,'function [targetRxnNum,targetRxnNum_i,targetRxnBounds,targetRxnBounds_i] = evalTargets(targetNum,bound,bound_i,rxnNum,rxnNum_i)\n');
% fprintf(fileID,['lenTargetNum = length(targetNum);\n',...
%                 'targetRxnNum_i = zeros(1,0);\n',...
%                 'for i=1:lenTargetNum\n\t',...
%                 'targetRxnNum_i = [targetRxnNum_i,rxnNum_i{targetNum(i)}];\n end\n',...
%                 'targetRxnNum_i = [rxnNum_i{targetNum}];\n',...
%                 'targetRxnNum = rxnNum(targetNum);\n',...
%                 'targetRxnBounds_i_cell  = cell(lenTargetNum,2);\n',...
%                 'for i=1:lenTargetNum\n\t',...
%                 'bound_save_lb = zeros(1,0);\n\t',...
%                 'bound_save_ub = zeros(1,0);\n\t',...
%                 'for j=1:length(bound_i{targetNum(i),1})\n\t\t',...
%                 'end\n',...
%                 'targetRxnBounds_i_cell{i,1} = bound_save_lb;\n\t',...
%                 'targetRxnBounds_i_cell{i,2} = bound_save_ub;\nend\n',...
%                 'targetRxnBounds_i = [[targetRxnBounds_i_cell{:,1}]'',[targetRxnBounds_i_cell{:,2}]''];\n',...
%                 'targetRxnBounds = bound(targetNum,:);\n',...
%                 ]);
% fprintf(fileID,'end');      
% fclose(fileID);
% 
% %                 'bound_save_lb = [bound_save_lb,bound_i{targetNum(i),1}(j)];\n\t\t',...
% %                 'bound_save_ub = [bound_save_ub,bound_i{targetNum(i),2}(j)];\n\t end\n\t',...
% 
% % translate into mex file in 'initializeFitFun.m'
% % i1  = zeros(1,K);
% % i2  = targets.bound;
% % i3  = num2cell(targets.bound);
% % i4  = targets.rxnNum;
% % i5  = num2cell(targets.rxnNum);
% % codegen -o evalTargets evalTargets.m -args {i1, i2, i3, i4, i5}
% % 
% % % delete m-file
% % % delete('evalTargets.m');

end