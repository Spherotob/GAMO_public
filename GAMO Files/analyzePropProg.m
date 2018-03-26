function [plateau_gen,x,resnorm] = analyzePropProg(popFit_mean)
% Analysis of population progression of past generations
% function evaluates if GA reached a plateau and stalls
% fit mean fitness progression to hyperbole function


%% conduct linear regressions of progressing data
% parameter
numGen  = length(popFit_mean);  % number of past generations
limit   = 0.04;     % maximum gradient indicating plateau
plateau_gen    = -1;

% normalize data by maximal mean
popFit_mean_norm    = popFit_mean./max(popFit_mean);

x       = zeros(numGen-1,2);
resnorm = ones(numGen-1,1);
for i=1:(numGen-3)
    if (numGen-i)<50
        v   = i:numGen;
    else
        v   = i:(i+50);
    end
    
    % linear fit
    [x(i,:),S,mu]       = polyfit(v',popFit_mean_norm(v,1),1);
    y                   = polyval([x(i,:)],v,S,mu);
    resnorm(i)          = (sum((y'-popFit_mean_norm(v,1)).^2))/length(v);
    
    % detect plateau region
    if i>50
        % negative gradient criteria
        if (x(i,1)<=0 && (numGen-i)>=50)
            plateau_gen     = i;
            disp('Mean fitness among population reached a plateau')
            break
        end
        % low gradient criteria
        if all(x((i-50):i,1)<limit)
            plateau_gen    = i-50;
            disp('Mean fitness among population reached a plateau')
            break
        end
    end
end




end