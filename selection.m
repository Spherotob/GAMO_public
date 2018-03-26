function [matingMat,popSort_slct] = selection(pop,pop_Tbin,pop_bin,popFit,Nslct,Ndel,numPairs,avoidDubl)
% conduct selection of the fittest parent chromosomes for creating mating
% pairs

%% INPUT
% pop:              Population containing intervention strategies as
%                   chromosomes
% pop_Tbin:         Population in binary intervention format
% pop_bin:          Population in binary format 
% popFit:           Fitness of each chromosome in the population
% Nslct:            Number of parent chromosomes to be kept during selection
% Ndel:             Number of parent chromosomes to be dismissed during selection
% numPairs:         Number of mating pairs (2 offsprings per pair)
% avoidDubl:        Flag if identical mating pairs should be avoided (1) or not (0)

%% OUTPUT
% matingMat:        Matrix containing pairs of selected parents for mating
%                   (rows).
% popSort_slct:     Indices of selected chromosomes for the next generation

%% General parameters
minProb     = 1/(10*Nslct);     % minimally guaranteed probability
maxSample   = 100;               % maximal number of resampling of second parent
%% Natural selection (remove most unfit parents)
numInt          = sum(pop_Tbin,2);  % number of unique interventions
[pop_I,pop_P]   = sortrows([popFit,-numInt]);     % sort population according to the fitness and number of interventions
pop_I           = flipud(pop_I(:,1));
pop_P           = flipud(pop_P);
% [pop_I,pop_P]   = sort(popFit,'descend');
% check diversity of population according to the fitness
popFit_slct     = pop_I(1:Nslct);
popSort_slct    = pop_P(1:Nslct);
% pop_slct        = pop(pop_P(1:Nslct),:);

    % if all zero           % if all entries are the same
if all(popFit_slct<=0) || all(popFit_slct==popFit_slct(1))
    % No non zero fitness, randomly choose parents for mating
    popFit_slct_norm    = ones(Nslct,1)./Nslct;
else
    % Roulette selection with a minimal non zero selection probability for
    % parents with zero fitness
    % Normalize fitness with best removed parent
    popFit_slct_norm    = popFit_slct-(pop_I(Nslct+1)*(pop_I(Nslct+1)>=0));
    if any(popFit_slct_norm<=0)
        % add fitness to guarantee minimal probability for zero-fitness parents
        popFit_slct_norm    = popFit_slct_norm+((sum(popFit_slct_norm)*minProb)/(1-(minProb*Nslct)));
    end
end

% calculate selection probability 
slctProb        = popFit_slct_norm./(sum(popFit_slct_norm));
slctProbSum1    = cumsum(slctProb,'reverse');
slctProbSum2    = zeros(Nslct,1);
slctProbSum2(1:(Nslct-1))   = slctProbSum1(2:end);

% calculate selection probability for second parent
slctProb_sec    = (ones(Nslct).*popFit_slct_norm')./((eye(Nslct)==0)*popFit_slct_norm);
slctProb_sec    = slctProb_sec.*(eye(Nslct)==0);    % exclude value at the position of the first parent
slctProb_sec    = slctProb_sec';
slctProb_sec_sum1   = cumsum(slctProb_sec,'reverse');
slctProb_sec_sum2   = zeros(Nslct);
slctProb_sec_sum2(1:(Nslct-1),:)    = slctProb_sec_sum1(2:end,:); 

% draw random numbers
RN  = rand(2,numPairs);

matingMat   = zeros(numPairs,2);    % matrix containing the mating partners (row)
% select first parents
[parents,~]     = find((RN(1,:)<=slctProbSum1) & (RN(1,:)>slctProbSum2));

matingMat(:,1)  = parents;

% select second parent while excluding first parent 
[parents2,~]    = find((RN(2,:)<=slctProb_sec_sum1(:,parents)) ...
                    & (RN(2,:)>slctProb_sec_sum2(:,parents)));
matingMat(:,2)  = parents2;


% Eventually try to avoid multiple same pairs
if avoidDubl==1
% select second parent while excluding first parent while avoiding pair
% dublicates
    matingMat_flip  = fliplr(matingMat);
    for i=1:numPairs
        dublicate   = all(matingMat==matingMat(i,:),2) | all(matingMat==matingMat_flip(i,:),2);
        if sum(dublicate)>1
            % duplicates detected

            RN      = rand(maxSample,1);    % draw random numbers
            sample  = 1;
            while sample <= maxSample
                parent1         = matingMat(i,1);   % first parent
                [parents2,~]    = find((RN(sample)<=slctProb_sec_sum1(:,parent1)) ...
                        & (RN(sample)>slctProb_sec_sum2(:,parent1)));
                 matingMat(i,2)     = parents2;

                 matingMat_act_flip     = fliplr(matingMat(i,:));
                 dublicate   = all(matingMat==matingMat(i,:),2) ...
                                | all(matingMat==matingMat_act_flip,2);

                 if sum(dublicate)>1
                     % still dublicates, resample
                     sample     = sample+1;
                 else
                     % jump to next pair
                     break;
                 end
            end
        end  
    end
end
% Transform to intervention identifiers 
if size(matingMat,1)==1
    matingMat   = popSort_slct(matingMat)';
else
    matingMat   = popSort_slct(matingMat);
end
   

end