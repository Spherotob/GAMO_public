function offspring_bin = mating(pop_bin,K_tot,matingMat,numPairs,numKntChr,crossType,enc)
% Conduct mating of mating pairs to create offsprings for the next
% generation. Mating includes crossover of parents genome

%% INPUT
% pop_bin:          Selected parent Population in binary format 
% K:                Number of genes per chromosome
% matingMat:        Matrix containing pairs of selected parents for mating
%                   (rows).
% numPairs:         Number of mating pairs (2 offsprings per pair)
% numKntChr:        Number of kinetochores (max. 2)
% crossType:        Type of crossover routine
% enc:              Struct containing encoding information for the
%                   conversion from and to the binary format


%% OUTPUT
% offspring_bin:    Offspring chromosomes in binary format

%% General parameters
numBits_tot     = enc.numBits_tot;

% Choose crossover type
% 0:  Crossover of complete genes
% 1:  Crossover at random position within the binary representation of a
%     chromosome.

%% Define crossover positions for each mating pair
% if numKntChr==1 -> single crossover
% if numKntChr==2 -> double crossover
if crossType==0
    pos_cross           = ceil(rand(numPairs,numKntChr).*(K_tot-1));    % randomly select crossover point
    pos_cross_unique    = [1:(K_tot-1)]';
    pos_cuts            = enc.posGene;
elseif crossType==1
    pos_cross           = ceil(rand(numPairs,numKntChr).*(numBits_tot-1));  % randomly select crossover point
    pos_cross_unique    = unique(pos_cross);
    pos_cuts            = 1:numBits_tot;
end


%% conduct crossover
offspring1  = zeros(numPairs,numBits_tot);
offspring2  = zeros(numPairs,numBits_tot);
if numKntChr==1
    
    % single kinetochor
    for i=1:length(pos_cross_unique)
        pos_i               = pos_cross==pos_cross_unique(i);
        offspring1(pos_i,:) = [pop_bin(matingMat(pos_i,1),1:(pos_cuts(pos_cross_unique(i)+1)-1))...
                               ,pop_bin(matingMat(pos_i,2),pos_cuts(pos_cross_unique(i)+1):numBits_tot)];
        offspring2(pos_i,:) = [pop_bin(matingMat(pos_i,2),1:(pos_cuts(pos_cross_unique(i)+1)-1))...
                               ,pop_bin(matingMat(pos_i,1),pos_cuts(pos_cross_unique(i)+1):numBits_tot)];               
    end
elseif numKntChr>1
    % multiple kinetochors
    [pos_cross,~]   = sort(pos_cross,2);    % sort crossover points
    pos_cross_add   = [zeros(numPairs,1),pos_cuts(pos_cross+1)-1,ones(numPairs,1).*numBits_tot];    % add start and end point for easy indexing
    mating_choice   = mod(1:(numKntChr+1),2);   % switches thw two parents for gene annotation
    for i=1:numPairs
        for j=1:(numKntChr+1)
            pos_chr     = (pos_cross_add(i,j)+1):pos_cross_add(i,j+1);  % pick genes to crossover
            offspring1(i,pos_chr)     = pop_bin(matingMat(i,mating_choice(j)+1),...
                                            pos_chr);
            offspring2(i,pos_chr)     = pop_bin(matingMat(i,~mating_choice(j)+1),...
                                            pos_chr);                            
        end
    end 
end
% add offspring to selected parents
offspring_bin     = [offspring1;offspring2];
% encode offspring and add to selected parents
% a           = decode([offspring1;offspring2],enc.numBits,enc.encodeVec);
% pop         = [pop;decode([offspring1;offspring2],enc.numBits,enc.encodeVec)];     



end