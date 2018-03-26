function [pop,pop_Tbin,pop_bin] = mutation(pop_bin,popFit_slct,Np,Nt,K,Nslct,enc,mutRate,elite,targets,chr_map,opt)
% mutate population by changing single bits of the binary format (0->1 or
% 1->0)

%% INPUT
% pop_bin:          Population in binary format 
% Np:               Number of chromosomes in the population
% Nt:               Number of possible targets
% Nslct:            Number of chromosomes selected for mating
% enc:              Struct containing encoding information for the
%                   conversion from and to the binary format

% mutRate:          Mutation rate (0-1) related to the whole population and
%                   its number of bits
% elite:            Number of elite chromosomes which are not to be mutated
% targets:          Struct with possible targets
% exclOld:          Flag if old solutions should be avoided

%% OUTPUT
% pop_bin:          Mutated population in binary format 

%% general parameter
maxSample   = 100;  % maximal sample iterations
Nt_tot      = targets.Nt_tot;
numBits_tot = enc.numBits_tot;

%% Conduct mutations
mut_bin             = [ones(elite,numBits_tot);rand(Np-elite,numBits_tot)]<mutRate;
pop_bin(mut_bin)    = abs(pop_bin(mut_bin)-1);


%% Decode new population
pop                 = decode(pop_bin,enc);
pop_Tbin            = zeros(Np,Nt_tot); 

for i=1:Np
    pop_Tbin(i,pop(i,:))    = 1;
end

%% mutate positions of reactions that occur multiple times and mutate new solutions that have occured before
% check for duplicates (reaction target) within each chromosome
% try
    rxnFlag     = pop_Tbin(:,1:Nt)*targets.map;
    rxnFlag     = any(rxnFlag>1,2);
% catch
%     save RESCUE
%     return
% end

rxnFlag_pos = find(rxnFlag);

if ~isempty(rxnFlag_pos)
    num_pos             = length(rxnFlag_pos);
    mut_num             = ceil(mutRate*numBits_tot); % recalculate mutation number for single chromosomes
    for i=1:num_pos
        pop_bin_s   = pop_bin(rxnFlag_pos(i),:);
        sample      = 1;
        while sample<maxSample
            RN          = rand(mut_num,1);  % draw random numbers
            mut_bit     = ceil(RN*numBits_tot);  % randomly choose bits to be mutated
            pop_bin_s(mut_bit)    = abs(pop_bin_s(mut_bit)-1);
            pop_s                 = decode(pop_bin_s,enc);
            % check if duplicates or old solutions exist after mutation
            test_flag           = zeros(1,Nt_tot);
            test_flag(pop_s)    = 1;
            
            rxnFlag     = any((test_flag(:,1:Nt)*targets.map) > 1);
            
                            
            if rxnFlag
                sample  = sample+1;
                if sample>=maxSample
                    warning('Mutation routine: Maximal sample size exceeded!')
                end
            else 
                sample  = maxSample+1;
            end
        end
        pop(rxnFlag_pos(i),:)       = pop_s;
        pop_Tbin(rxnFlag_pos(i),:)  = test_flag;
        pop_bin(rxnFlag_pos(i),:)   = pop_bin_s;
    end
end


end