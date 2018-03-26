function pop = decode(pop_bin,enc)
% decodes a population from a binary format to actual target numbers

% INPUT
% pop_bin:    population in binary format
% enc:        encoding parameter struct

% OUTPUT
% pop:        population 

%% General parameter
Np          = size(pop_bin,1);  % population size and number bits per chromosome
K_tot       = enc.K_tot;      % number of interventions
numBits     = enc.numBits;
K           = enc.K;

%% decode whole population
% number of indices for deletion targets in encoding vector
numIdx      = (2^numBits)*(numBits>0);

pop         = zeros(Np,K_tot);
pos         = 1;
for i=1:K
    pos2        = pos+numBits-1;
    pop(:,i)    = bi2de(pop_bin(:,pos:pos2));
    pos         = pos2+1;
end
for i=(K+1):K_tot
    pos2        = pos+enc.numBits_hri-1;
    pop(:,i)    = bi2de(pop_bin(:,pos:pos2))+numIdx;
    pos         = pos2+1;
end

%% match discretization with intervention indices
pop     = enc.encodeVec(pop+1);

end

