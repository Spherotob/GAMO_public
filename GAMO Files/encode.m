function pop_bin = encode(pop,numBits,numBits_hri,encodeVec_midPos,K,K_hri)
% encodes a population to a binary format with a certain number of bits per
% gene

% INPUT
% pop:                  population containing intervention indices of each chromosome (rows)
% numBits:              Number of bits per deletion gene (intervention index)
% numBits_hri:          number of bits per insertion target
% encodeVec_midPos:     middle position of a target in encoding vector
% K:                    Number of deletion targets
% K_hri                 Number of insertion targets

% OUTPUT
% pop_bin:    population in binary format

%% General parameter
Np  = size(pop,1);  % population size

%% Rewrite population according to discretization scheme
pop_discr       = encodeVec_midPos(pop);

%% encode whole population
% number of indices for deletion targets in encoding vector
numIdx      = (2^numBits)*(numBits>0);
%
pop_bin     = zeros(Np,(K*numBits)+(K_hri*numBits_hri));
% deletion targets
pos             = 1;
for i=1:K
    pos2                 = pos+numBits-1;
    pop_bin(:,pos:pos2)  = de2bi(pop_discr(:,i),numBits); 
    pos                  = pos2+1;
end
% insertion targets
for i=(K+1):(K+K_hri)
    pos2                 = pos+numBits_hri-1;
    pop_bin(:,pos:pos2)  = de2bi(pop_discr(:,i)-numIdx,numBits_hri); 
    pos                  = pos2+1;
end

end