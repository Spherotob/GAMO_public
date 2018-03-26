function [mutRate,hd] = adaptMutRate(pop_Tbin,hd_prev,staticGen_c,opt_mut)
% adapt mutation rate according to the hamming distance and the change in
% hamming distance comapred to previous generation

% calculate actual hamming distance
hd  = sum(pdist(pop_Tbin,'hamming'))*opt_mut.hd_fac;
% adapt mutation rate
mutRate     = -opt_mut.adaptMut_P...
                +(opt_mut.adaptMut_P*hd)...
                +(opt_mut.adaptMut_I*staticGen_c);

end