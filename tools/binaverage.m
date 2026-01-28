function [Pa,fa] = binaverage(P,f,N)

% [Pa,fa] = binaverage(P,f,N)
% Bin-average power spectrum P(f) using N frequency bins per decade.
% The averaged power spectrum is Pa(fa), 
% where fa are the mean frequencies of the bins.

% Cédric Chavanne, 2026-01-26.
% cedric.chavanne@ensta.org
% Adapted from Gao et al. (2021), Scaling Analysis of the China France Oceanography 
% Satellite Along-Track Wind and Wave Data, Journal of Geophysical Research: Oceans.


if nargin ~= 3
    error('Three inputs are required');
end
df = 1/N;
f1 = log10(f);
bin = min(f1):df:max(f1)+df;
NB = length(bin);
fa = zeros(1,NB-1)*nan;
Pa = zeros(1,NB-1)*nan;
for i = 1:NB-1
    I = find(f1 >= bin(i) & f1 <= bin(i+1));
    if ~isempty(I)
        Pa(i) = nanmean(P(I));
        fa(i) = 10^(nanmean(f1(I)));
    end
end
% remove empty values
I = find(isnan(Pa));
Pa(I) = [];
fa(I) = [];
end