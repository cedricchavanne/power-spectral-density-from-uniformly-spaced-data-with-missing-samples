function [w1,w2] = synthetic_holes(N,M,p)

% [w1,w2] = synthetic_holes(N,M,p)
%
% generate random synthetic holes with Bernoulli and batch-Bernoulli 
% distributions with an exact number of missing data.
%
% Inputs:
% N: number of data samples
% M: number of realizations
% p: percent of missing data (0 to 100)
%
% Outputs:
% w1: MxN sampling matrix (1 for existing data, 0 for missing data)
%     with Bernoulli distributed holes
% w2: MxN sampling matrix (1 for existing data, 0 for missing data)
%     with batch-Bernoulli distributed holes

% Cédric Chavanne, 2026-01-26.
% cedric.chavanne@ensta.org

Nh = round(N*p/100); % number of missing points

% generate Bernoulli distributed holes
w1 = ones(M,N);
for i=1:M
    I = randsample(N-2,Nh); % avoid hole being last point (see also below)
    I = I+1; % avoid hole being first point
    w1(i,I) = 0;
end

% generate batch-Bernoulli distributed holes
w2 = ones(M,N);
for i=1:M
    % lengths of holes:
    L = [];
    while sum(L)<Nh
        l = randsample(Nh-sum(L),1);
        L = [L,l];
    end
    L = L(randperm(length(L))); % reorder lengths of holes randomly
    % lengths of good data:
    G = [];
    while length(G) ~= length(L)+1  % ensure there are good data before and after the gaps
        G = [];
        while sum(G)<N-Nh
            g = randsample(N-Nh-sum(G),1);
            G = [G,g];
        end
    end
    G = G(randperm(length(G))); % reorder lengths of good data randomly
    % position holes:
    j0 = 1;
    for j = 1:length(L)
        w2(i,j0+G(j):j0+G(j)+L(j)-1) = 0;
        j0 = j0+G(j)+L(j);
    end   
end
