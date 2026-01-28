% generate synthetic holes for different percentages of missing data
% Cédric Chavanne, 2026-01-27.
% cedric.chavanne@ensta.org

%% for fractional Brownian motion data
clear all
% parameters
p = [33]; % percent of missing data
load ../../data/fbm_1over3
[M,N] = size(y)
clear y

% generate and save holes
for i=1:length(p)
    [w1,w2] = synthetic_holes(M,N,p(i));
    % transpose matrices and convert to logical
    w1 = logical(w1');
    w2 = logical(w2');
    save(['../../data/fbm_holes_',num2str(p(i)),'percent'],'w1','w2')
end

%% for Kang et al. (2003) data
clear all
% parameters
M = 2^15; % number of data points
N = 2195; % number of realizations
p = [50]; % percent of missing data

% generate and save holes
for i=1:length(p)
    [w1,w2] = synthetic_holes(M,N,p(i));
    % transpose matrices and convert to logical
    w1 = logical(w1');
    w2 = logical(w2');
    save(['../../data/Kang2003_holes_',num2str(p(i)),'percent'],'w1','w2')
end

