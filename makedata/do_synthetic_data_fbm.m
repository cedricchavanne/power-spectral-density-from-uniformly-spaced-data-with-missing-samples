% generate synthetic data for fractional Brownian motion with 
% different Hurst parameters
% Cédric Chavanne, 2026-01-27.
% cedric.chavanne@ensta.org

%% H = 1/3 
clear all

% parameters
M = 2^10; % number of data points
N = 10000; % number of realizations
H = 1/3; % Hurst parameter

% generate data 
y = synthetic_data_fbm(M,N,H);

% save data to be able to redraw same figures
save('../../data/fbm_1over3','y')

%% H = 1/2 
clear all

% parameters
M = 2^10; % number of data points
N = 10000; % number of realizations
H = 1/2; % Hurst parameter

% generate data 
y = synthetic_data_fbm(M,N,H);

% save data to be able to redraw same figures
save('../../data/fbm_1over2','y','H')

