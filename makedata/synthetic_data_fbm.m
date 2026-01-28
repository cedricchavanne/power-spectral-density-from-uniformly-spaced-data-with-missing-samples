function y = synthetic_data_fbm(M,N,H)

% y = synthetic_data_fbm(M,N,H)
%
% generate random synthetic data for fractional Brownian motion
% with zero mean value and linearly-detrended.
%
% Inputs:
% M: number of data samples
% N: number of realizations
% H: Hurst parameter (0 <= H <= 1) 
%
% Outputs:
% y: MxN data matrix

% Cédric Chavanne, 2026-01-26.
% cedric.chavanne@ensta.org

y = repmat(nan,M,N);
for j = 1:N
    y(:,j) = wfbm(H,M);
end

% detrend for Fourier analysis
y = detrend(y);
