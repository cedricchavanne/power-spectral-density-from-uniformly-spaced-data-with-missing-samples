function [P,f] = psd(x,varargin)

% [P,f] = psd(x,y)
% Computes the power cross-spectral density P of the data series x and y,
% and also returns the frequencies f corresponding to P,
% assuming a data spacing dt = 1. 
% If dt differs from 1, you should divide f by dt.
% Missing data (NaN) are NOT allowed.
% x and y can be 2D matrices (of the same size), 
% with rows corresponding to data samples
% and columns corresponding to different realizations.
%
% [P,f] = psd(x)
% Computes the power auto-spectral density P of the data series x.

% Cédric Chavanne, 2026-01-26.
% cedric.chavanne@ensta.org

error(nargchk(1,2,nargin));
if nargin == 1
    y = x;
else
    y = varargin{1};
end

[M,N] = size(x);
if size(y,1) ~= M | size(y,2) ~= N
    error('inputs must have same size')
end

% if horizontal vectors, make them vertical
if M==1
    [x,nshift] = shiftdim(x);
    y = shiftdim(y);
    M = N;
else
    nshift = 0;
end

% compute power cross-spectral density
X = fft(x);
Y = fft(y);
P = X.*conj(Y)/M;

% compute frequencies
if (mod(M,2)==0)
    f=[0 [1:M/2]/M [-M/2+1:-1]/M]';
else
    f=[0 [1:(M-1)/2]/M [(-M+1)/2:-1]/M]';
end

% If first vector is a row, return a row
P = shiftdim(P,-nshift);
f = shiftdim(f,-nshift);
