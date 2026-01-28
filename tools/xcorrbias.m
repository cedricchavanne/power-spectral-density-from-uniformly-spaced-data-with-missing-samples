function [c,l,n]=xcorrbias(x,varargin)

% [c,l,n]=xcorrbias(x,y)
% Computes the biased standard cross-correlation function c 
% of the data series x and y.
% Missing data should be NaN's.
% l gives the lags.
% n gives the number of data pairs available for computation for each lag.
% c is normalized by n(l=0).
% x and y can be 2D matrices (of the same size), 
% with rows corresponding to data samples
% and columns corresponding to different realizations.
% In that case, each column of c is the cross-correlation function
% of the corresponding column of x and y.
%
% [c,l,n]=xcorrbias(x)
% Computes the biased standard auto-correlation function c 
% of the data series x.

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

% replace missing data by zeros
I = find(isnan(x));
x(I) = 0;
J = find(isnan(y));
y(J)=0;

% sampling functions
zx=ones(size(x));
zx(I)=0;
zy=ones(size(x));
zy(J)=0;

% number of data pairs available to compute the correlation
% rounded to nearest integer
n = round(real(ifft(fft(zx,2^nextpow2(2*M-1)).*conj(fft(zy,2^nextpow2(2*M-1))))));

% compute the biased standard correlation, normalized by n(l=0)
c = ifft(fft(x,2^nextpow2(2*M-1)).*conj(fft(y,2^nextpow2(2*M-1))))./repmat(n(1,:),size(n,1),1);

% move negative lags before positive lags
c = [c(end-M+2:end,:);c(1:M,:)];
n = [n(end-M+2:end,:);n(1:M,:)];

% lags
l = [-M+1:M-1]';

% If first vector is a row, return a row
c = shiftdim(c,-nshift);
l = shiftdim(l,-nshift);
n = shiftdim(n,-nshift);

