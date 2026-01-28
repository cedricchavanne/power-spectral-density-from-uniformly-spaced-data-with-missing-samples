function [c,l,n]=xcorrunbi(x,varargin)

% [c,l,n]=xcorrunbi(x,y)
% Computes the unbiased standard cross-correlation function c 
% of the data series x and y.
% Missing data should be NaN's.
% l gives the lags.
% n gives the number of data pairs available for computation for each lag.
% If n = 0 for a lag, c is set to 0 for this lag.
% x and y can be 2D matrices (of the same size), 
% with rows corresponding to data samples
% and columns corresponding to different realizations.
% In that case, each column of c is the cross-correlation function
% of the corresponding column of x and y.
%
% [c,l,n]=xcorrunbi(x)
% Computes the unbiased standard auto-correlation function c 
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

% compute the biased standard correlation
c = ifft(fft(x,2^nextpow2(2*M-1)).*conj(fft(y,2^nextpow2(2*M-1))));

% number of data pairs available to compute the correlation
% rounded to nearest integer
n = round(real(ifft(fft(zx,2^nextpow2(2*M-1)).*conj(fft(zy,2^nextpow2(2*M-1))))));

% compute the unbiased standard correlation
c = c./n;
% if no data pair is available for a lag, set c to zero
c(find(n==0)) = 0;
% % if no data pair is available for a lag, make c not-a-number
% c(find(n==0)) = nan;

% move negative lags before positive lags
c = [c(end-M+2:end,:);c(1:M,:)];
n = [n(end-M+2:end,:);n(1:M,:)];

% lags
l = [-M+1:M-1]';

% If first vector is a row, return a row
c = shiftdim(c,-nshift);
l = shiftdim(l,-nshift);
n = shiftdim(n,-nshift);

