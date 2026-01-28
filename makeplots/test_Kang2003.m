% test case: decaying turbulence data from Kang et al. (2003)
% Cédric Chavanne, 2026-01-27.
% cedric.chavanne@ensta.org

clear all
close all
addpath ../tools
load('../../data/Kang2003.mat');
load('../../data/Kang2003_holes_50percent.mat');
m = 2^15; % segments length
n = floor(size(u,1)/m); % number of contiguous segments
n = n+n-1; % number of half-overlapping segments

% rearrange in half-overlapping segments
U = repmat(nan,m,n);
for k = 1:n
    U(:,k) = u(1+(k-1)*m/2:(k-1)*m/2+m);
end
u = U;
clear U

% detrend for Fourier analysis
u = detrend(u);

% replace missing data with nans
u1 = u;
u1(find(w1==0)) = nan;
u2 = u;
u2(find(w2==0)) = nan;

% window
%w = bartlett(m);
w = tukeywin(m,0); % no window
w = repmat(w,1,n);

% PSD Wiener-Khinchin with circular unbiased correlation
[Pu,f] = psd_WKcirc(w.*u);
Pu1 = psd_WKcirc(w.*u1);
Pu2 = psd_WKcirc(w.*u2);
f = f/dt;
% one-sided averaged spectra:
I = find(f>0);
Pu = 2*mean(Pu(I,:),2);
Pu1c = 2*mean(Pu1(I,:),2);
Pu2c = 2*mean(Pu2(I,:),2);
f = f(I);

% PSD Wiener-Khinchin with standard unbiased correlation
[Pu1,fs] = psd_WKunbi(w.*u1);
Pu2 = psd_WKunbi(w.*u2);
fs = fs/dt;
% one-sided averaged spectra:
I = find(fs>0);
Pu1u = 2*mean(Pu1(I,:),2);
Pu2u = 2*mean(Pu2(I,:),2);
fs = fs(I);
% Gao estimator
Pu1G = 2*mean(abs(Pu1(I,:)),2);
Pu2G = 2*mean(abs(Pu2(I,:)),2);

% PSD Wiener-Khinchin with standard biased correlation
[Pu1,fs] = psd_WKbias(w.*u1);
Pu2 = psd_WKbias(w.*u2);
fs = fs/dt;
% one-sided averaged spectra:
I = find(fs>0);
Pu1b = 2*mean(Pu1(I,:),2);
Pu2b = 2*mean(Pu2(I,:),2);
fs = fs(I);

% PSD Lomb-Scargle
ofac = 1;
fmax = max(f);
[Pu1L,fL] = plomb(w.*u1,1/dt,fmax,ofac);
[Pu2L,fL] = plomb(w.*u2,1/dt,fmax,ofac);
% one-sided averaged spectra:
I = find(fL>0);
Pu1L = 2*mean(Pu1L(I,:),2);
Pu2L = 2*mean(Pu2L(I,:),2);
fL = fL(I);

% save results
save ../../data/Kang2003_spectra Pu Pu1b Pu1c Pu1G Pu1L Pu1u Pu2b Pu2c Pu2G Pu2L Pu2u f fL fs n

