% test Wiener-Khinchine's theorem on fractional Brownian motion data with
% missing samples (compute spectra)
% Cédric Chavanne, 2026-01-27.
% cedric.chavanne@ensta.org
close all
clear all
addpath ../tools

% parameters
p = [33]; % percent of missing data (must be the same as in do_synthetic_holes)
hurst = {'1over3','1over2'}; % Hurst parameter (name)

for j = 1:length(hurst)
    load(['../../data/fbm_',char(hurst(j))])
    [M,N] = size(y);
    % window for Fourier analysis
    w = tukeywin(M,0); % no window
    w = repmat(w,1,N);

    for i=1:length(p)
        load(['../../data/fbm_holes_',num2str(p(i)),'percent'])

        % replace missing data by nans
        y1 = y;
        y1(find(w1 == 0)) = nan;
        y2 = y;
        y2(find(w2 == 0)) = nan;

        % compute psd of sampling functions with periodogram
        [Pw1,f] = psd(w1);
        Pw2 = psd(w2);

        % compute psd of original data with periodogram
        Py = psd(w.*y);

        % compute psd via Wiener-Khinchin with unbiased circular autocorrelation:
        Py1c = psd_WKcirc(w.*y1);
        Py2c = psd_WKcirc(w.*y2);

        % compute psd via Wiener-Khinchin with unbiased standard autocorrelation:
        [Py1u,fs] = psd_WKunbi(w.*y1);
        Py2u = psd_WKunbi(w.*y2);

        % compute psd via Wiener-Khinchin with biased standard autocorrelation:
        Py1b = psd_WKbias(w.*y1);
        Py2b = psd_WKbias(w.*y2);

        % compute psd with Lomb-Scargle:
        ofac = 1;
        fmax = max(f);
        [Py1L,fL] = plomb(w.*y1,1,fmax,ofac);
        Py2L = plomb(w.*y2,1,fmax,ofac);

        % save results
        save(['../../data/fbm_',char(hurst(j)),'_missing_',num2str(p(i)),'percent_spectra'],'Pw1','Pw2','Py','Py1c','Py1b','Py1u','Py1L','Py2c','Py2b','Py2u','Py2L','f','fs','fL')
        
    end
end
