% test Wiener-Khinchine's theorem on fractional Brownian motion data
% Cédric Chavanne, 2026-01-27.
% cedric.chavanne@ensta.org

close all
clear all
addpath ../tools
load ../../data/fbm_1over3
ya = y;
load ../../data/fbm_1over2
yb = y;
clear y

[N,M] = size(ya);
x = 0:N-1;
w = tukeywin(N,0); % no window
w = repmat(w,1,M);

figure(1)
clf
subplot(211)
plot(x,ya(:,1))
subplot(212)
plot(x,yb(:,1))

figure(2)
clf
plot(x,w(:,1))

% compute power spectral density with periodogram
[pa0,k0] = psd(w.*ya);
[pb0,k0] = psd(w.*yb);
J = find(k0>0);
k0 = k0(J);
pa0 = 2*mean(pa0(J,:),2); % average one-sided spectrum
pb0 = 2*mean(pb0(J,:),2); % average one-sided spectrum

% compute psd via Wiener-Khinchin with circular unbiased correlation
[pa1,k1] = psd_WKcirc(w.*ya);
[pb1,k1] = psd_WKcirc(w.*yb);
J = find(k1>0);
k1 = k1(J);
pa1 = 2*mean(pa1(J,:),2); % average one-sided spectrum
pb1 = 2*mean(pb1(J,:),2); % average one-sided spectrum

% compute psd via Wiener-Khinchin with standard biased correlation
[pa2,k2] = psd_WKbias(w.*ya);
[pb2,k2] = psd_WKbias(w.*yb);
J = find(k2>0);
k2 = k2(J);
pa2 = 2*mean(pa2(J,:),2); % average one-sided spectrum
pb2 = 2*mean(pb2(J,:),2); % average one-sided spectrum

% compute psd via Wiener-Khinchin with standard unbiased correlation
[pa3,k3] = psd_WKunbi(w.*ya);
[pb3,k3] = psd_WKunbi(w.*yb);
J = find(k3>0);
k3 = k3(J);
pa3 = 2*mean(pa3(J,:),2); % average one-sided spectrum
pb3 = 2*mean(pb3(J,:),2); % average one-sided spectrum

% plot spectra
fs = 12; % font size

figure(2)
clf
subplot(122)
h2 = loglog(k0,pb0,'--k','linewidth',2);
axis tight
Ax = axis;
Ax(3) = Ax(3)/2;
Ax(4) = 2*Ax(4);
axis(Ax)
hold on
h4 = loglog(k3,pb3,'c','linewidth',2);
h3 = loglog(k2,pb2,'b','linewidth',2);
h5 = loglog(k1,pb1,'r','linewidth',2);
h2 = loglog(k0,pb0,'--k','linewidth',2);
h1 = loglog(k0,pb0(4)*(k0/k0(4)).^(-2),'--k');
legend([h1],'{\itf}^{-2}','location','southwest')
legend boxoff
%text(1.5*Ax(1),10*Ax(3),['{\itM} = ',num2str(M)],'fontsize',fs)
text(2*Ax(1),0.7*Ax(4),'b)','fontsize',fs)
xlabel('{\itf} [tu^{-1}]','fontsize',fs)
set(gca,'xtick',[1e-3,1e-2,1e-1,5e-1],'fontsize',fs)
title('Fractional Brownian motion with H = 1/2') 

subplot(121)
h2 = loglog(k0,pa0,'--k','linewidth',2);
hold on
h3 = loglog(k2,pa2,'b','linewidth',2);
h4 = loglog(k3,pa3,'c','linewidth',2);
h5 = loglog(k1,pa1,'r','linewidth',2);
h2 = loglog(k0,pa0,'--k','linewidth',2);
h1 = loglog(k2,pa2(4)*(k2/k2(4)).^(-5/3),'--k');
axis(Ax)
legend([h1 h2 h3 h4 h5],'{\itf}^{-5/3}','{\itS_{xx}^p}','{\itS_{xx}^b}','{\itS_{xx}^u}','{\itS_{xx}^c}','location','southwest')
legend boxoff
text(2*Ax(1),0.3*Ax(4),['{\itM} = ',num2str(M)],'fontsize',fs)
text(2*Ax(1),0.7*Ax(4),'a)','fontsize',fs)
ylabel('<{\itS_{xx}}> [au^2 tu]','fontsize',fs)
xlabel('{\itf} [tu^{-1}]','fontsize',fs)
set(gca,'xtick',[1e-3,1e-2,1e-1,5e-1],'fontsize',fs)
title('Fractional Brownian motion with H = 1/3') 

print -f2 -depsc2 ../../figs/test_WK_spectra

