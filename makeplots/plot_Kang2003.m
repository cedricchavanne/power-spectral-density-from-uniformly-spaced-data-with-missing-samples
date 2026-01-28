% plot test case: decaying turbulence data from Kang et al. (2003)
% before executing this script, you must execute:
% test_Kang2003
% Cédric Chavanne, 2026-01-27.
% cedric.chavanne@ensta.org

clear all
close all
addpath ../tools
load('../../data/Kang2003_spectra')

% bin-average spectra
Nb = 25; % number of bins per decade for bin averaging
[Pub,fb] = binaverage(Pu,f,Nb); % bin average
Pu1cb = binaverage(Pu1c,f,Nb); % bin average
[Pu1bb,fsb] = binaverage(Pu1b,fs,Nb); % bin average
Pu1ub = binaverage(Pu1u,fs,Nb); % bin average
Pu1Gb = binaverage(Pu1G,fs,Nb); % bin average
[Pu1Lb,fLb] = binaverage(Pu1L,fL,Nb); % bin average
Pu2cb = binaverage(Pu2c,f,Nb); % bin average
Pu2bb = binaverage(Pu2b,fs,Nb); % bin average
Pu2ub = binaverage(Pu2u,fs,Nb); % bin average
Pu2Gb = binaverage(Pu2G,fs,Nb); % bin average
Pu2Lb = binaverage(Pu2L,fL,Nb); % bin average

% figures
Fs = 10; % font size
Jb = find(fb>5);

figure(1)
clf

subplot(121)
h2 = loglog(f,Pu,'--k','linewidth',2);
axis tight
Ax = axis;
Ax(4) = 2*Ax(4);
axis(Ax)
hold on
h3 = loglog(fsb,Pu1bb,'b','linewidth',2);
h4 = loglog(fsb,Pu1ub,'c','linewidth',2);
h5 = loglog(fb,Pu1cb,'r','linewidth',2);
%h6 = loglog(fLb,2e4*Pu1Lb,'m');
h7 = loglog(fsb,Pu1Gb,'--g','linewidth',2);
h2 = loglog(f,Pu,'--k','linewidth',2);
h1 = loglog(fb(Jb),2e5*fb(Jb).^(-5/3),'--k');
set(gca,'xtick',[1,1e1,1e2,1e3,1e4],'fontsize',Fs)
text(5e-2*Ax(2),5e-1*Ax(4),'(a)','fontsize',Fs)
text(5e-2*Ax(2),2e-1*Ax(4),['{\itp} = 50%'],'fontsize',Fs)
text(5e-2*Ax(2),1e-1*Ax(4),['{\itM} = ',num2str(n)],'fontsize',Fs)
text(5e-2*Ax(2),4e-2*Ax(4),['{\itN_b} = ',num2str(Nb)],'fontsize',Fs)
hl = legend([h1 h2 h3 h4 h5 h7 ],'{\itf}^{-5/3}','complete','{\itS_{uu}^b}','{\itS_{uu}^u}','{\itS_{uu}^c}','{\itS_{uu}^G}','location','southwest');
%hl = legend([h1 h2 h3 h4 h5 h6 h7 ],'{\itf}^{-5/3}','complete','{\itS_{xx}^b}','{\itS_{xx}^u}','{\itS_{xx}^c}','{\itS_{xx}^L}','{\itS_{xx}^G}','location','southwest');
legend boxoff
set(hl,'fontsize',Fs)
title(['Bernoulli missing ({\itw_a})'],'fontsize',Fs)
ylabel('<{\itS_{uu}}> [m^2 s^{-1}]','fontsize',Fs)
xlabel('{\itf} [s^{-1}]','fontsize',Fs)

subplot(122)
h2 = loglog(f,Pu,'--k','linewidth',2);
hold on
h3 = loglog(fsb,Pu2bb,'b','linewidth',2);
h4 = loglog(fsb,Pu2ub,'c','linewidth',2);
h5 = loglog(fb,Pu2cb,'r','linewidth',2);
%h6 = loglog(fLb,2e4*Pu2Lb,'m');
h7 = loglog(fsb,Pu2Gb,'--g','linewidth',2);
h2 = loglog(f,Pu,'--k','linewidth',2);
h1 = loglog(fb(Jb),2e5*fb(Jb).^(-5/3),'--k');
axis(Ax)
set(gca,'xtick',[1,1e1,1e2,1e3,1e4],'fontsize',Fs)
text(5e-2*Ax(2),5e-1*Ax(4),'(b)','fontsize',Fs)
title(['batch-Bernoulli missing ({\itw_b})'],'fontsize',Fs)
xlabel('{\itf} [s^{-1}]','fontsize',Fs)

print -f1 -depsc2 ../../figs/Kang2003_spectra


