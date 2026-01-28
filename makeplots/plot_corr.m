%% plot auto-correlations for fractional Brownian motion with H = 1/3
% Cédric Chavanne, 2026-01-27.
% cedric.chavanne@ensta.org
close all
clear all
load ../../data/fbm_1over3

% window
[N,M] = size(y);
w = tukeywin(N,0); % no window
w = repmat(w,1,M);

figure(1)
clf
plot(w(:,1))

% circular unbiased correlation
[c1,l1] = xcorrcirc(w.*y);
c1 = mean(c1,2); % average cross-correlation

% standard biased correlation
[c2,l2] = xcorrbias(w.*y);
c2 = mean(c2,2); % average cross-correlation

% standard unbiased correlation
[c3,l3] = xcorrunbi(w.*y);
c3 = mean(c3,2); % average cross-correlation

% plot auto-correlations
fs = 12; % font size
I1 = find(l1 >= 0);
I2 = find(l2 >= 0);
I3 = find(l3 >= 0);

figure(3)
clf

subplot(211)
plot(l1(I1),c1(I1),'r','linewidth',4)
hold on
plot(l2(I2),c2(I2),'b','linewidth',2)
plot(l3(I3),c3(I3),'c','linewidth',2)
axis tight
Ax = axis;
plot(Ax(1:2),[0 0],'k--')
text(Ax(1)+0.05*(Ax(2)-Ax(1)),0.9*Ax(4),'(a)','fontsize',fs)
text(Ax(1)+0.1*(Ax(2)-Ax(1)),0.9*Ax(4),['{\itM} = ',num2str(M)],'fontsize',fs)
legend('<R_{xx,n}^c>','<R_{xx,n}^b>','<R_{xx,n}^u>','location','north')
legend boxoff
ylabel('<{\itR_{xx}}> [au^2]','fontsize',fs)
set(gca,'fontsize',fs)
title('Autocorrelation of fractional Brownian motion with H = 1/3') 

subplot(212)
plot(l1(I1),c1(I1),'r','linewidth',4)
hold on
h1 = plot(l2(I2),c2(I2),'b','linewidth',2);
h2 = plot(l2(I2(2:end)),flipud(c2(I2(2:end))),'b--','linewidth',2);
h3 = plot(l2(I2),[c2(I2(1)); c2(I2(2:end)) + flipud(c2(I2(2:end)))],'g','linewidth',2);
axis(Ax)
plot(Ax(1:2),[0 0],'k--')
text(Ax(1)+0.05*(Ax(2)-Ax(1)),0.9*Ax(4),'(b)','fontsize',fs)
legend([h2 h3],'<R_{xx,N-n}^b>','<R_{xx,n}^b> + <R_{xx,N-n}^b>','location','north')
legend boxoff
ylabel('<{\itR_{xx}}> [au^2]','fontsize',fs)
xlabel('Lag {\itn} [tu]','fontsize',fs)    
set(gca,'fontsize',fs)

print -f3 -depsc2 ../../figs/corr_fbm_1over3


%% plot auto-correlations for fractional Brownian motion with H = 1/2
close all
clear all
load ../../data/fbm_1over2

% window
[N,M] = size(y);
w = tukeywin(N,0); % no window
w = repmat(w,1,M);

figure(1)
clf
plot(w(:,1))

% circular unbiased correlation
[c1,l1] = xcorrcirc(w.*y);
c1 = mean(c1,2); % average cross-correlation

% standard biased correlation
[c2,l2] = xcorrbias(w.*y);
c2 = mean(c2,2); % average cross-correlation

% standard unbiased correlation
[c3,l3] = xcorrunbi(w.*y);
c3 = mean(c3,2); % average cross-correlation

% plot auto-correlations
fs = 12; % font size
I1 = find(l1 >= 0);
I2 = find(l2 >= 0);
I3 = find(l3 >= 0);

figure(3)
clf

subplot(211)
plot(l1(I1),c1(I1),'r','linewidth',4)
hold on
plot(l2(I2),c2(I2),'b','linewidth',2)
plot(l3(I3),c3(I3),'c','linewidth',2)
axis tight
Ax = axis;
plot(Ax(1:2),[0 0],'k--')
text(7,0.9*Ax(4),'(a)','fontsize',fs)
legend('<R_{xx,n}^c>','<R_{xx,n}^b>','<R_{xx,n}^u>','location','north')
legend boxoff
ylabel('<{\itR_{xx}}> [au^2]','fontsize',fs)
set(gca,'fontsize',fs)
title('Autocorrelation of fractional Brownian motion with H = 1/2') 

subplot(212)
plot(l1(I1),c1(I1),'r','linewidth',4)
hold on
h1 = plot(l2(I2),c2(I2),'b','linewidth',2);
h2 = plot(l2(I2(2:end)),flipud(c2(I2(2:end))),'b--','linewidth',2);
h3 = plot(l2(I2),[c2(I2(1)); c2(I2(2:end)) + flipud(c2(I2(2:end)))],'g','linewidth',2);
axis(Ax)
plot(Ax(1:2),[0 0],'k--')
text(7,0.9*Ax(4),'(b)','fontsize',fs)
legend([h2 h3],'<R_{xx,N-n}^b>','<R_{xx,n}^b> + <R_{xx,N-n}^b>','location','north')
legend boxoff
ylabel('<{\itR_{xx}}> [au^2]','fontsize',fs)
xlabel('Lag {\itn} [tu]','fontsize',fs)    
set(gca,'fontsize',fs)

print -f3 -depsc2 ../../figs/corr_fbm_1over2


%% plot auto-correlations for decaying turbulence data from Kang et al. (2003)
close all
clear all
load('../../data/Kang2003.mat');
m = 2^15; % segments length
n = floor(size(u,1)/m); % number of contiguous segments
n = n+n-1; % number of half-overlapping segments

% rearrange in half-overlapping segments
y = repmat(nan,m,n);
for k = 1:n
    y(:,k) = u(1+(k-1)*m/2:(k-1)*m/2+m);
end
clear u

% detrend 
y = detrend(y);

% window
[N,M] = size(y);
w = tukeywin(N,0); % no window
w = repmat(w,1,M);

figure(1)
clf
plot(w(:,1))

% circular unbiased correlation
[c1,l1] = xcorrcirc(w.*y);
c1 = mean(c1,2); % average cross-correlation

% standard biased correlation
[c2,l2] = xcorrbias(w.*y);
c2 = mean(c2,2); % average cross-correlation

% standard unbiased correlation
[c3,l3] = xcorrunbi(w.*y);
c3 = mean(c3,2); % average cross-correlation

% plot auto-correlations
fs = 12; % font size
I1 = find(l1 >= 0);
I2 = find(l2 >= 0);
I3 = find(l3 >= 0);

figure(3)
clf

subplot(211)
plot(l1(I1),c1(I1),'r','linewidth',4)
hold on
plot(l2(I2),c2(I2),'b','linewidth',2)
plot(l3(I3),c3(I3),'c','linewidth',2)
axis tight
Ax = axis;
plot(Ax(1:2),[0 0],'k--')
text(7,0.9*Ax(4),'(a)','fontsize',fs)
legend('<R_{xx,n}^c>','<R_{xx,n}^b>','<R_{xx,n}^u>','location','north')
legend boxoff
ylabel('<{\itR_{xx}}> [au^2]','fontsize',fs)
set(gca,'fontsize',fs)
title('Autocorrelation of decaying turbulence data from Kang et al. (2003)') 

subplot(212)
plot(l1(I1),c1(I1),'r','linewidth',4)
hold on
h1 = plot(l2(I2),c2(I2),'b','linewidth',2);
h2 = plot(l2(I2(2:end)),flipud(c2(I2(2:end))),'b--','linewidth',2);
h3 = plot(l2(I2),[c2(I2(1)); c2(I2(2:end)) + flipud(c2(I2(2:end)))],'g','linewidth',2);
axis(Ax)
plot(Ax(1:2),[0 0],'k--')
text(7,0.9*Ax(4),'(b)','fontsize',fs)
legend([h2 h3],'<R_{xx,N-n}^b>','<R_{xx,n}^b> + <R_{xx,N-n}^b>','location','north')
legend boxoff
ylabel('<{\itR_{xx}}> [au^2]','fontsize',fs)
xlabel('Lag {\itn} [tu]','fontsize',fs)    
set(gca,'fontsize',fs)

print -f3 -depsc2 ../../figs/corr_Kang2003
