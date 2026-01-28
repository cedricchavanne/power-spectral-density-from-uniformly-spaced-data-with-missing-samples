% test Wiener-Khinchine's theorem on fractional Brownian motion data with
% missing samples (plot spectra)
% before executing this script, you must execute:
% test_fbm
% Cédric Chavanne, 2026-01-27.
% cedric.chavanne@ensta.org

close all
clear all

% parameters
p = [33]; % percent of missing data (must be the same as in do_synthetic_holes)
N = 10000; % number of realizations (must be less or equal than size(y,2))
%N = 100; % number of realizations (must be less or equal than size(y,2))
f0 = 1e-1; % particular frequency at which to plot the PSD histograms
Fs = 10; % font size
ms = 5; % marker size

hurst = {'1over3','1over2'}; % Hurst parameter (name)
hurst_num = {'1/3','1/2'}; % Hurst parameter (numeric)

for j = 1:length(hurst)
    load(['../../data/fbm_',char(hurst(j))])
    
%     % restrict number of realizations
%     y = y(:,1:N);

    % initialization
    [M,N] = size(y);

    for i=1:length(p)
        load(['../../data/fbm_holes_',num2str(p(i)),'percent'])
        load(['../../data/fbm_',char(hurst(j)),'_missing_',num2str(p(i)),'percent_spectra'])

%         % restrict number of realizations
%         w1 = w1(:,1:N);
%         w2 = w2(:,1:N);

        % check total length of holes:
        figure(10)
        clf
        subplot(121)
        hist(M-sum(w1),100)
        title('w1')
        subplot(122)
        hist(M-sum(w2),100)
        title('w2')

        % visualize distribution of holes
        figure(11)
        clf
        imagesc(double(w1))
        colorbar

        figure(12)
        clf
        imagesc(double(w2))
        colorbar

        % replace missing data by nans
        y1 = y;
        y1(find(w1 == 0)) = nan;
        y2 = y;
        y2(find(w2 == 0)) = nan;

        % plot time series example
        x = [0:M-1]'; % independent variable

        figure(1)
        clf
        subplot(511)
        plot(x,y(:,1),'.k','markersize',ms)
        hold on
        plot(x,0*x,'k--')
        axis tight
        Ax = axis;
        text(10,Ax(3)+0.8*(Ax(4)-Ax(3)),'(a)','fontsize',Fs)
        title(['Single realization of fractional Brownian motion with H = ',char(hurst_num(j)),' (linearly-detrended)'],'fontsize',Fs)
        ylabel('{\itx} [au]','fontsize',Fs)
        set(gca,'xticklabel','','fontsize',Fs)

        subplot(512)
        plot(x,w1(:,1),'.k','markersize',ms)
        set(gca,'xticklabel','','ytick',[0 1],'fontsize',Fs)
        axis tight
        Ax = axis;
        text(10,Ax(3)+0.8*(Ax(4)-Ax(3)),'(b)','fontsize',Fs)
        text(50,Ax(3)+0.8*(Ax(4)-Ax(3)),['{\itp} = ',num2str(p(i)),'%'],'fontsize',Fs)
        ylabel('\itw_a','fontsize',Fs)

        subplot(513)
        plot(x,y1(:,1),'.k','markersize',ms)
        hold on
        plot(x,0*x,'k--')
        axis tight
        Ax = axis;
        text(10,Ax(3)+0.8*(Ax(4)-Ax(3)),'(c)','fontsize',Fs)
        ylabel('{\ity_a} [au]','fontsize',Fs)
        set(gca,'xticklabel','','fontsize',Fs)

        subplot(514)
        plot(x,w2(:,1),'.k','markersize',ms)
        set(gca,'xticklabel','','ytick',[0 1],'fontsize',Fs)
        axis tight
        Ax = axis;
        text(10,Ax(3)+0.8*(Ax(4)-Ax(3)),'(d)','fontsize',Fs)
        text(50,Ax(3)+0.8*(Ax(4)-Ax(3)),['{\itp} = ',num2str(p(i)),'%'],'fontsize',Fs)
        ylabel('\itw_b','fontsize',Fs)

        subplot(515)
        plot(x,y2(:,1),'.k','markersize',ms)
        hold on
        plot(x,0*x,'k--')
        axis tight
        Ax = axis;
        text(10,Ax(3)+0.8*(Ax(4)-Ax(3)),'(e)','fontsize',Fs)
        ylabel('{\ity_b} [au]','fontsize',Fs)
        set(gca,'fontsize',Fs)
        xlabel('Time {\itt} [tu]','fontsize',Fs)

        print(1,'-deps2',['../../figs/fbm_',char(hurst(j)),'_series_M',num2str(N),'_',num2str(p(i)),'percent'])

        % average power spectral density of sampling functions
        J = find(f>0);
        f = f(J);
        Pw1 = 2*Pw1(J,:);
        Pw2 = 2*Pw2(J,:);
        % average spectra
        Pw1a = mean(Pw1,2); 
        Pw2a = mean(Pw2,2); 

        figure(2)
        clf
        loglog(f,Pw1a,'k--','linewidth',2)
        hold on
        loglog(f,Pw2a,'k','linewidth',2)
        axis tight
        Ax = axis;
        set(gca,'xtick',[1e-2,1e-1,5e-1],'xticklabel',{'0.01','0.1','0.5'},'ytick',[1e-1,1,1e1],'yticklabel',{'0.1','1','10'},'fontsize',Fs)
        hl = legend('\itw_a','\itw_b','location',[0.7 0.8 0.2 0.1]);
        legend boxoff
        text(0.4*Ax(2),0.1*Ax(4),['{\itp} = ',num2str(p(i)),'%'],'fontsize',Fs)
        text(0.4*Ax(2),0.05*Ax(4),['{\itM} = ',num2str(N)],'fontsize',Fs)
        title(['Sampling functions spectra'],'fontsize',Fs)
        ylabel('<{\itS_{ww}}> [tu]','fontsize',Fs)
        xlabel('{\itf} [tu^{-1}]','fontsize',Fs)

        print(2,'-deps2',['../../figs/sampling_spectra_M',num2str(N),'_',num2str(p(i)),'percent'])

        % average all power spectral densities
        Py = 2*Py(J,:);
        Pya = mean(Py,2); 
        Py1c = 2*Py1c(J,:);
        Py1ca = mean(Py1c,2); 
        Py2c = 2*Py2c(J,:);
        Py2ca = mean(Py2c,2); 
        J = find(fs>0);
        fs = fs(J);
        Py1u = 2*Py1u(J,:);
        Py1ua = mean(Py1u,2);
        Py1Ga = mean(abs(Py1u),2);
        Py2u = 2*Py2u(J,:);
        Py2ua = mean(Py2u,2); 
        Py2Ga = mean(abs(Py2u),2);
        Py1b = 2*Py1b(J,:);
        Py1ba = mean(Py1b,2); 
        Py2b = 2*Py2b(J,:);
        Py2ba = mean(Py2b,2); 
        Py1La = mean(Py1L,2);
        Py2La = mean(Py2L,2);

        % plot data spectra

        figure(3)
        clf
        subplot(121)
        if j==1
            h1 = loglog(f,Pya(2)*(f/f(2)).^(-5/3),'--k');
        else
            h1 = loglog(f,Pya(2)*(f/f(2)).^(-2),'--k');
        end
        hold on
        axis tight
        Ax = axis;
        h3 = loglog(fs,Py1ba,'b','linewidth',2);
        h4 = loglog(fs,Py1ua,'c','linewidth',2);
        h5 = loglog(f,Py1ca,'r','linewidth',2);
        h6 = loglog(fL,Py1La,'--m','linewidth',2);
        h7 = loglog(fs,Py1Ga,'--g','linewidth',2);
        h2 = loglog(f,Pya,'--k','linewidth',2);
        if j==1
            h1 = loglog(f,Pya(2)*(f/f(2)).^(-5/3),'--k');
        else
            h1 = loglog(f,Pya(2)*(f/f(2)).^(-2),'--k');
        end
        if (j == 1 & i == 1 & N == 10000)
            loglog(f0*[1 1],Ax(3:4),'k:')
        end
        text(2e-2*Ax(2),5e-1*Ax(4),'(a)','fontsize',Fs)
        text(2e-2*Ax(2),2e-1*Ax(4),['{\itp} = ',num2str(p(i)),'%'],'fontsize',Fs)
        text(2e-2*Ax(2),1e-1*Ax(4),['{\itM} = ',num2str(N)],'fontsize',Fs)
        set(gca,'xtick',[1e-2,1e-1,5e-1],'xticklabel',{'0.01','0.1','0.5'},'fontsize',Fs-1)
        if j==1
            hl = legend([h1 h2 h3 h4 h5 h6 h7],'{\itf}^{-5/3}','complete','{\itS_{xx}^b}','{\itS_{xx}^u}','{\itS_{xx}^c}','{\itS_{xx}^L}','{\itS_{xx}^G}','location','southwest');
        else
            hl = legend([h1 h2 h3 h4 h5 h6 h7],'{\itf}^{-2}','complete','{\itS_{xx}^b}','{\itS_{xx}^u}','{\itS_{xx}^c}','{\itS_{xx}^L}','{\itS_{xx}^G}','location','southwest');
        end
        set(hl,'fontsize',Fs)
        legend boxoff
        title(['Bernoulli missing ({\itw_a})'],'fontsize',Fs)
        ylabel('<{\itS_{xx}}> [au^2 tu]','fontsize',Fs)
        xlabel('{\itf} [tu^{-1}]','fontsize',Fs)

        subplot(122)
        if j==1
            loglog(f,Pya(2)*(f/f(2)).^(-5/3),'--k')
        else
            loglog(f,Pya(2)*(f/f(2)).^(-2),'--k')
        end
        hold on
        loglog(fs,Py2ba,'b','linewidth',2)
        loglog(fs,Py2ua,'c','linewidth',2)
        loglog(f,Py2ca,'r','linewidth',2)
        loglog(fL,Py2La,'--m','linewidth',2)
        loglog(fs,Py2Ga,'--g','linewidth',2)
        loglog(f,Pya,'--k','linewidth',2)
        if j==1
            loglog(f,Pya(2)*(f/f(2)).^(-5/3),'--k')
        else
            loglog(f,Pya(2)*(f/f(2)).^(-2),'--k')
        end
        axis(Ax)
        text(2e-2*Ax(2),5e-1*Ax(4),'(b)','fontsize',Fs)
        set(gca,'xtick',[1e-2,1e-1,5e-1],'xticklabel',{'0.01','0.1','0.5'},'fontsize',Fs)
        %legend boxoff
        title(['batch-Bernoulli missing ({\itw_b})'],'fontsize',Fs)
        xlabel('{\itf} [tu^{-1}]','fontsize',Fs)

        print(3,'-depsc2',['../../figs/fbm_',char(hurst(j)),'_spectra_M',num2str(N),'_',num2str(p(i)),'percent'])

        % plot psd histograms for a particular frequency
        [junk,k] = min(abs(f-f0));
        [junk,ks] = min(abs(fs-f0));
        [junk,kL] = min(abs(fL-f0));
        pmin = floor(min([min(Py(k,:)) min(Py1b(ks,:)) min(Py1u(ks,:)) min(Py1c(k,:)) min(Py1L(kL,:))]));
        pmax = ceil(max([max(Py(k,:)) max(Py1b(ks,:)) max(Py1u(ks,:)) max(Py1c(k,:)) max(Py1L(kL,:))]));
        pval = pmin:0.05:pmax;
        pmax = 1.5*max([Pya(k) Py1ba(ks) Py1ua(ks) Py1ca(k) Py1La(kL) Py1Ga(ks)]);
        pmin = -pmax/2;
        Nmax = N/10;

        figure(4)
        clf

        subplot(611)
        Ny = histc(Py(k,:),pval);
        bar(pval,Ny,'histc')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',0.7*[1 1 1],'EdgeColor',0.7*[1 1 1])
        set(gca,'xlim',[pmin pmax])
        hold on
        Ax = axis;
        text(Ax(1)+0.1,0.9*Ax(4),'(a)','fontsize',Fs)
        text(Ax(1)+0.1*(Ax(2)-Ax(1)),0.7*Ax(4),['{\itM} = ',num2str(N)],'fontsize',Fs)
        h1 = plot(Pya(k)*[1 1],Ax(3:4),'k','linewidth',4);
        plot([0 0],Ax(3:4),'k--')
        set(gca,'xticklabel','','fontsize',Fs)
        hl = legend(h1,'complete','Location','NorthEast');
        legend boxoff
        title(['fractional Brownian motion with H = ',char(hurst_num(j)),' and Bernoulli missing ({\itw_a})'],'fontsize',Fs-1)

        subplot(612)
        Nb = histc(Py1b(ks,:),pval);
        bar(pval,Nb,'histc')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',0.7*[1 1 1],'EdgeColor',0.7*[1 1 1])
        set(gca,'xlim',[pmin pmax])
        hold on
        Ax = axis;
        text(Ax(1)+0.1,0.9*Ax(4),'(b)','fontsize',Fs)
        text(Ax(1)+0.1*(Ax(2)-Ax(1)),0.7*Ax(4),['{\itp} = ',num2str(p(i)),'%'],'fontsize',Fs)
        plot(Pya(k)*[1 1],Ax(3:4),'k','linewidth',4);
        h2 = plot(Py1ba(ks)*[1 1],Ax(3:4),'b','linewidth',2);
        plot([0 0],Ax(3:4),'k--')
        set(gca,'xticklabel','','fontsize',Fs)
        hl = legend(h2,'{\itS_{xx}^b}','Location','NorthEast');
        legend boxoff

        subplot(613)
        Nu = histc(Py1u(k,:),pval);
        bar(pval,Nu,'histc')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',0.7*[1 1 1],'EdgeColor',0.7*[1 1 1])
        set(gca,'xlim',[pmin pmax])
        hold on
        Ax = axis;
        text(Ax(1)+0.1,0.9*Ax(4),'(c)','fontsize',Fs)
        plot(Pya(k)*[1 1],Ax(3:4),'k','linewidth',4)
        h4 = plot(Py1ua(ks)*[1 1],Ax(3:4),'c','linewidth',2);
        plot([0 0],Ax(3:4),'k--')
        set(gca,'xticklabel','','fontsize',Fs)
        hl = legend(h4,'{\itS_{xx}^u}','Location','NorthEast');
        legend boxoff
        ylabel('Realisations','fontsize',Fs)
  
        subplot(614)
        Nc = histc(Py1c(k,:),pval);
        bar(pval,Nc,'histc')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',0.7*[1 1 1],'EdgeColor',0.7*[1 1 1])
        set(gca,'xlim',[pmin pmax])
        hold on
        Ax = axis;
        text(Ax(1)+0.1,0.9*Ax(4),'(d)','fontsize',Fs)
        plot(Pya(k)*[1 1],Ax(3:4),'k','linewidth',4)
        h4 = plot(Py1ca(k)*[1 1],Ax(3:4),'r','linewidth',2);
        plot([0 0],Ax(3:4),'k--')
        set(gca,'xticklabel','','fontsize',Fs)
        hl = legend(h4,'{\itS_{xx}^c}','Location','NorthEast');
        legend boxoff
  
        subplot(615)
        NL = histc(Py1L(kL,:),pval);
        bar(pval,NL,'histc')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',0.7*[1 1 1],'EdgeColor',0.7*[1 1 1])
        set(gca,'xlim',[pmin pmax])
        hold on
        Ax = axis;
        text(Ax(1)+0.1,0.9*Ax(4),'(e)','fontsize',Fs)
        plot(Pya(k)*[1 1],Ax(3:4),'k','linewidth',4)
        h3 = plot(Py1La(kL)*[1 1],Ax(3:4),'m','linewidth',2);
        plot([0 0],Ax(3:4),'k--')
        set(gca,'xticklabel','','fontsize',Fs);
        hl = legend(h3,'{\itS_{xx}^L}','Location','NorthEast');
        legend boxoff

        subplot(616)
        NG = histc(abs(Py1u(k,:)),pval);
        bar(pval,NG,'histc')
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',0.7*[1 1 1],'EdgeColor',0.7*[1 1 1])
        set(gca,'xlim',[pmin pmax])
        hold on
        Ax = axis;
        text(Ax(1)+0.1,0.9*Ax(4),'(f)','fontsize',Fs)
        plot(Pya(k)*[1 1],Ax(3:4),'k','linewidth',4)
        h4 = plot(Py1Ga(ks)*[1 1],Ax(3:4),'g','linewidth',2);
        plot([0 0],Ax(3:4),'k--')
        set(gca,'fontsize',Fs)
        hl = legend(h4,'{\itS_{xx}^G}','Location','NorthEast');
        legend boxoff
        xlabel(['{\itS_{yy}}(',num2str(f0),') [au^2 tu]'],'fontsize',Fs)

        print(4,'-depsc2',['../../figs/fbm_',char(hurst(j)),'_spectra_hist_M',num2str(N),'_',num2str(p(i)),'percent'])

    end
end
