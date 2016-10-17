%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% MakeEmilysPlot.m
%
% Can you make a scatter plot with N2 (or dTdz) on one axis, log10(chi) on a second axis, and color by the gamma (or log10 gamma) needed to have eps_chi=eps? Sorry I am stuck on the stratification piece. (Patches clumping so strongly at 1 m does seem like it might be problematic.)
%
% * Note: I think the scatter plots are misleading, because there are too
% many points. See the 2D histograms further down, these should be more
% accurate.
%
%
%-------------------
%  10/05/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% load EQ14 data
dir_base='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed';
clear cham
load( fullfile( dir_base, '/Cstar=0_01366/sum/eq14_sum.mat') )

% throw out epsilon below noise floor
ib=find(log10(cham.EPSILON)<-8);
cham.EPSILON(ib)=nan;

% find points not in mixed layer or near bottom of profiles
ig=find(cham.P>60 & cham.P<180);

n2=cham.N2(ig);
dtdz=cham.DTDZ_RHOORDER(ig);
chi=cham.CHI(ig);
eps=cham.EPSILON(ig);

% gamma using all data points
%gam_cham=cham.N2 .* cham.CHI ./2 ./ cham.EPSILON ./ (cham.DTDZ.^2);
gam_cham=n2 .* chi ./2 ./ eps ./ (dtdz.^2);
% 



%% Plot N2 vs chi, color by gamma ?
% 
% % There's too many points and they are being overplotted, try taking a
% % random subsample to plot
% ir=randi(numel(n2),[round(0.05*numel(n2)),1] );
% %ir=randi(numel(n2),[40000,1] );
% 
% close all
% figure(1);clf
% set(gca,'DefaultAxesFontsize',15)
% %h=scatter( log10(n2(:)), log10(chi(:)), 10, log10(gam_cham(:)),'o', 'filled' ) 
% h=scatter( log10(n2(ir)), log10(chi(ir)), 10, log10(gam_cham(ir)),'o', 'filled' ) 
% 
% % Make transparent (note only works in vers >= 2015b)
% h.MarkerFaceAlpha=0.2
% 
% xlabel('log_{10}N^2','fontsize',16)
% ylabel('log_{10}\chi','fontsize',16)
% xlim([-6 -2.5])
% ylim([-12 -4])
% grid on
% hc=colorbar;
% hc.Label.String='log_{10}\Gamma';
% hc.Label.FontSize=14;
% caxis([-2.5 -0.5])
% 
% 
% %%
% 
% figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
% print( fullfile( figdir,'N2vsChi_byGamma'), '-dpng')
% 
% %% Plot N2 vs Epsilon, colored by gamma
% close all
% figure(1);clf
% %h=scatter( log10(n2(:)), log10(eps(:)), 10, log10(gam_cham(:)), 'o', 'filled' ) 
% h=scatter( log10(n2(ir)), log10(eps(ir)), 10, log10(gam_cham(ir)), 'o', 'filled' ) 
% h.MarkerFaceAlpha=0.2
% xlabel('log_{10}N^2','fontsize',16)
% ylabel('log_{10}\epsilon','fontsize',16)
% xlim([-6 -2.5])
% ylim([-8 -4])
% grid on
% hc=colorbar;
% hc.Label.String='log_{10}\Gamma';
% hc.Label.FontSize=14;
% %caxis([-2.25 -.85])
% caxis([-2.5 -0.5])
% %colormap(jet)
% 
% figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
% print( fullfile( figdir,'N2vsEps_byGamma'), '-dpng')
% 
% %%
% 
% figure(3);clf
% h=scatter(log10(eps),log10(gam_cham),'.','MarkerFaceAlpha',0.25)
% xlabel('\epsilon','fontsize',16)
% ylabel('\Gamma','fontsize',16)
% grid on
% ylim([-5 1])
% %%
% figure(5);clf
% histogram2(log10(eps),log10(gam_cham),200,'DisplayStyle','tile')
% xlabel('log_{10}\epsilon','fontsize',16)
% ylabel('log_{10}\Gamma','fontsize',16)
% colorbar
% ylim([-5 2])
% xlim([-8 -4])
% title('2D hist')
% 
% %%
% 
% figure(4);clf
% h=scatter(log10(chi),log10(gam_cham),'.','MarkerFaceAlpha',0.25)
% xlabel('\chi','fontsize',16)
% ylabel('\Gamma','fontsize',16)
% grid on
% ylim([-5 1])
% 
% 
% %%
% figure(5);clf
% histogram2(log10(chi),log10(gam_cham),200,'DisplayStyle','tile')
% xlabel('log_{10}\chi','fontsize',16)
% ylabel('log_{10}\Gamma','fontsize',16)
% colorbar
% ylim([-6 3])
% xlim([-13 -3])
% title('2D hist')
% 
% %%
% 
% figure(5);clf
% histogram2(real(log10(n2)),log10(gam_cham),200,'DisplayStyle','tile')
% xlabel('log_{10}N^2','fontsize',16)
% ylabel('log_{10}\Gamma','fontsize',16)
% colorbar
% ylim([-6 3])
% xlim([-5.5 -3])
% title('2D hist')
% 
% %%
% 
% figure(5);clf
% histogram2(real(log10(dtdz)),log10(gam_cham),200,'DisplayStyle','tile')
% xlabel('log_{10}dT/dz','fontsize',16)
% ylabel('log_{10}\Gamma','fontsize',16)
% colorbar
% ylim([-6 3])
% xlim([-4 -0])
% title('2D hist')
% 

%% Plot 2D histograms of gamma vs each variable

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(221);
histogram2(real(log10(n2)),log10(gam_cham),200,'DisplayStyle','tile')
xlabel('log_{10}N^2','fontsize',16)
ylabel('log_{10}\Gamma','fontsize',16)
colorbar
ylim([-6 3])
xlim([-6 -2])

ax2=subplot(222);
histogram2(real(log10(dtdz)),log10(gam_cham),200,'DisplayStyle','tile')
xlabel('log_{10}dT/dz','fontsize',16)
ylabel('log_{10}\Gamma','fontsize',16)
colorbar
ylim([-6 3])
xlim([-4 -0])

ax3=subplot(223);
histogram2(real(log10(chi)),log10(gam_cham),200,'DisplayStyle','tile')
xlabel('log_{10}dT/dz','fontsize',16)
ylabel('log_{10}\Gamma','fontsize',16)
colorbar
ylim([-6 3])
xlim([-13 -3])

ax4=subplot(224);
histogram2(real(log10(eps)),log10(gam_cham),200,'DisplayStyle','tile')
xlabel('log_{10}dT/dz','fontsize',16)
ylabel('log_{10}\Gamma','fontsize',16)
colorbar
ylim([-6 3])
xlim([-8 -3])

linkaxes([ax1 ax2 ax3 ax4],'y')

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
print( fullfile( figdir,'2X2_2dhist_Gamma'), '-dpng')

%% try binning gamma by N^2 and chi

% brute force method (there has to be some better way of doing this?)

dchi=logspace(-12,-4,100); % make different lengths so we can't confuse axes
dN=logspace(-6,-2.5,95);
mean_gamma=nan*ones(length(dchi)-1,length(dN)-1);
median_gamma=mean_gamma;
for m=1:length(dchi)-1;
        for n=1:length(dN-1)-1;
            clear ind
                ind=find( n2>dN(n) & n2<dN(n+1) & chi>dchi(m) & chi<dchi(m+1) );
                mean_gamma(m,n)=nanmean(gam_cham(ind));
                median_gamma(m,n)=nanmedian(gam_cham(ind));
        end
end

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(211)
pcolor(log10(dN(1:end-1)),log10(dchi(1:end-1)),mean_gamma)
%pcolor(log10(dN(1:end-1)),log10(dchi(1:end-1)),median_gamma)
shading flat
caxis([0 0.1])
cb=colorbar;
cb.Label.String='\Gamma';
ylabel('log_{10}\chi','fontsize',16)
xlabel('log_{10}N^2','fontsize',16)
grid on
title('mean gamma') 

subplot(212)
%pcolor(log10(dN(1:end-1)),log10(dchi(1:end-1)),mean_gamma)
pcolor(log10(dN(1:end-1)),log10(dchi(1:end-1)),median_gamma)
shading flat
caxis([0 0.11])
cb=colorbar;
cb.Label.String='\Gamma';
ylabel('log_{10}\chi','fontsize',16)
xlabel('log_{10}N^2','fontsize',16)
grid on
title('median gamma')

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
print( fullfile( figdir,'Gamma_binnedBy_N2_chi'), '-dpng')
%%