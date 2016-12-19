%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% MakePlotsGammaAnalysis.m
%
% Used to be part of Compute_N2_dTdz_ChamProfiles_V2.m
%
%-----------
% 12/12/16 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Make plots

clear ; close all

patch_size_min = 0.15
usetemp = 1
saveplots = 1

load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

%figdir = '/Users/Andy/Cruises_Research/ChiPod/Analyses/Patch_n2_dTdz'
figdir=['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) ]
ChkMkDir(figdir)
%%

figure(1);clf
agutwocolumn(0.5)
wysiwyg

h1 = histogram(log10(patches.gam1(:)),200,'EdgeColor','none') ;
hold on
histogram(log10(patches.gam2(:)),h1.BinEdges,'EdgeColor','none') ;
histogram(real(log10(patches.gam3(:))),h1.BinEdges,'EdgeColor','none') ;
histogram(real(log10(patches.gam4(:))),h1.BinEdges,'EdgeColor','none') ;
xlim([-3 1])
grid on
xlabel('log_{10}\Gamma')
ylabel('count')

iz1 = find( patches.gam1>0 & patches.gam1<1 ) ;
iz2 = find( patches.gam2>0 & patches.gam2<1 ) ;
iz3 = find( patches.gam3>0 & patches.gam3<1 ) ;
iz4 = find( patches.gam4>0 & patches.gam4<1 ) ;

figure(2);clf
agutwocolumn(0.6)
wysiwyg

h1 = histogram( patches.gam1(iz1), 50 , 'Normalization', 'pdf') ;
hold on
h2 = histogram( patches.gam2(iz2) , h1.BinEdges , 'Normalization', 'pdf') ;
h3 = histogram( patches.gam3(iz2) , h1.BinEdges , 'Normalization', 'pdf') ;
h4 = histogram( patches.gam4(iz2) , h1.BinEdges , 'Normalization', 'pdf') ;
freqline(nanmedian(patches.gam1))
freqline(nanmedian(patches.gam2))
freqline(nanmedian(patches.gam3))
freqline(nanmedian(patches.gam4))
text(nanmedian(patches.gam1),6,'\Gamma 1')
text(nanmedian(patches.gam2),5.5,'\Gamma 2')
text(nanmedian(patches.gam3),5,'\Gamma 3')
text(nanmedian(patches.gam4),4.5,'\Gamma 4')
legend([h1 h2 h3 h4],'\Gamma 1','\Gamma 2','\Gamma 3','\Gamma 4')
xlim([0 0.5])
grid on
xlabel('\Gamma')
ylabel('count')

if saveplots==1
    print( fullfile( figdir, ['eq14_cham_patch_gammas'] ), '-dpng' )
end


%% Plot 2D hist/scatter of different dt/dz methods

figure(1);clf
agutwocolumn(1)
wysiwyg
set(gcf,'defaultaxesfontsize',15)

subplot(211)
%loglog(abs(dtdz1(:)),abs(dtdz2(:)),'.')
histogram2(log10(abs(patches.dtdz1(:))),log10(abs(patches.dtdz2(:))),50,'displaystyle','tile')
xlabel('log_{10}[dtdz1]')
ylabel('log_{10}[dtdz2]')
%grid on
grid on

subplot(212)
%loglog(abs(dtdz1(:)),abs(dtdz2(:)),'.')
histogram2(log10(abs(patches.dtdz1(:))),log10(abs(patches.dtdz3(:))),50,'displaystyle','tile')
xlabel('log_{10}[dtdz1]')
ylabel('log_{10}[dtdz3]')
%grid on
grid on

if saveplots==1
print( fullfile( figdir, ['eq14_cham_patch_dTdzs_scatter'] ), '-dpng' )
end

%% Make histogram of different dt/dzs

figure(1);clf
agutwocolumn(0.7)
wysiwyg
set(gcf,'defaultaxesfontsize',14)

h1 = histogram(log10(abs(patches.dtdz1(:))),'DisplayStyle','stair') ;%,'edgecolor','none');
hold on
h2 = histogram(log10(abs(patches.dtdz2(:))),h1.BinEdges,'DisplayStyle','stair') ;%,'edgecolor','none');
h3 = histogram(log10(abs(patches.dtdz3(:))),h1.BinEdges,'DisplayStyle','stair') ;%,'edgecolor','none');
grid on
xlim([-4 -0.5])
xlabel('log_{10}[dT/dz]')
ylabel('count')
legend([h1 h2 h3],'dtdz1','dtdz2','dtdz3')

if saveplots==1
print( fullfile( figdir, ['eq14_cham_patch_dTdzs'] ), '-dpng' )
end

%% Plot 2D histograms of different N2 methods

figure(2);clf
agutwocolumn(1)
wysiwyg
set(gcf,'defaultaxesfontsize',15)

xl=[-5.5 -3]
yl=[-6.75 -2.5]
m=3
n=1

subplot(m,n,1)
%loglog(n1(:),nb(:),'.')
histogram2(real(log10(patches.n1(:))),real(log10(patches.nb(:))),50,'displaystyle','tile')
grid on
xlim(xl)
ylim(yl)

xlabel('log_{10}[N^2_1]')
ylabel('log_{10}[N^2_2]')

subplot(m,n,2)
%loglog(n1(:),n3(:),'.')
histogram2(real(log10(patches.n1(:))),real(log10(patches.n3(:))),50,'displaystyle','tile')
grid on
xlabel('log_{10}[N^2_1]')
ylabel('log_{10}[N^2_3]')
xlim(xl)
ylim(yl)

subplot(m,n,3)
%loglog(n1(:),n3(:),'.')
histogram2(real(log10(patches.n1(:))),real(log10(patches.n4(:))),50,'displaystyle','tile')
grid on
xlabel('log_{10}[N^2_1]')
ylabel('log_{10}[N^2_4]')
xlim([xl])
ylim([yl])

if saveplots==1
print( fullfile( figdir, ['eq14_cham_patch_N2s_scatter'] ), '-dpng' )
end

%% Plot histogram of different N2s

figure(2);clf
agutwocolumn(0.7)
wysiwyg
set(gcf,'defaultaxesfontsize',14)

h1 = histogram(real(log10(patches.n1(:)))) ;
hold on
h2 = histogram(real(log10(patches.nb(:))),h1.BinEdges) ;
h3 = histogram(real(log10(patches.n3(:))),h1.BinEdges) ;
h4 = histogram(real(log10(patches.n4(:))),h1.BinEdges) ;
freqline(nanmedian(h1.Data))
freqline(nanmedian(h2.Data))
freqline(nanmedian(h3.Data))
freqline(nanmedian(h4.Data))
grid on
xlim([-6 -2.5])
xlabel('log_{10}[N^2]')
ylabel('count')
legend([h1 h2 h3 h4],'1  ','2 ','3 ','4 ','location','best')

if saveplots==1
print( fullfile( figdir, ['eq14_cham_patch_N2s'] ), '-dpng' )
end

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1 = subplot(221);
h1 = histogram(real(log10(patches.nb(:))),'Normalization','pdf','edgecolor','none');
hold on
h2 = histogram(real(log10(patches.n2_bin(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
xlim([-6.5 -2.5])
grid on
xlabel('log_{10}[N^2]')
legend([h1 h2],'patch','bin','location','best')
freqline(nanmedian(h1.Data),'b--')
freqline(nanmedian(h2.Data),'r--')

ax2 = subplot(222);
h1 = histogram(real(log10(patches.dtdz2(:))),'Normalization','pdf','edgecolor','none');
hold on
h2 = histogram(real(log10(patches.dtdz_bin(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
xlim([-5 0])
grid on
xlabel('log_{10}[dT/dz]')
legend([h1 h2],'patch','bin','location','best')
freqline(nanmedian(h1.Data),'b--')
freqline(nanmedian(h2.Data),'r--')

ax3 = subplot(223);
h1 = histogram(real(log10(patches.chi(:))),'Normalization','pdf','edgecolor','none');
hold on
h2 = histogram(real(log10(patches.chi_bin(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
xlim([-12 -3])
grid on
xlabel('log_{10}[\chi]')
legend([h1 h2],'patch','bin','location','best')
freqline(nanmedian(h1.Data),'b--')
freqline(nanmedian(h2.Data),'r--')

ax4 = subplot(224);
h1 = histogram(real(log10(patches.eps(:))),'Normalization','pdf','edgecolor','none');
hold on
h2 = histogram(real(log10(patches.eps_bin(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
%xlim([-12 -3])
grid on
xlabel('log_{10}[\epsilon]')
legend([h1 h2],'patch','bin','location','best')
freqline(nanmedian(h1.Data),'b--')
freqline(nanmedian(h2.Data),'r--')

if saveplots==1
print( fullfile( figdir, ['eq14_cham_gamma_binVspatch_hists'] ), '-dpng' )
end

%%

ib=find(patches.gam_bin>1);

figure(1);clf
agutwocolumn(0.5)
wysiwyg
%scatter(log10(patches.gam_bin(:)),log10(patches.gam3(:)),'o','filled','MarkerFaceAlpha',0.05)
scatter(patches.gam_bin(:), patches.gam3(:),'o','filled','MarkerFaceAlpha',0.1) ; xlim([0 0.5]) ; ylim([0 0.5])
%xlim([-3 1]);ylim([-3 1])
xlabel('\Gamma bin')
ylabel('\Gamma patch')

% ig=find(~isnan(patches.gam_bin) & ~isnan(patches.gam2));
% P=polyfit(patches.gam_bin(ig),patches.gam3(ig),1);
% hold on
% plot(0:0.001:1,polyval(P,0:0.001:1),'r','linewidth',2)
grid on

if saveplots==1
print( fullfile( figdir, ['eq14_cham_gamma_binVspatch_scatter'] ), '-dpng' )
end


%%

figure(2);clf
agutwocolumn(0.6)
wysiwyg

h1 = histogram(log10(patches.gam_bin(:)),'Normalization','pdf');
hold on
h2 = histogram(log10(patches.gam2(:)),h1.BinEdges,'Normalization','pdf');
freqline(nanmedian(h1.Data),'b--')
freqline(nanmedian(h2.Data),'r--')
xlim([-4 2])
grid on
xlabel('log_{10}[\Gamma]','fontsize',16)
ylabel('pdf','fontsize',16)
legend([h1 h2],'bin','patch')
nanmedian(patches.gam_bin(:))
nanmedian(patches.gam3(:))
hf=freqline(log10(0.2),'k--')
set(hf,'linewidth',2)

if saveplots==1
print( fullfile( figdir, ['eq14_cham_LOGgamma_binVspatch'] ), '-dpng' )
end


%% Scatter plot binned vs patch values

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
%loglog(patches.chi_bin(:),patches.chi(:),'.')
%scatter( log10(patches.n2_bin(:)), log10(patches.nb(:)),'o','filled','MarkerFaceAlpha',0.05) ; 
histogram2(real(log10(patches.n2_bin(:))), real(log10(patches.nb(:))), 100 ,'DisplayStyle','Tile')
grid on
%xl=[1e-11 1e-4]
xl=[-6.5 -2.5]
xlim(xl)
ylim(xl)
xlabel('N^2 binned','fontsize',16)
ylabel('N^2 patch' ,'fontsize',16)
xvec=linspace(xl(1),xl(2),100);
hold on
plot(xvec,xvec,'k--','linewidth',2)


subplot(222)
%loglog(patches.chi_bin(:),patches.chi(:),'.')
%scatter( log10(patches.dtdz_bin(:)), log10(patches.dtdz2(:)),'o','filled','MarkerFaceAlpha',0.05) ; 
histogram2(real(log10(patches.dtdz_bin(:))), real(log10(patches.dtdz2(:))), 100 ,'DisplayStyle','Tile')
grid on
%xl=[1e-11 1e-4]
xl=[-4.5 0]
xlim(xl)
ylim(xl)
xlabel('dT/dz binned','fontsize',16)
ylabel('dT/dz patch' ,'fontsize',16)
xvec=linspace(xl(1),xl(2),100);
hold on
plot(xvec,xvec,'k--','linewidth',2)

subplot(223)
%loglog(patches.chi_bin(:),patches.chi(:),'.')
%scatter( log10(patches.chi_bin(:)), log10(patches.chi(:)),'o','filled','MarkerFaceAlpha',0.05) ; 
histogram2(real(log10(patches.chi_bin(:))), real(log10(patches.chi(:))), 100 ,'DisplayStyle','Tile')
grid on
%xl=[1e-11 1e-4]
xl=[-12 -3]
xlim(xl)
ylim(xl)
xlabel('\chi binned','fontsize',16)
ylabel('\chi patch' ,'fontsize',16)

subplot(224)
%loglog(patches.chi_bin(:),patches.chi(:),'.')
%scatter( log10(patches.eps_bin(:)), log10(patches.eps(:)),'o','filled','MarkerFaceAlpha',0.05) ; 
histogram2(real(log10(patches.eps_bin(:))), real(log10(patches.eps(:))), 100 ,'DisplayStyle','Tile')
grid on
xl=[-8.5 -5]
xlim(xl)
ylim(xl)
xlabel('\epsilon binned','fontsize',16)
ylabel('\epsilon patch' ,'fontsize',16)
xvec=linspace(xl(1),xl(2),100);
hold on
plot(xvec,xvec,'k--','linewidth',2)

if saveplots==1
    print( fullfile( figdir, ['eq14_binVspatch_N2dTdzChiEps'] ), '-dpng' )
end

%%


figure(1);clf
agutwocolumn(0.6)
wysiwyg

histogram2(real(log10(patches.gam_bin)), real(log10(patches.gam2)), 200 ,'DisplayStyle','Tile')
grid on
xl=[-3.5 1]
xlim(xl)
ylim(xl)
xlabel('\Gamma binned','fontsize',16)
ylabel('\Gamma patch' ,'fontsize',16)
xvec=linspace(xl(1),xl(2),100);
hold on
plot(xvec,xvec,'k--','linewidth',2)
%plot(log10(0.2),log10(0.2),'rx','markersize',15,'linewidth',3)
vline(log10(0.2),'r')
hline(log10(0.2),'r')

if saveplots==1
    print( fullfile( figdir, ['eq14_binVspatch_gamma'] ), '-dpng' )
end

%%