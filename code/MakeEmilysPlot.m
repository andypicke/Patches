%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% MakeEmilysPlot.m
%
% Can you make a scatter plot with N2 (or dTdz) on one axis, log10(chi) on a second axis, and color by the gamma (or log10 gamma) needed to have eps_chi=eps? Sorry I am stuck on the stratification piece. (Patches clumping so strongly at 1 m does seem like it might be problematic.)
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

% There's too many points and they are being overplotted, try taking a
% random subsample to plot
ir=randi(numel(n2),[round(0.15*numel(n2)),1] );
%ir=randi(numel(n2),[40000,1] );

close all
figure(1);clf
set(gca,'DefaultAxesFontsize',15)
h=scatter( log10(n2(:)), log10(chi(:)), 10, log10(gam_cham(:)),'o', 'filled' ) 
%h=scatter( log10(n2(ir)), log10(chi(ir)), 7, log10(gam_cham(ir)),'o', 'filled' ) 
h.MarkerFaceAlpha=0.2
xlabel('log_{10}N^2','fontsize',16)
ylabel('log_{10}\chi','fontsize',16)
xlim([-6 -2.5])
ylim([-12 -4])
grid on
hc=colorbar;
hc.Label.String='log_{10}\Gamma';
hc.Label.FontSize=14;
caxis([-2.5 -0.5])


%%

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
print( fullfile( figdir,'N2vsChi_byGamma'), '-dpng')

%% Plot N2 vs Epsilon, colored by gamma
close all
figure(1);clf
h=scatter( log10(n2(:)), log10(eps(:)), 10, log10(gam_cham(:)), 'o', 'filled' ) 
%h=scatter( log10(n2(ir)), log10(eps(ir)), 7, log10(gam_cham(ir)), 'o', 'filled' ) 
h.MarkerFaceAlpha=0.2
xlabel('log_{10}N^2','fontsize',16)
ylabel('log_{10}\epsilon','fontsize',16)
xlim([-6 -2.5])
ylim([-8 -4])
grid on
hc=colorbar;
hc.Label.String='log_{10}\Gamma';
hc.Label.FontSize=14;
%caxis([-2.25 -.85])
caxis([-2.5 -0.5])
%colormap(jet)

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
print( fullfile( figdir,'N2vsEps_byGamma'), '-dpng')
%%