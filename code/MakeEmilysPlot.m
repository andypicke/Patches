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

% gamma using all data points
gam_cham=cham.N2 .* cham.CHI ./2 ./ cham.EPSILON ./ (cham.DTDZ.^2);

% find points not in mixed layer or near bottom of profiles
i4=find(cham.P>60 & cham.P<180);

%% Plot N2 vs epsilon, color by gamma ?

close all
figure(1);clf
%$loglog(cham.N2(:),cham.CHI(:),'.')
%scatter( log10(cham.N2(:)), log10(cham.CHI(:)), 1, log10(cham.EPSILON(:)) ) 
scatter( log10(cham.N2(:)), log10(cham.CHI(:)), 0.5, gam_cham(:) ) 
xlabel('log_{10}N^2')
ylabel('log_{10}\chi')
xlim([-6.5 -2.5])
ylim([-14 -2])
grid on
hc=colorbar;
hc.Label.String='\Gamma'
caxis([0 0.2])

%%
figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
print( fullfile( figdir,'N2vsChi_byGamma'), '-dpng')

%%
