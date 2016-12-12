%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Misc_Dec_5.m
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot epsilon vs. gamma

clear ; close all

% load the patch data ( from Compute_N2_dTdz_ChamProfiles_V2.m)
load(fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
    'eq14_cham_patches_diffn2dtdzgamma.mat'))

figure(1);clf
%loglog(patches.gam2(:),patches.eps(:),'.')
%histogram2(log10(patches.gam2(:)),log10(patches.eps(:)),100,'DisplayStyle','tile')

xl=[-4 1]

subplot(121)
hs=scatter(log10(patches.gam2(:)),log10(patches.eps(:)),'filled','MarkerFaceAlpha',0.05)
grid on
axis tight
freqline(log10(0.2),'k--')
xlabel('log_{10}[\Gamma]')
ylabel('log_{10}[\epsilon]')
xlim(xl)
title('patches')

subplot(122)
hs=scatter(log10(patches.gam_bin(:)),log10(patches.eps(:)),'filled','MarkerFaceAlpha',0.05)
grid on
axis tight
freqline(log10(0.2),'k--')
xlabel('log_{10}[\Gamma]')
ylabel('log_{10}[\epsilon]')
xlim(xl)
title('bin')

figdir='/Users/Andy/Cruises_Research/ChiPod/Analyses/Patch_n2_dTdz'
print( fullfile( figdir, ['eq14_cham_gamma_epsVsgam_scatter'] ), '-dpng' )


%%