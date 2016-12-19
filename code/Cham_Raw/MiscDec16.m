%%

% looking for relation ships between patch gamma and other variables that
% will let us predict patch gamma when we don't measure epsilon

%%
clear ; close all

patch_size_min = 0.25
usetemp = 1

load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

figure(1);clf
%loglog(patches.n3,patches.gam3,'.')
histogram2( log10(patches.n3), log10(patches.gam3), 'Displaystyle','tile');
%histogram2( log10(patches.dtdz3), log10(patches.gam3), 'Displaystyle','tile');
xl=[-4 1]
xl=[-6 1]
xlim(xl)
ylim(xl)

%%

figure(1);clf
scatter(patches.Lt(:),patches.gam3(:),'filled','MarkerFaceAlpha',0.025)
ylim([0 1])
xlim([0 2])

%%

figure(1);clf
scatter(patches.gam_bin(:),patches.gam3(:),'filled','MarkerFaceAlpha',0.025)
xlim([0 1])
ylim([0 1])

%%
figure(1);clf
%loglog(patches.gam_bin(:),patches.gam3(:),'.')
histogram2( log10(patches.gam_bin(:)), log10(patches.gam3(:)),200,'Displaystyle','tile')
xlim([-4 1])
ylim([-4 1])
%xlim([1e-4 1e1])
%ylim([1e-4 1e1])

%%

figure(1);clf
loglog(patches.n3 ./ patches.dtdz3.^2,patches.gam3,'.')
xl=[10^(-2.5) 10^(1.5)]
xlim(xl);ylim(xl)
grid on

figure(2);clf
histogram2( log10(patches.n3 ./ patches.dtdz3.^2),log10(patches.gam3),'Displaystyle','tile')
xl=[-3 2]
xlim(xl);ylim(xl)
grid on

%%

