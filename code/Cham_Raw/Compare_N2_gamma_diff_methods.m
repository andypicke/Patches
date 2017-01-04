%~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_N2_gamma_diff_methods.m
%
% This code makes plots of (1) N2/dTdz^2 and (2) gamma vs eps for two
% different methods of computing N2 and dTdz (the poly fit method, and
% Bill's 'bulk' method.
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

patch_size_min = 0.25
usetemp = 1

load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )


%%

ratio1=patches.nb./(patches.dtdz2.^2);
ratio2=patches.n3./(patches.dtdz3.^2);

figure(1);clf
%loglog(ratio2,ratio1,'.')
histogram2(real(log10(ratio2)), real(log10(ratio1)),100,'DisplayStyle','tile')
grid on
xlabel('log_{10}[N^2/dTdz^2] bulk')
ylabel('log_{10}[N^2/dTdz^2] line')
xlim([-2 2])
ylim([-2 3])
xvec=linspace(-2,2,100);
hold on
plot(xvec,xvec,'k--')
title('ratio of [N^2/dTdz^2] for different methods')
%%

print(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/compare_N2dTdz_ratios'],'-dpng')

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(211)
histogram2(real(log10(patches.eps)), real(log10(patches.gam2)),100,'DisplayStyle','tile')
%histogram2(real(log10(patches.eps)), real(log10(patches.gam3)),100,'DisplayStyle','tile')
ylim([-4 1])
hline(log10(0.2),'k--')
xlabel('log_{10}\epsilon','fontsize',16)
ylabel('log_{10}\Gamma','fontsize',16)
title('straight-line N^2,dTdz')

subplot(212)
%histogram2(real(log10(patches.eps)), real(log10(patches.gam2)),100,'DisplayStyle','tile')
histogram2(real(log10(patches.eps)), real(log10(patches.gam3)),100,'DisplayStyle','tile')
ylim([-4 1])
hline(log10(0.2),'k--')
xlabel('log_{10}\epsilon','fontsize',16)
ylabel('log_{10}\Gamma','fontsize',16)
title('bulk N2,dTdz')

%%

print(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/gam_vs_eps_diff_methods'],'-dpng')

%%

