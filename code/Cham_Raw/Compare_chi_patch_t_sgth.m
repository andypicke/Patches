%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_chi_patch_t_sgth.m
%
% Compare chi-pod estimates for patches, for patches found using
% temperature vs density
%
% 1/23/17 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

whN2dTdz='bulk'
%whN2dTdz='bulk2'
Params.gamma=0.2;
patch_size_min = 0.25

usetemp = 1
sav_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches'
sav_name1 = ['epsilons_N2dTdz_' num2str(whN2dTdz) '_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_chipodmethods.mat'] ; 

clear usetemp
usetemp=0
sav_name2 = ['epsilons_N2dTdz_' num2str(whN2dTdz) '_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_chipodmethods.mat'] ; 


load( fullfile( sav_dir, sav_name1) )
Pt=AllEps; clear AllEps

load( fullfile( sav_dir, sav_name2) )
Pd=AllEps; clear AllEps

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(211)
h1=histogram( real(log10(Pt.chi_patchN2dTdz_constGam(:))),'Normalization','pdf','Edgecolor','none')
hold on
h2=histogram( real(log10(Pd.chi_patchN2dTdz_constGam(:))),h1.BinEdges,'Normalization','pdf','Edgecolor','none')
xlim([-13 -2])
grid on

subplot(212)
h1=histogram( real(log10(Pt.eps_patchN2dTdz_constGam(:))),'Normalization','pdf','Edgecolor','none')
hold on
h2=histogram( real(log10(Pd.eps_patchN2dTdz_constGam(:))),h1.BinEdges,'Normalization','pdf','Edgecolor','none')
xlim([-13 -2])
grid on

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
histogram2( real(log10(Pt.chi_patch(:))), real(log10(Pt.chi_patchN2dTdz_constGam(:))), 50,'DisplayStyle','tile')

subplot(222)
histogram2( real(log10(Pt.eps_patch(:))), real(log10(Pt.eps_patchN2dTdz_constGam(:))), 50,'DisplayStyle','tile')

subplot(223)
histogram2( real(log10(Pd.chi_patch(:))), real(log10(Pd.chi_patchN2dTdz_constGam(:))), 50,'DisplayStyle','tile')

subplot(224)
histogram2( real(log10(Pd.eps_patch(:))), real(log10(Pd.eps_patchN2dTdz_constGam(:))), 50,'DisplayStyle','tile')

%%