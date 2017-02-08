%~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Misc_Jan19.m
%
%
%
%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches/epsilons_N2dTdz_bulk_chipodmethods.mat')
A1=AllEps; clear AllEps

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches/epsilons_N2dTdz_bulk2_chipodmethods.mat')
A2=AllEps; clear AllEps

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(2,2,1)
histogram2( log10(A1.chi_patch(:)), log10(A1.chi_patchN2dTdz_constGam(:)),  100, 'DisplayStyle','Tile')
xlim([-12 -3])
ylim([-12 -3])

subplot(2,2,2)
histogram2( log10(A2.chi_patch(:)), log10(A2.chi_patchN2dTdz_constGam(:)), 100, 'DisplayStyle','Tile')
xlim([-12 -3])
ylim([-12 -3])

subplot(2,2,3)
histogram2( log10(A1.eps_patch(:)), log10(A1.eps_patchN2dTdz_constGam(:)), 100, 'DisplayStyle','Tile')
xlim([-9 -4])
ylim([-9 -4])

subplot(2,2,4)
histogram2( log10(A2.eps_patch(:)), log10(A2.eps_patchN2dTdz_constGam(:)) , 100, 'DisplayStyle','Tile')
xlim([-9 -4])
ylim([-9 -4])

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(2,2,1)
histogram2( log10(A1.chi_bin(:)), log10(A1.chi_patchN2dTdz_constGam(:)),  100, 'DisplayStyle','Tile')
xlim([-12 -3])
ylim([-12 -3])

subplot(2,2,2)
histogram2( log10(A2.chi_bin(:)), log10(A2.chi_patchN2dTdz_constGam(:)), 100, 'DisplayStyle','Tile')
xlim([-12 -3])
ylim([-12 -3])

subplot(2,2,3)
histogram2( log10(A1.eps_bin(:)), log10(A1.eps_patchN2dTdz_constGam(:)), 100, 'DisplayStyle','Tile')
xlim([-9 -4])
ylim([-9 -4])

subplot(2,2,4)
histogram2( log10(A2.eps_bin(:)), log10(A2.eps_patchN2dTdz_constGam(:)) , 100, 'DisplayStyle','Tile')
xlim([-9 -4])
ylim([-9 -4])

%%


figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(2,2,1)
histogram2( log10(A1.chi_patch(:)), log10(A1.chi_patchN2dTdzGam(:)),  100, 'DisplayStyle','Tile')
xlim([-12 -3])
ylim([-12 -3])

subplot(2,2,2)
histogram2( log10(A2.chi_patch(:)), log10(A2.chi_patchN2dTdzGam(:)), 100, 'DisplayStyle','Tile')
xlim([-12 -3])
ylim([-12 -3])

subplot(2,2,3)
histogram2( log10(A1.eps_patch(:)), log10(A1.eps_patchN2dTdzGam(:)), 100, 'DisplayStyle','Tile')
xlim([-9 -4])
ylim([-9 -4])

subplot(2,2,4)
histogram2( log10(A2.eps_patch(:)), log10(A2.eps_patchN2dTdzGam(:)) , 100, 'DisplayStyle','Tile')
xlim([-9 -4])
ylim([-9 -4])