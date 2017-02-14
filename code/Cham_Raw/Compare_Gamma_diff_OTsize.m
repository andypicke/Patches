%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_Gamma_diff_OTsize.m
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

usetemp = 1

patch_size_min = 0.15
load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']))
patches1=patches; clear patches

patch_size_min = 0.25
load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']))
patches2=patches; clear patches

patch_size_min = 0.5
load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']))
patches3=patches; clear patches

figure(1);clf
histogram( log10(patches1.gam_bulk(:)), 'Normalization','pdf')
hold on
histogram( log10(patches2.gam_bulk(:)), 'Normalization','pdf')
histogram( log10(patches3.gam_bulk(:)), 'Normalization','pdf')
xlim([-4 0.5])
freqline(log10(0.2))
grid on
xlabel('log_{10}\Gamma')

%%

figure(1);clf
histogram( log10(patches1.gam_line(:)), 'Normalization','pdf')
hold on
histogram( log10(patches2.gam_line(:)), 'Normalization','pdf')
histogram( log10(patches3.gam_line(:)), 'Normalization','pdf')
xlim([-3.5 0.5])
freqline(log10(0.2))
grid on
xlabel('log_{10}\Gamma')

%%


figure(1);clf
histogram( log10(patches1.gam4(:)), 'Normalization','pdf')
hold on
histogram( log10(patches2.gam4(:)), 'Normalization','pdf')
histogram( log10(patches3.gam4(:)), 'Normalization','pdf')
xlim([-3.5 0.5])
freqline(log10(0.2))
grid on
xlabel('log_{10}\Gamma')

%%