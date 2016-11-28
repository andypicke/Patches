%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_Patches_EQ14_raw.m
%
% Examine/plot patches identified in FindPatches_EQ14_Raw.m
%
%
%------------
% 11/09/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% plot patch locations over pcolor of epsilon

clear ; close all

datadir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/'

% load patches
patch_size_min=0.5
usetemp=1
fname=['EQ14_raw_patches_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '.mat']
load( fullfile( datadir, fname) )

% load Chameleon data (1m avg)
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

patch1=new_patch_data;

figure(1);clf
ezpc(cham.castnumber,cham.P,log10(cham.EPSILON))
hold on
plot(patch1(:,1),patch1(:,2),'k.')
hold on
plot(patch1(:,1),patch1(:,3),'r.')
axis ij
caxis([-11 -5])
cb=colorbar;
cb.Label.String='log_{10}\epsilon';
xlim([2000 2500])
xlabel('castnumber')
ylabel('p [db]')

%% Compare to same plot for overturns computed from the averaged/binned data

clear ; close all


datadir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/'

% load patches from raw data
patch_size_min=1
usetemp=1
fname=['EQ14_raw_patches_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '.mat']
load( fullfile( datadir, fname) )
patch1=new_patch_data;

% load patches from 1m avg data
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/EQ14_raw_patches_minOT_1_join_0_sep_50.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/eq14_1m_patches.mat')
patch2=patches_all;
% load Chameleon data
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')


figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(211);
ezpc(cham.castnumber,cham.P,log10(cham.EPSILON))
hold on
plot(patch1(:,1),patch1(:,2),'k.')
hold on
plot(patch1(:,1),patch1(:,3),'r.')
axis ij
caxis([-11 -5])
cb=colorbar;
cb.Label.String='log_{10}\epsilon';
xlim([2000 2500])
xlabel('castnumber')
ylabel('p [db]')

ax2=subplot(212);
ezpc(cham.castnumber,cham.P,log10(cham.EPSILON))
hold on
plot(patch2(:,1),patch2(:,2),'g.')
hold on
plot(patch2(:,1),patch2(:,3),'r.')
axis ij
caxis([-11 -5])
cb=colorbar;
cb.Label.String='log_{10}\epsilon';
xlim([2000 2500])
xlabel('castnumber')
ylabel('p [db]')

linkaxes([ax1 ax2])


%% Plot histogram of patch sizes

clear ; close all

datadir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/'

% load patches from raw data
patch_size_min=1
usetemp=1
fname=['EQ14_raw_patches_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '.mat']
load( fullfile( datadir, fname) )
%patch1=new_patch_data;

figure(1);clf
h1=histogram(new_patch_data(:,4));
xlim([0 20])
freqline(nanmedian(h1.Data))
xlabel('patch size')
grid on

%%

