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

% load patches
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/EQ14_raw_patches_minOT_1_join_0_sep_50.mat')

% load Chameleon data
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

patch1=new_patch_data;

figure(1);clf
ezpc(cham.castnumber,cham.P,log10(cham.EPSILON))
hold on
plot(patch1(:,1),patch1(:,2),'g.')
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

% load patches
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/EQ14_raw_patches_minOT_1_join_0_sep_50.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/eq14_1m_patches.mat')

% load Chameleon data
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

patch1=patches_all;

figure(1);clf
ezpc(cham.castnumber,cham.P,log10(cham.EPSILON))
hold on
plot(patch1(:,1),patch1(:,2),'g.')
hold on
plot(patch1(:,1),patch1(:,3),'r.')
axis ij
caxis([-11 -5])
cb=colorbar;
cb.Label.String='log_{10}\epsilon';
xlim([2000 2500])
xlabel('castnumber')
ylabel('p [db]')



%% Plot histogram of patch sizes

clear ; close all

% load patches
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/EQ14_raw_patches_minOT_1_join_0_sep_50.mat')

figure(1);clf
h1=histogram(new_patch_data(:,4));
xlim([0 20])
freqline(nanmedian(h1.Data))
xlabel('patch size')
grid on

%%

