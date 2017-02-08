%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotEpsilonProfiles.m
%
% Plot profiles of binned epsilon w/ patch epsilon overplotted. Is epsilon
% only large within patches, or in other places?
%
%------------
% 12/9/16 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% load the patch data ( from Compute_N2_dTdz_ChamProfiles_V2.m)
%load(fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
%   'eq14_cham_patches_diffn2dtdzgamma.mat'))
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/eq14_cham_minOT_25_usetemp_1_patches_diffn2dtdzgamma.mat')


addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/


for cnum=1:50:3100
    
    try
        
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/mat/eq14_' sprintf('%04d',cnum) '.mat'])
        
        clear ig
        ig=find(patches.cnum==cnum);
        
        %
        figure(1);clf
        agutwocolumn(1)
        wysiwyg
        
        
        semilogx(avg.EPSILON,avg.P,'.-')
        axis ij
        grid on
        xlim([1e-11 1e-4])
        hold on
        
        % plot patch epsilons
        semilogx(patches.eps(ig),patches.p1(ig),'d')
        semilogx(patches.eps(ig),patches.p2(ig),'d')
        
        pause(0.5)
    end
end



%% Compare histograms of epsilon from (1) all binned chameleon to (2) just patches

% load binned chameleon data
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')


figure(1);clf
histogram(log10(cham.EPSILON(:)),'Normalization','pdf')
hold on
histogram(log10(patches.eps(:)),'Normalization','pdf')
xlabel('log_{10}[\epsilon]')
ylabel('pdf')
grid on


%% Try comparing time-averaged profiles of epsilon

figure(1);clf
semilogx(nanmean(cham.EPSILON,2),cham.P(:,1))
axis ij
grid on

%% grid patches epsilon so we can plot a profile?

figure(2);clf
semilogx(patches.eps,patches.p1,'.')
grid on
axis ij
%%

pmn=nanmean([patches.p1(:) patches.p2(:)],2);
[Y,I] = sort(pmn,1,'ascend');
eps_sort=patches.eps(I);

figure(1);clf
semilogx(eps_sort,Y,'.','color',0.5*[1 1 1])
axis ij

[xout, hm, NN]=BinProfiles(eps_sort,Y,5,0);
hold on
semilogx(xout,Y,'k','linewidth',3)
grid on
ylabel('P')
xlabel('log_{10}[\epsilon]')
ylim([0 230])
title('patch epsilons')

semilogx(nanmean(cham.EPSILON,2),cham.P(:,1),'r','linewidth',3)

%%

figure(1);clf
semilogx(xout,patches.p1)



%%