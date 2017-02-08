%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareProfiles_bin_patch_cham.m
%
% Compare profiles of chi/epsilon for Chameleon (shear probe), chi-pod
% method applied to binned data, and chi-pod method applied to patches
%
% - Chipod method is applied to patches in: xxx.m
%
%
%-------------
% 12/15/16 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% loop through and look at single profiles

clear ; close all

for cnum=60%1:50:31000
    
%    try
        
        % patch N^2,dTdz w/ constant gamma
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches/N2dTdz_bulk_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        avg1=avg;clear avg
        
        % patch N^2,dTdz w/ patch gamma
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches/N2dTdz_bulk_fmax7Hz_respcorr0_fc_99hz_gammaPATCH_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        avg1a=avg;clear avg
        
        % regular chi-pod method on binned data
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        avg2=avg;clear avg
        
        % chameleon (using shear probes)
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/mat/eq14_' sprintf('%04d',cnum) '.mat'])
        
        figure(1);clf
        agutwocolumn(0.7)
        wysiwyg
        
        subplot(121)
        h1=semilogx(avg1.chi1,avg1.P,'o');
        hold on
        h1a=semilogx(avg1a.chi1,avg1a.P,'go');
        h2=semilogx(avg2.chi1,avg2.P,'r-') ;
        h3=semilogx(avg.CHI,avg.P,'k-','linewidth',1);
        axis ij
        grid on
        xlim([1e-11 1e-3])
        ylim([0 200])
        xlabel('\chi')
        ylabel('P')
        set(gca,'XTick',[1e-11 1e-9 1e-7 1e-5 1e-3])
        
        subplot(122)
        h1=semilogx(avg1.eps1,avg1.P,'o');
        hold on
        h1a=semilogx(avg1a.eps1,avg1a.P,'go');
        h2=semilogx(avg2.eps1,avg2.P,'r-') ;
        h3=semilogx(avg.EPSILON,avg.P,'k-','linewidth',1);
        axis ij
        grid on
        xlim([1e-11 1e-4])
        ylim([0 200])
        %legend([h1 h2 h3],'patch','bin','cham')
        xlabel('\epsilon')
        ylabel('P')
        set(gca,'XTick',[1e-11 1e-9 1e-7 1e-5 1e-3])
        
        pause(0.5)
        
%    end % try
    
end % cnum



%% 
% Make a 3X1 pcolor figure:
% 1) chipod method on patches, using n2,dTdz but constant gamma=0.2
% 2) same as (1) but use patch gamma (should get exact answer)
% 3) old chipod method on binned data

clear ; close all

whN2dTdz='bulk'
sav_name = ['epsilons_N2dTdz_' (whN2dTdz) '_chipodmethods.mat'] ; 
sav_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches'
load(fullfile(sav_dir,sav_name))

figure(1);clf
agutwocolumn(1)
wysiwyg

xl=[-8.5 -5]
yl=[-11 -4]

subplot(311)
histogram2(log10(AllEps.eps_bin),log10(AllEps.eps_patchN2dTdz_constGam),100,'DisplayStyle','tile')
hold on
grid on
xlim(xl)
ylim(yl)
xvec=linspace(xl(1),xl(2),100);
loglog(xvec,xvec,'k--')
xlabel('log_{10}[\epsilon Cham]')
ylabel('log_{10}[\epsilon]')
title('patch N^2 and dT/dz, gam=0.2')

subplot(312)
histogram2(log10(AllEps.eps_bin),log10(AllEps.eps_patchN2dTdzGam),100,'DisplayStyle','tile')
%histogram2(log10(eps_patch),log10(eps_2),75,'DisplayStyle','tile')
hold on
grid on
xlim(xl)
ylim(yl)
loglog(xvec,xvec,'k--')
xlabel('log_{10}[\epsilon Cham]')
ylabel('log_{10}[\epsilon]')
title('patch N2,dTdz,and gam')

subplot(313)
histogram2(log10(AllEps.eps_bin),log10(AllEps.eps_chipod_binned),100,'DisplayStyle','tile')
hold on
grid on
xlim(xl)
ylim(yl)
loglog(xvec,xvec,'k--')
xlabel('log_{10}[\epsilon Cham]')
ylabel('log_{10}[\epsilon]')
title('regular chi-pod method, gam=0.2')

% save figure

%print('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/eps_scatter_compare','-dpng')
%print(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/eps_scatter_compare_N2dTdz_' num2str(whN2dTdz) ],'-dpng')

%%

figure(1);clf
agutwocolumn(0.6)
wysiwyg
histogram(log10(eps_bin ./ eps_1 ) )
xlim([-4 4])
freqline(nanmedian(log10(eps_bin ./ eps_1 ) ))

%%

clear ; close all

cnum=1100%:50:31000

%    try

% patch N^2,dTdz w/ constant gamma
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw/zsm0m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
avg1=avg;clear avg

% patch N^2,dTdz w/ patch gamma
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw/zsm0m_fmax7Hz_respcorr0_fc_99hz_gammaPATCH_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
avg1a=avg;clear avg

% regular chi-pod method on binned data
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
avg2=avg;clear avg

% chameleon (using shear probes)
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/mat/eq14_' sprintf('%04d',cnum) '.mat'])


figure(1);clf
agutwocolumn(0.5)
wysiwyg
loglog(avg1.eps_bin(:),avg1.eps1(:),'.')
hold on
loglog(avg1a.eps_bin(:),avg1a.eps1(:),'.')
grid on

%%
figure(1);clf
semilogx(avg.EPSILON,avg.P)
hold on
semilogx(avg1.eps_bin,avg1.P,'o')
semilogx(avg1a.eps_bin,avg1a.P,'d')
grid on
axis ij
%%