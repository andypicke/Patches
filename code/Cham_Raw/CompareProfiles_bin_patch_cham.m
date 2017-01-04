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

for cnum=1:50:31000
    
    try
        
        % patch N^2,dTdz w/ constant gamma
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw/zsm0m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        avg1=avg;clear avg
        
        % patch N^2,dTdz w/ patch gamma
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw/zsm0m_fmax7Hz_respcorr0_fc_99hz_gammaPATCH_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
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
        
    end % try
    
end % cnum


%% Combine results from all profiles (so we can make scatter plot etc) and save

clear ; close all

eps_patchN2dTdz_constGam = [] ;
eps_patchN2dTdzGam    = [] ;
eps_bin  = [] ;
eps_patch= [] ;
eps_chipod_binned    = [] ;

whN2dTdz=2
dir1=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw/')
%,...
 %       ['N2dTdz_' num2str(whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gammaPATCH_nfft_' num2str(Params.nfft)]);

for cnum=1:3100
    
    try
        
        % patch N^2,dTdz w/ constant gamma
        %load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw/zsm0m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        load( fullfile( dir1, ['N2dTdz_' num2str(whN2dTdz) '_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128'],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        avg_patchN2dTdz_constGam=avg;clear avg
        
        % patch N^2,dTdz w/ patch gamma
        load( fullfile( dir1, ['N2dTdz_' num2str(whN2dTdz) '_fmax7Hz_respcorr0_fc_99hz_gammaPATCH_nfft_128'],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        %load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw/zsm0m_fmax7Hz_respcorr0_fc_99hz_gammaPATCH_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        %load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw/zsm0m_fmax7Hz_respcorr1_fc_32hz_gammaPATCH_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        avg_patchN2dTdzGam=avg;clear avg
                
        % regular chi-pod method on binned data
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        %load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma5_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        avg_chipod_binned=avg;clear avg   
        
        % get these values at patch locations
        eps3=nan*ones(size(avg_patchN2dTdz_constGam.eps1));
        for ip=1:length(avg_patchN2dTdz_constGam.P)
           [val,I]=nanmin(abs(avg_chipod_binned.P-avg_patchN2dTdz_constGam.P(ip)));
           eps3(ip)=avg_chipod_binned.eps1(I);
        end
        
        % regular chi-pod method epsilons (10m smooth, gam=0.2 etc)
        eps_chipod_binned  =[eps_chipod_binned   ; eps3(:)        ];

        eps_patchN2dTdz_constGam = [eps_patchN2dTdz_constGam ; avg_patchN2dTdz_constGam.eps1(:) ];
        eps_patchN2dTdzGam       = [eps_patchN2dTdzGam       ; avg_patchN2dTdzGam.eps1(:)       ];
    
        % chameleon epsilon binned data at patch locations
        eps_bin=[eps_bin ; avg_patchN2dTdz_constGam.eps_bin(:)];
                
        % chameleon epsilon values computed over patches only
        eps_patch=[eps_patch ; avg_patchN2dTdzGam.eps_patch(:) ] ;
        
    end % try
    
end % cnum

%

ib=find(log10(eps_bin)<-8.5);
eps_bin(ib)=nan;

AllEps=struct('eps_chipod_binned',eps_chipod_binned,'eps_patchN2dTdz_constGam',...
    eps_patchN2dTdz_constGam,'eps_patchN2dTdzGam',eps_patchN2dTdzGam,...
    'eps_bin',eps_bin,'eps_patch',eps_patch)
AllEps.MakeInfo=['Made ' datestr(now) ' w/ CompareProfiles_bin_patch_cham.m']
% save data
sav_name = ['epsilons_N2dTdz_' num2str(whN2dTdz) '_chipodmethods.mat'] ; 
sav_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw'
save(fullfile(sav_dir,sav_name),'AllEps')

%%

clear ; close all

whN2dTdz=3
sav_name = ['epsilons_N2dTdz_' num2str(whN2dTdz) '_chipodmethods.mat'] ; 
sav_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/Cham_Raw'
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
print(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/eps_scatter_compare_N2dTdz_' num2str(whN2dTdz) ],'-dpng')

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