%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareProfiles_bin_patch_cham.m
%
% Compare profiles of chi/epsilon for Chameleon (shear probe), chi-pod
% method applied to binned data, and chi-pod method applied to patches
%
%-------------
% 12/15/16 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

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
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
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
%%