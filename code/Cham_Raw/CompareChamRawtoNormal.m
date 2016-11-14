%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareChamRawtoNormal.m
%
%
% Compare chameleon data computed just for patches to the normal 1m averaged
% data to see if they mostly agree.
%
% Patch data are computed in run_eq14_for_PATCHES
%
%--------------
% 11/01/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

minOTsize=1


for cnum=1:50:3100
    
    try
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/minOT_' num2str(10*minOTsize) '/eq14_' sprintf('%04d',cnum) '.mat'])
        avg1=avg;clear avg
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/mat/eq14_' sprintf('%04d',cnum) '.mat'])
        avg2=avg;clear avg
        
        %
        figure(1);clf
        agutwocolumn(1)
        wysiwyg
        
        subplot(221)
        semilogx(avg2.N2,avg2.P)
        axis ij
        hold on
        semilogx(avg1.N2,avg1.P,'kd','linewidth',3)
        semilogx(avg1.n2_OT,avg1.P,'rd','linewidth',3)
        grid on
        ylabel('P [db]')
        axis tight
        xlabel('log_{10}N^2','fontsize',16)
        title(['profile ' num2str(cnum)])
        
        subplot(222)
        semilogx(avg2.DTDZ,avg2.P)
        axis ij
        hold on
        semilogx(avg1.dTdz,avg1.P,'kd','linewidth',3)
        grid on
        xlabel('log_{10}dT/dz','fontsize',16)
        ylabel('P [db]')
        axis tight
        
        subplot(223)
        semilogx(avg2.CHI,avg2.P)
        axis ij
        hold on
        semilogx(avg1.CHI,avg1.P,'kd','linewidth',3)
        grid on
        ylabel('P [db]')
        axis tight
        xlabel('log_{10}\chi','fontsize',16)
        
        subplot(224)
        semilogx(avg2.EPSILON,avg2.P)
        axis ij
        hold on
        semilogx(avg1.EPSILON,avg1.P,'kd','linewidth',3)
        grid on
        ylabel('P [db]')
        axis tight
        xlabel('log_{10}\epsilon','fontsize',16)
        
        
    end
    
    pause(0.5)
    
end
%%