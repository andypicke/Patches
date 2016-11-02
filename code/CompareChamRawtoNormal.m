%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CompareChamRawtoNormal.m
%
% Compare chameleon data computed just for patches to the normal 1m data to
% make sure they mostly agree.
%
% 11/01/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

for cnum=1:50:3100
    
    try
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/eq14_' sprintf('%04d',cnum) '.mat'])
        avg1=avg;clear avg
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/mat/eq14_' sprintf('%04d',cnum) '.mat'])
        avg2=avg;clear avg
        
        %
        figure(1);clf
        semilogx(avg2.EPSILON,avg2.P)
        axis ij
        hold on
        semilogx(avg1.EPSILON,avg1.P,'o')
        title(cnum)
    end
    
    pause(0.5)
    
end
%%