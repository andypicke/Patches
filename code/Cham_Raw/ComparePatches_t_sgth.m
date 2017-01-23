%%
%
% ComparePatches_t_sgth.m
%
% See if patches from t/density line up
%
%
%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/EQ14_raw_patches_minOT_25_usetemp_1.mat')
P1=new_patch_data; clear new_patch_data

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/EQ14_raw_patches_minOT_25_usetemp_0.mat')
P2=new_patch_data; clear new_patch_data


%%

cnum=55;


datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal'

load(fullfile(datdir,['eq14_' sprintf('%04d',cnum) '.mat']))


i1=find(P1(:,1)==cnum);
i2=find(P2(:,1)==cnum);

cal2.sgth=sw_pden(cal2.SAL,cal2.T1,cal2.P,0);

figure(1);clf

ax1=subplot(121);
plot(cal2.T1,cal2.P,'color',0.5*[1 1 1],'linewidth',2)
hold on
axis ij
axis tight
grid on

for ip=1:length(i1)
    iz=isin(cal2.P,[P1(i1(ip),2) P1(i1(ip),3)]);
    plot(cal2.T1(iz),cal2.P(iz),'k','linewidth',3)
end

ax2=subplot(122);
plot(cal2.sgth,cal2.P,'color',0.5*[1 1 1],'linewidth',2)
hold on
axis ij
axis tight
grid on

for ip=1:length(i2)
    iz=isin(cal2.P,[P2(i2(ip),2) P2(i2(ip),3)]);
    plot(cal2.sgth(iz),cal2.P(iz),'k','linewidth',3)
end


linkaxes([ax1 ax2],'y')
