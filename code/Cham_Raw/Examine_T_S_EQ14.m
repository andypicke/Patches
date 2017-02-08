%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Examine_T_S_EQ14.m
%
%
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/

cnum=2024
cham_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal';
load( fullfile( cham_dir, ['eq14_' sprintf('%04d',cnum) '.mat'] ) )
%            cnum_loaded=cnum;
cal2.ptmp=sw_ptmp(cal2.SAL,cal2.T1,cal2.P,0);
cal2.sgth=sw_pden(cal2.SAL,cal2.T1,cal2.P,0);

figure(1);clf
loglog(cal2.ptmp,cal2.sgth,'.')
grid on
% axis equal
P=polyfit(cal2.ptmp,cal2.sgth,1)
hold on
plot(cal2.ptmp,polyval(P,cal2.ptmp),'linewidth',2)
axis tight
xlabel('pot. temp')
ylabel('pot. dens')
%%


sgth=sw_pden(cal.SAL,cal.T1,cal.P,0);

figure(1);clf
plot(cal.SIGMA+1000,cal.P)
hold on
plot(sgth,cal.P,'--')
axis ij

%%

ptmp=sw_ptmp(cal.SAL,cal.T1,cal.P,0);

figure(1);clf
plot(cal.theta,cal.P)
hold on
plot(ptmp,cal.P,'--')
axis ij

%%