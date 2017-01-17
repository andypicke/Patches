%~~~~~~~~~~~~~~~~~~~~~~~
%
% LookAtN2Profiles.m
%
%
% 1/9/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

cnum=1355

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/


load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal/EQ14_' sprintf('%04d',cnum) '.mat'])
clear cal
cal=cal2;
clear cal2

% smooth salinity a bit because it is noisy
cal.SALsm=smooth(cal.SAL,30);
figure(1);clf
plot(cal.SAL,cal.P,'.')
hold on
plot(cal.SALsm,cal.P)
axis ij
grid on

%%
% calc 1m binned T,S
% Average temp and sal in 1m bins like we normally do for CTD data
clear zmin dz zmax tbin zbin sbin
zmin=nanmin(cal.P);
dz=1;
zmax=nanmax(cal.P);
minobs=2;
[tbin zbin Nobs] = binprofile(cal.T1 ,cal.P, zmin, dz, zmax,minobs);
[sbin zbin Nobs] = binprofile(cal.SAL,cal.P, zmin, dz, zmax,minobs);
clear zmin dz zmax minobs
cal.T1bin=tbin;
cal.SALbin=sbin;
cal.Pbin=zbin;

figure(1);clf
ax1=subplot(131);
h1=plot(cal.T1,cal.P);
hold on
h2=plot(cal.T1bin,cal.Pbin);
axis ij
xlabel('T')
grid on
ylabel('P [db]')
axis tight
legend([h1 h2],'raw','bin','location','best')

ax2=subplot(132);
plot(cal.SALsm,cal.P)
hold on
plot(cal.SALbin,cal.Pbin)
axis ij
axis tight
xlabel('S')
grid on

sgth=sw_pden(cal.SALsm,cal.T1,cal.P,0);
sgthbin=sw_pden(cal.SALbin,cal.T1bin,cal.Pbin,0);

ax3=subplot(133);
plot(sgth,cal.P)
hold on
plot(sgthbin,cal.Pbin)
axis ij
xlabel('sgth')
grid on
axis tight

linkaxes([ax1 ax2 ax3],'y')

%%
figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures'
print(fullfile(figdir,['cnum_' num2str(cnum) '_T_S_sgth']),'-dpng')

%%

% calc N^2 from raw (high-res) profile
n2=sw_bfrq(cal.SAL,cal.T1,cal.P,0.3);
n2=[n2 ; nan];

% calc N^2 from raw (high-res) profile
dTdz=diffs(cal.T1) ./ diffs(cal.P);

% calc N^2 from 1m binned 
n2bin=sw_bfrq(cal.SALbin,cal.T1bin,cal.Pbin,0.3);
n2bin=[n2bin ; nan];

% calc dTdz from 1m binned 
dTdzbin=diffs(cal.T1bin) ./ diffs(cal.Pbin);


figure(2);clf
agutwocolumn(0.7)
wysiwyg

ax1=subplot(131);
h1=plot(real(log10(n2)),cal.P,'.');
hold on
h2=plot(real(log10(n2bin)),cal.Pbin,'-','linewidth',2);
axis ij
xlabel('log_{10}[N^2]')
ylabel('P [db]')
grid on
xlim([-6 -1.5])
legend([h1 h2],'raw','1mbin','location','best')

ax2=subplot(132);
h1=plot(real(log10(dTdz)),cal.P,'.');
hold on
h2=plot(real(log10(dTdzbin)),cal.Pbin,'-','linewidth',2);
axis ij
xlabel('log_{10}[dT/dz]')
%ylabel('P [db]')
grid on
xlim([-4 0.5])
legend([h1 h2],'raw','1mbin','location','best')
ytloff

% compute ratio N^2/dTdz^2 for both
rat= n2 ./ (dTdz.^2);
ratbin= n2bin ./ (dTdzbin.^2);

ax3=subplot(133);
plot(log10(rat),cal.P,'.')
hold on
plot(log10(ratbin),cal.Pbin,'-','linewidth',2)
axis ij
grid on
axis tight
xlabel('log_{10}[N^2/dTdz^2]')

linkaxes([ax1 ax2 ax3],'y')

%%
figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures'
print(fullfile(figdir,['cnum_' num2str(cnum) '_N2_dTdz_ratio']),'-dpng')

%% try sorting profiles by sgth first

[Y,I] = sort(sgth,1,'ascend') ;

figure(1) ; clf

ax1=subplot(121);
plot(sgth,cal.P)
hold on
plot(sgth(I),cal.P)
axis ij

T1sort=cal.T1(I);

ax2=subplot(122);
plot(cal.T1,cal.P)
hold on
plot(T1sort,cal.P)
axis ij
grid on

linkaxes([ax1 ax2],'y')

% calc N^2 from raw (high-res) profile
n2sort=sw_bfrq(cal.SALsm,T1sort,cal.P,0.3);
n2sort=[n2sort ; nan];

% calc N^2 from raw (high-res) profile
dTdzsort=diffs(T1sort) ./ diffs(cal.P);

% calc N^2 from 1m binned 
n2bin=sw_bfrq(cal.SALbin,cal.T1bin,cal.Pbin,0.3);
n2bin=[n2bin ; nan];

% calc dTdz from 1m binned 
dTdzbin=diffs(cal.T1bin) ./ diffs(cal.Pbin);


figure(2);clf
agutwocolumn(0.7)
wysiwyg

ax1=subplot(131);
h1=plot(real(log10(n2sort)),cal.P,'.');
hold on
h2=plot(real(log10(n2bin)),cal.Pbin,'-','linewidth',2);
axis ij
xlabel('log_{10}[N^2]')
ylabel('P [db]')
grid on
xlim([-6 -1.5])
legend([h1 h2],'raw','1mbin','location','best')

ax2=subplot(132);
h1=plot(real(log10(dTdzsort)),cal.P,'.');
hold on
h2=plot(real(log10(dTdzbin)),cal.Pbin,'-','linewidth',2);
axis ij
xlabel('log_{10}[dT/dz]')
grid on
xlim([-4 0.5])
legend([h1 h2],'raw','1mbin','location','best')
ytloff

% compute ratio N^2/dTdz^2 for both
rat= n2sort ./ (dTdzsort.^2);
ratbin= n2bin ./ (dTdzbin.^2);

ax3=subplot(133);
plot(log10(rat),cal.P,'.')
hold on
plot(log10(ratbin),cal.Pbin,'-','linewidth',2)
axis ij
grid on
%axis tight
xlim([-4 2])
linkaxes([ax1 ax2 ax3],'y')

%%

figure(1);clf
plot(log10(n2),cal.P,'.')
hold on
plot(log10(n2sort),cal.P,'.')

figure(2);clf
loglog(n2(:),n2sort(:),'.')

to_plot=real(log10(n2sort./n2));
figure(3);clf
histogram(to_plot)
freqline(nanmean(to_plot),'r--')
%%

figure(3);clf
h1=histogram(real(log10(rat(:))),30,'Normalization','pdf');
hold on
h2=histogram(real(log10(ratbin(:))),h1.BinEdges,'Normalization','pdf');
grid on
legend([h1 h2],'raw','bin')
xlim([-4 1])
xlabel('log_{10}[N^2/dTdz^2]')
ylabel('pdf')
title(['Profile ' num2str(cnum) ])

%% what if we then interpolate each to same 1m points?

[C,IA,IC] = unique(cal.P);
n2i=interp1(cal.P(IA),n2(IA),cal.Pbin);

figure(3);clf
%loglog(n2i(:),n2bin(:),'.')
plot(log10(n2i),cal.Pbin,'.-')
hold on
plot(log10(n2bin),cal.Pbin,'.-')
axis i
%%