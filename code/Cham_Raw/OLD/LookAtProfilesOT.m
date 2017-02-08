%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% LookAtProfilesOT.m
%
% Look at individual profiles and overturns in chameleon eq14 data
%
% Trying to get a better idea of what hte overturns look like, if they are
% correct, and also looking at different methods of computing N2 and dT/dz.
%
%
%---------------
% 11/16/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

cnum=1050; % choose a cast # to look at

% load overturn data

% load chameleon cast
cham_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal'
load( fullfile( cham_dir, ['eq14_' sprintf('%04d',cnum) '.mat'] ) )

figure(1);clf

plot(cal2.T1,cal2.P)
axis ij
grid on

ptmp=sw_ptmp(cal2.SAL,cal2.T1,cal2.P,100);
hold on
plot(ptmp,cal2.P)

figure(2);clf
plot(cal2.SAL,cal2.P)
axis ij

%%
clear cal
cal=cal2;

clear s t p lat
s=cal.SAL(1:end-1); % (end-1) b/c last 2 values are same;
t=cal.T1(1:end-1);
p=cal.P(1:end-1);

clear idot lat1 lat2
idot=strfind(head.lat.start,'.');
lat1=str2num(head.lat.start(1:idot-3));
lat2=str2num(head.lat.start(idot-2:end))/60;
lat=nanmean([lat1 lat2]);

% find overturns
%        clear pstarts pstops
%       [pstarts pstops]=IdentifyPatches(s,t,p,lat);

clear Params
Params.lat=lat;
Params.plotit=1;
Params.sigma=1e-4;
Params.runlmin=0;
Params.minotsize=1;
Params.usetemp=1;
addpath /Users/Andy/Standard-Mixing-Routines/ThorpeScales/
clear OT
OT=compute_overturns_discrete_AP(p,t,s,Params);

pstarts=OT.pstarts_each;
pstops=OT.pstops_each;


%% compute N^2 and dT/dz for patch
% 1) fit a line over patch
% 2) use max-min value of rho or dT/dz
% 3) use Bill's 'bulk' method


for iot=1:length(OT.pstarts_each)
    figure(1);clf
    plot(t,p,'o-')
    ylim( [ OT.pstarts_each(iot) OT.pstops_each(iot) ] + [-1 1])
    hold on
    hline(OT.pstarts_each(iot),'--')
    hline(OT.pstops_each(iot),'--')
    axis ij
    pause
    
end

%% Look at a single patch

iot=8

figure(1);clf
plot(t,p,'k-','linewidth',2)
ylim( [ OT.pstarts_each(iot) OT.pstops_each(iot) ] + 1*[-1 1])
hold on
hline(OT.pstarts_each(iot),'k--')
hline(OT.pstops_each(iot),'k--')
axis ij
grid on
xlabel('T')
ylabel('P [db]')

iz=isin(p,[ OT.pstarts_each(iot) OT.pstops_each(iot) ]);
t_ot=t(iz);
s_ot=s(iz);
p_ot=p(iz);
hold on
plot(t_ot,p_ot,'r-','linewidth',2)

[t_sort I]=sort(t_ot,1,'descend');
h1=plot(t_sort,p_ot,'b-','linewidth',2)


dT=nanmax(t_ot)-nanmin(t_ot);
dz=nanmax(p_ot)-nanmin(p_ot);
dTdz=-dT/dz


b=p_ot(end) - (1/dTdz)*t_ot(end);
tvec=linspace(nanmin(t_ot),nanmax(t_ot),100);
pvec=linspace(nanmin(p_ot),nanmax(p_ot),100);

yvec=dTdz*pvec - dTdz*b;
h2=plot(yvec,pvec,'m','linewidth',2)

% fit a line
P=polyfit(p_ot,t_sort,1);
h3=plot(polyval(P,pvec),pvec,'g','linewidth',2)

legend([h1 h2 h3],'T sort','dT/dz','linefit')
%%

dTdz
P(1)

%% Now do similar for density

sgth=sw_pden(s,t,p,0);
sgth_ot=sw_pden(s_ot,t_ot,p_ot,0);

figure(2);clf
plot(sgth,p,'k','linewidth',2)
hold on
plot(sgth_ot,p_ot,'r','linewidth',2)
ylim( [ OT.pstarts_each(iot) OT.pstops_each(iot) ] + [-1 1])
hold on
hline(OT.pstarts_each(iot),'k--')
hline(OT.pstops_each(iot),'k--')
axis ij
grid on
xlabel('pot. dens')
ylabel('pres')

[sgth_sort I]=sort(sgth_ot,1,'ascend');
plot(sgth_sort,p_ot,'b','linewidth',2)

% compute N^2
n2=sw_bfrq(s_ot,t_ot,p_ot,1);
nanmean(n2)

% fit a line to sgth
P1=polyfit(p_ot,sgth_sort,1);
plot(polyval(P1,pvec),pvec,'g')

% calculate N^2 from this fit and compare to other N2 estimates
drhodz=P1(1);
n2_2=9.81/nanmean(sgth_ot)*drhodz;

% try the range/dz method
drho=nanmax(sgth_ot)-nanmin(sgth_ot)
dz=nanmax(p_ot)-nanmin(p_ot)
n2_3=9.81/nanmean(sgth_ot)*drho/dz

%%
nanmean(n2)
n2_2
n2_3

%%
% bill's method ('bulk')




%%