%~~~~~~~~~~~~~~~~~~~~~~~~~
%
% LookAtFlux91Data_Bill.m
%
% Look at flux 91 data that *Bill* shared. Will compute gamma etc. similar to
% EQ14 and compare results..
%
% I downloaded the data to my laptop at /Chipod/Flux91/
%
% * I think this data is only for patches?
%
%-------------
% 10/18/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

addpath /Users/Andy/Cruises_Research/GenMatlabFunctions/

load('/Users/Andy/Cruises_Research/ChiPod/Flux91/events_FLX91.mat')

% get variables
n2=A.N2;
dtdz=A.tgrad;
chi=A.chi;
eps=A.eps;

% compute gamma from this data
gam =  n2 .* chi ./2 ./ eps ./ (dtdz.^2);


%% plot distrbutions of variables

figure(1);clf
agutwocolumn(1)
wysiwyg

h=Plot4hist(n2,dtdz,chi,eps)

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures'
print( fullfile( figdir,'flx91_patches_Bill_hist'), '-dpng')

%% Plot histogram of gamma computed from these values

figure(1);clf
agutwocolumn(0.5)
wysiwyg

histogram(gam(:),30);
freqline(nanmedian(gam));
title(['median=' num2str(roundx(nanmedian(gam),2))])
xlabel('\Gamma')
grid on

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures'
print( fullfile( figdir,'flx91_patches_Bill_hist_gam'), '-dpng')

%% Scatter Plot of gamma vs each variable
figure(1);clf
agutwocolumn(1)
wysiwyg

%yl=[-6 2];
yl=[-1.5 1];
%yl=[-3 0]

ax1=subplot(221);
%histogram2(real(log10(n2)),log10(gam),200,'DisplayStyle','tile')
scatter(real(log10(n2)),log10(gam),'filled','MarkerFaceAlpha',0.15)
xlabel('log_{10}N^2','fontsize',16)
ylabel('log_{10}\Gamma','fontsize',16)
ylim(yl)
grid on
hold on

xvec=linspace(-6,-2,100);
plot(xvec,xvec+3,'k--')


ax2=subplot(222);
%histogram2(real(log10(dtdz)),log10(gam),200,'DisplayStyle','tile')
scatter(real(log10(dtdz)),log10(gam),'filled','MarkerFaceAlpha',0.15)
xlabel('log_{10}dT/dz','fontsize',16)
ylabel('log_{10}\Gamma','fontsize',16)
ylim(yl)
grid on
hold on
xvec=linspace(-4,0,100);
plot(xvec,xvec+1,'k--')


ax3=subplot(223);
%histogram2(real(log10(chi)),log10(gam),200,'DisplayStyle','tile')
scatter(real(log10(chi)),log10(gam),'filled','MarkerFaceAlpha',0.15)
xlabel('log_{10}\chi','fontsize',16)
ylabel('log_{10}\Gamma','fontsize',16)
ylim(yl)
xlim([-12 -4])
grid on
hold on
xvec=linspace(-12,-4,100);
plot(xvec,xvec+7,'k--')


ax4=subplot(224);
%h=histogram2(real(log10(eps)),log10(gam),200,'DisplayStyle','tile')
scatter(real(log10(eps)),log10(gam),'filled','MarkerFaceAlpha',0.15);
xlabel('log_{10}\epsilon','fontsize',16)
ylabel('log_{10}\Gamma','fontsize',16)
ylim(yl)
grid on
xlim([-10 -5])
hold on
xvec=linspace(-10,-5,100);
plot(xvec,xvec+7,'k--')

linkaxes([ax1 ax2 ax3 ax4],'y')

%%