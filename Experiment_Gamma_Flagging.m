%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% See how estimates of gamma are affected by differnet flagging of N^2,
% epsilon etc..
%
%
%-----------------
% 10/5/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/TestChiMethod/

% load EQ14 data
dir_base='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed';

clear cham
load( fullfile( dir_base, '/Cstar=0_01366/sum/eq14_sum.mat') )
clear idb

chi=cham.CHI(:);
eps=cham.EPSILON(:);
n2=cham.N2(:);
dtdz=cham.DTDZ(:);

% no flaggin/threholding
%gam1=est_gam(cham);
%clear cham
gam1=n2 .* chi ./2 ./ eps ./ (dtdz.^2);

%ib=find(gam1>10)

%%


%% Plot full distributions from chameleon

figure(1);clf
agutwocolumn(1)
wysiwyg
set(gcf,'defaultaxesfontsize',15)

subplot(221)
h1=histogram(log10(eps),'edgecolor','none','Normalization','pdf')
grid on
xlabel('\epsilon')
xlim([-12 -3])

subplot(222)
h1=histogram(log10(chi),'edgecolor','none','Normalization','pdf')
grid on
xlabel('\chi')
xlim([-13 -2])

subplot(223)
h1=histogram(real(log10(n2)),'edgecolor','none','Normalization','pdf')
grid on
xlabel('N2')
xlim([-7 -2])

subplot(224)
h1=histogram(real(log10(dtdz)),'edgecolor','none','Normalization','pdf')
hold on
grid on
xlabel('dTdz')
xlim([-5 -0])

%% apply some thresholds based on these distributions (assume they should be normal?)

ib=find(log10(eps)<-8);
eps2=eps;
eps2(ib)=nan;

clear ib
ib=find(log10(n2)<-5);
n22=n2;
n22(ib)=nan;

clear ib
ib=find(log10(dtdz)<-2.5);
dtdz2=dtdz;
dtdz2(ib)=nan;

chi2=chi;

gam2=n22 .* chi2 ./2 ./ eps2 ./ (dtdz2.^2);


%%
figure(1);clf
agutwocolumn(1)
wysiwyg
set(gcf,'defaultaxesfontsize',15)

subplot(221)
h1=histogram(log10(eps2),'edgecolor','none','Normalization','pdf')
grid on
xlabel('\epsilon')
xlim([-12 -3])

subplot(222)
h1=histogram(log10(chi2),'edgecolor','none','Normalization','pdf')
grid on
xlabel('\chi')
xlim([-13 -2])

subplot(223)
h1=histogram(real(log10(n22)),'edgecolor','none','Normalization','pdf')
grid on
xlabel('N2')
xlim([-7 -2])

subplot(224)
h1=histogram(real(log10(dtdz2)),'edgecolor','none','Normalization','pdf')
hold on
grid on
xlabel('dTdz')
xlim([-5 -0])

%%


ig1=find(gam1<2);
ig2=find(gam2<2);

figure(2);clf
histogram(gam1(ig1),'edgecolor','none','Normalization','pdf')
hold on
histogram(gam2(ig2),'edgecolor','none','Normalization','pdf')
xlim([0 1])

%%
nanmean(gam1(:))
nanmean(gam2(:))

figure(3);clf

ax1=subplot(121);
boxplot(log10(gam1(:)))
ylabel('log_{10}\Gamma')
grid on

ax2=subplot(122);
boxplot(log10(gam2(:)))
ylabel('log_{10}\Gamma')
grid on

linkaxes([ax1 ax2])
%%
figure(1);clf

subplot(221)
h1=histogram(log10(cham.EPSILON(:)),'edgecolor','none','Normalization','pdf')
hold on
h2=histogram(log10(cham.EPSILON(ib)),'Normalization','pdf')
grid on
legend([h1 h2],'full','big \Gamma')
xlabel('\epsilon')

subplot(222)
h1=histogram(log10(cham.CHI(:)),'edgecolor','none','Normalization','pdf')
hold on
h2=histogram(log10(cham.CHI(ib)),'Normalization','pdf')
grid on
xlabel('\chi')

subplot(223)
h1=histogram(real(log10(cham.N2(:))),'edgecolor','none','Normalization','pdf')
hold on
h2=histogram(log10(cham.N2(ib)),'Normalization','pdf')
grid on
xlabel('N2')

subplot(224)
h1=histogram(real(log10(cham.DTDZ(:))),'edgecolor','none','Normalization','pdf')
hold on
h2=histogram(real(log10(cham.DTDZ(ib))),'Normalization','pdf')
grid on
xlabel('dTdz')
legend([h1 h2],'full','big \Gamma')

% From above plot, very large gammas seem to be associated mostly with
% small dTdz, and larger epsilon

%% remove values where epsilon is <noise floor
clear ib
ib=find(log10(cham.EPSILON)<-8);
gam2=gam1;
gam2(ib)=nan;

nanmax(gam1(:))
nanmax(gam2(:))

nanmean(log10(gam1(:)))
nanmean(log10(gam2(:)))

%%

clear ib
ib=find(gam1>10)

figure(1);clf
h1=histogram(gam1(:),40,'edgecolor','none','Normalization','pdf')
hold on
h2=histogram(gam2(:),'edgecolor','none','Normalization','pdf')
%xlim([0 1])

%%

clear ib gam3
ib=find(log10(cham.DTDZ)<-3);

gam3=gam1;
gam3(ib)=nan;

nanmax(gam3(:))
nanmean(gam3(:))

figure(1);clf
histogram(gam3(:))
freqline(nanmean(gam3(:)))
xlim([0 1])

%%

clear ib gam4
ib=find(log10(cham.DTDZ)<-3 | log10(cham.N2)>3 );

gam4=gam1;
gam4(ib)=nan;

nanmax(gam4(:))
nanmean(gam4(:))
nanmedian(gam4(:))
%%

ig=find(gam4<2);

figure(1);clf
histogram(gam4(ig),500,'edgecolor','none','Normalization','pdf');
freqline(nanmean(gam4(ig)));
%freqline(nanmedian(gam4(:)));
xlim([0 1])

clear ib 
ib=find(gam4>10);
length(ib)
%%