%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Find_patches_TIWE_cham.m
%
% Goal is to identify patches (via overturns codes etc.), and then compute
% gamma using values averaged over those patches.
%
% * (1) Loop through all files and find overturns, make table of cast#,
% overturn depths, patch size etc.
%
% Modified from Find_patches_EQ14_cham.m
%
%----------------------
% 10/24/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% First let's look at gamma as computed from 1m binned chameleon data

clear ; close all

% load data
% dir_base='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed';
% clear cham
% load( fullfile( dir_base, '/Cstar=0_01366/sum/eq14_sum.mat') )

load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches','tiwe_comb_AP.mat') )
cham=tiwe;

% gamma using all data points
gam_cham=cham.N2 .* cham.CHI ./2 ./ cham.EPSILON ./ (cham.DTDZ.^2);

% find points not in mixed layer or near bottom of profiles
i4=find(cham.P>40 & cham.P<180);

% i4 contains about half as many values
numel(i4)/numel(cham.P)

%%
figure(1);clf
h1=histogram(real(log10(gam_cham(:))),65,'edgecolor','none','Normalization','pdf');
hold on
h2=histogram(log10(gam_cham(i4)),h1.BinEdges,'edgecolor','none','Normalization','pdf');
xlim([-6 3])
grid on
legend([h1 h2],'all','60>P<180')
xlabel('log_{10}\Gamma')
ylabel('pdf')

freqline(nanmedian(real(log10(gam_cham(:)))),'b')
freqline(nanmedian(log10(gam_cham(i4))),'r')

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
print(fullfile(figdir,'gam_cham_hist'),'-dpng')

% but the distribution doens't change. So it doesn't seem like including
% the mixed layer dirunal cycle data is affecting our estimates.

%%
% the median gamma from cham data is about 0.02, about 10 times smaller 
% than the normally assumed value of 0.2

nanmedian(gam_cham(:))
nanmedian(gam_cham(i4))


%% Identify patches in 1-m binned EQ14 profiles

clear ; close all

patches=[];

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/
% load an EQ14 chameleon profile
%dir_base='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/mat/';

load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches','tiwe_comb_AP.mat') )
cham=tiwe;

plots=0
warning off
for ip=1:length(cham.ip)
        
        % below is modified from compute_overturns_discrete_AP.m **
        clear s t p refd lat pden ptmp
        s=cham.s(:,ip);
        t=cham.t(:,ip);
        p=cham.P(:);
        refd=0;
        lat=1.7;%str2num(head.lat.start);%1.7;%nanmean(head.lat.start)

        [pstarts pstops]=IdentifyPatches(s,t,p,lat);
                
        
        for i=1:length(pstarts)
            patches=[patches ; ip pstarts(i) pstops(i)  ( pstops(i) - pstarts(i) ) ];
        end
        
        %
    %end % try
    
end % ip

warning on
%%

figure(1);clf
histogram(patches(:,4))

%%

figure(2);clf
histogram(patches(:,2),'Normalization','pdf')
hold on
histogram(patches(:,3),'Normalization','pdf')

%%

% don't use patches near the surface/mixed layer or bottom of profiles
ig=find(patches(:,2)>60 & patches(:,3)<180 );
length(ig)
length(patches)
disp([num2str(length(ig)/length(patches)*100) ' percent of total patches were not near surface or bottom'])


patches_2=patches(ig,:);

% save new variable w/ just the data in good depth range
%patches_60_180=dat_sum(ig,:);

% this reduces the number by a lot since majority were in mixed layer where
% diurnal convection dominates

figure
histogram(patches(:,4),'Normalization','pdf')
hold on
histogram(patches_2(:,4),'Normalization','pdf')

%% Plot patch locataions on pcolor plot of data


patch=patches_2;

close all
figure(2);clf
ezpc(cham.ip,cham.P,log10(cham.EPSILON));caxis([-8 -5])
%ezpc(cham.castnumber,cham.P,cham.SIGMA);caxis([20 28])
cb=colorbar;
cmap=colormap;cmap=[0.75*[1 1 1] ; cmap  ];
colormap(cmap)
cb.Label.String='log_{10}\epsilon'
hold on
plot(patch(:,1),patch(:,2),'g.','markersize',10)
plot(patch(:,1),patch(:,3),'r.','markersize',10)
ylim([0 200])
xlim([500 3000])
xlabel('castnumber')
ylabel('P')

% Plot the really large patch locations also
hold on
i2=find(patch(:,4)>40);
plot(patch(i2,1),patch(i2,2),'kd')
plot(patch(i2,1),patch(i2,3),'kp')
hline(180,'k')
hline(60,'k')
title('TIWE Chameleon 1m binned profiles')



%%

%% save results

datadir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
%save( fullfile( datadir, 'eq14_1m_Template.mat'), 'patches_all','patches_60_180')


%% Plot histogram of patch sizes

figure(1);clf
%histogram(dat_sum(ig,4))
histogram(patches(:,4),80)
xlabel('Patch size [m]')
title('Distribution of patch sizes from EQ14 Cham profiles')
xlim([0 20])

%% Plot boxplot of patch sizes

figure(1);clf
boxplot(dat_sum(ig,4))
%boxplot(dat_sum(:,4))
ylabel('patch size [db]')
grid on
title('Distribution of patch sizes from EQ14 Cham profiles')

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
print(fullfile(figdir,'GamByPatch_patchsize_boxplot'),'-dpng')


% ** most patches are only 1-2 m

%% pcolor epsilon or chi, and plot patch locations on top

% load EQ14 data
dir_base='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed';
clear cham
load( fullfile( dir_base, '/Cstar=0_01366/sum/eq14_sum.mat') )
%
close all
figure(2);clf
ezpc(cham.castnumber,cham.P,log10(cham.EPSILON));caxis([-8 -5])
%ezpc(cham.castnumber,cham.P,cham.SIGMA);caxis([20 28])
cb=colorbar;
cmap=colormap;cmap=[0.75*[1 1 1] ; cmap  ];
colormap(cmap)
cb.Label.String='log_{10}\epsilon'
hold on
plot(dat_sum(:,1),dat_sum(:,2),'g.')
plot(dat_sum(:,1),dat_sum(:,3),'r.')
ylim([0 200])
xlim([500 3000])
xlabel('castnumber')
ylabel('P')

% Plot the really large patch locations also
hold on
i2=find(Template(:,4)>20);
plot(patches(i2,1),patches(i2,2),'kd')
plot(patches(i2,1),patches(i2,3),'kd')
hline(180,'k')
hline(60,'k')
title('EQ14 Chameleon 1m binned profiles')

% save figure
print(fullfile(figdir,'GamByPatch_eps_pcolor_wpatches'),'-dpng')

%%