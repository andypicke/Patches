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
%% Identify patches in 1-m binned TIWE profiles

clear ; close all

patches=[];

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/

% load combined TIWE data
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
        
end % ip

warning on

% save data
save( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/','tiwe_patches.mat'), 'patches')

%%

clear ; close all

load(fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/','tiwe_patches.mat'))

figure(1);clf
histogram(patches(:,4))
xlabel('patch size')
ylabel('count')
title('tiwe')
xlim([0 150])
grid on

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print( fullfile( figdir,'tiwe_patchsize_hist'), '-dpng')


%% Plot depth of patches; most are near surface or bottom

figure(2);clf
histogram(patches(:,2),'Normalization','pdf')
hold on
histogram(patches(:,3),'Normalization','pdf')

%% don't use patches near the surface/mixed layer or bottom of profiles

ig=find(patches(:,2)>60 & patches(:,3)<150 );
length(ig)
length(patches)
disp([num2str(length(ig)/length(patches)*100) ' percent of total patches were not near surface or bottom'])

patches_trim=patches(ig,:);

save( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/','tiwe_patches_trim.mat'), 'patches_trim')

% plot patch size distribution again for the trimmed set
figure(2);clf
histogram(patches_trim(:,4))
figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print( fullfile( figdir,'tiwe_patchsize_mid_hist'), '-dpng')

%% Plot patch locataions on pcolor plot of data

% load combined TIWE data
load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches','tiwe_comb_AP.mat') )
cham=tiwe;

patch=patches;

close all
figure(2);clf
ezpc(cham.ip,cham.P,log10(cham.EPSILON));caxis([-8 -5])
cb=colorbar;
cmap=colormap;cmap=[0.75*[1 1 1] ; cmap  ];
colormap(cmap)
cb.Label.String='log_{10}\epsilon'
hold on
plot(patch(:,1),patch(:,2),'g.','markersize',10)
plot(patch(:,1),patch(:,3),'r.','markersize',10)
xlim([500 3000])
ylim([0 200])
xlabel('castnumber')
ylabel('P')

title('TIWE Chameleon 1m binned profiles')


%%

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print( fullfile( figdir,'tiwe_pcolor_wpatches'), '-dpng')


%% Plot trimmed patch locataions on pcolor plot of data

patch=patches_trim;

close all
figure(2);clf
ezpc(cham.ip,cham.P,log10(cham.EPSILON));caxis([-8 -5])
cb=colorbar;
cmap=colormap;cmap=[0.75*[1 1 1] ; cmap  ];
colormap(cmap)
cb.Label.String='log_{10}\epsilon'
hold on
plot(patch(:,1),patch(:,2),'g.','markersize',10)
plot(patch(:,1),patch(:,3),'r.','markersize',10)
xlim([500 3000])
ylim([0 200])
xlabel('castnumber')
ylabel('P')

hline(150,'k')
hline(60,'k')
title('TIWE Chameleon 1m binned profiles')


%%

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print( fullfile( figdir,'tiwe_pcolor_wpatches_mid'), '-dpng')


%% Now get all data from just within patches

n2_patches  =[];
dtdz_patches=[];
chi_patches =[];
eps_patches =[];

the_patches=patches; % choose which set of patches to use

for ipatch=1:length(the_patches)

    clear prange iz ip
    prange=[the_patches(ipatch,2) the_patches(ipatch,3)]; % prange of patch
    iz=isin(tiwe.P,prange+[-1 1]); % indices in patch
    ip=the_patches(ipatch,1); % profile #
    
    n2_patches  =[n2_patches   ; tiwe.N2(iz,ip)     ];
    dtdz_patches=[dtdz_patches ; tiwe.DTDZ(iz,ip)   ];
    chi_patches =[chi_patches  ; tiwe.CHI(iz,ip)    ];
    eps_patches =[eps_patches  ; tiwe.EPSILON(iz,ip)];
    
end % ipatch

%% save this data

%save( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/','tiwe_patches_trim_data.mat'), 'n2_patches','dtdz_patches','chi_patches','eps_patches')
save( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/','tiwe_patches_data.mat'), 'n2_patches','dtdz_patches','chi_patches','eps_patches')

%%

clear ; close all

% load data from patches
load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/','tiwe_patches_data.mat') )

% load combined TIWE data
load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches','tiwe_comb_AP.mat') )
cham=tiwe;


%% Plot distributions of these data

h=Plot4hist(n2_patches,dtdz_patches,chi_patches,eps_patches);

%% Compare distributions from all points in patch with patch-average values

nbins=50
m=2;
n=2;

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(m,n,1)
h1=histogram(real(log10(cham.N2(:))),nbins,'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(n2_patches(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
grid on
xlim([-7 -2])
legend([h1 h2],'All data','patches','location','best')
xlabel('log_{10}N^2','fontsize',15)
ylabel('pdf')
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')

subplot(m,n,2)
h1=histogram(real(log10(cham.DTDZ(:))),nbins,'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(dtdz_patches(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
grid on
xlim([-5 0])
legend([h1 h2],'All data','patches','location','best')
xlabel('log_{10}dT/dz','fontsize',15)
ylabel('pdf')
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')


subplot(m,n,3)
h1=histogram(real(log10(cham.CHI(:))),nbins,'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(chi_patches(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
grid on
xlim([-13 -4])
legend([h1 h2],'All data','patches','location','best')
xlabel('log_{10}\chi','fontsize',15)
ylabel('pdf')
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')


subplot(m,n,4)
h1=histogram(real(log10(cham.EPSILON(:))),nbins,'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(eps_patches(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
grid on
xlim([-11 -4])
legend([h1 h2],'All data','patches','location','best')
xlabel('log_{10}\epsilon','fontsize',15)
ylabel('pdf')
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')

%%
figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures'
print(fullfile(figdir,'tiwe_hists_allVspatch'),'-dpng')


%% Compare my patch data to Bill's patch data

clear ; close all

% load my patch data
load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/','tiwe_patches_trim_data.mat') )
%load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/','tiwe_patches_data.mat') )

% load Bill's patch data
load('/Users/Andy/Cruises_Research/ChiPod/TIWE/events_TIWE.mat')

% only use values above epsilon noise floor?
%ig=find(log10(A.eps)>-8);

% get variables
n2=A.N2;
dtdz=A.tgrad;
chi=A.chi;
eps=A.eps;

[h h1 h2]=Plot4hist_2(n2,dtdz,chi,eps,n2_patches,dtdz_patches,chi_patches,eps_patches)
legend([h1 h2],'Bill patches','AP patches','location','best')

%% Compare gamma

% compute gamma from this data
gam_bill =  n2 .* chi ./2 ./ eps ./ (dtdz.^2);
gam_AP= n2_patches .* chi_patches ./2 ./eps_patches ./ (dtdz_patches.^2);

ig1=find(gam_bill<5);
ig2=find(gam_AP<5);

figure(1);clf
agutwocolumn(0.5)
wysiwyg

h1=histogram(gam_bill(ig1),'Normalization','pdf');
hold on
h2=histogram(gam_AP(ig2),'Normalization','pdf');
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')
grid on
text(2,2.5,['median=' num2str(nanmedian(gam_bill))],'color','b','fontsize',16)
text(2,2.2,['median=' num2str(nanmedian(gam_AP))],'color','r','fontsize',16)
legend([h1 h2],'Bill patches','AP patches 1m')
xlabel('\Gamma','fontsize',16)
ylabel('pdf','fontsize',16)
title('\Gamma from tiwe patches')

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures'
print(fullfile(figdir,'tiwe_gam_billvsAP'),'-dpng')

%%