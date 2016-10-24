%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute_gamma_by_patches.m
%
%
%-----------
% 10/18/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% load identified patches (see Find_patches_EQ14_cham.m)
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/eq14_1m_patches.mat')
patches=patches_60_180;

% load EQ14 data
dir_base='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed';
clear cham
load( fullfile( dir_base, '/Cstar=0_01366/sum/eq14_sum.mat') )


%% Now try computing gamma by patch. 
% 1st I will try the simpler but not really correct method of using the 1m
% binned chameleon data averaged over the patch. (really should use spectra
% from entire patch to compute chi and eps?)

% Make empty arrays for results
gam_all=[];       % gamma using all points in patches
gam_patch_all=[]; % average all points in patch, then compute gamma

% also store values of n2,dtdz,chi,eps within patches so we can compare to
% full distributions
n2_patches=[];
dtdz_patches=[];
chi_patches=[];
eps_patches=[];

dir_base='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/mat/';
% loop over all identified patches
for ip=1:length(patches)

    clear cnum avg iz n2 dtdz chi eps gam gam_patch
    
    % load appropriate chameleon profile
    cnum=patches(ip,1);
    fname=['eq14_' sprintf('%04d',cnum) '.mat'];    
    load( fullfile(dir_base, fname))

    % find data points in range of patch (+/-1 m?)
    iz=isin(avg.P,[patches(ip,2) patches(ip,3)]);
            
    % get data in patch
    n2=avg.N2(iz);
%   dtdz=avg.DTDZ(iz);
    dtdz=avg.DTDZ_RHOORDER(iz);
    chi=avg.CHI(iz);
    eps=avg.EPSILON(iz);
    
    % save data within patch
    n2_patches  =[n2_patches(:)   ; n2(:)  ];
    dtdz_patches=[dtdz_patches(:) ; dtdz(:)];
    chi_patches =[chi_patches(:)  ; chi(:) ];
    eps_patches =[eps_patches(:)  ; eps(:) ];
    
    % compute gamma for each value in patch
    gam=n2.*chi./2./eps./(dtdz.^2);
       
    % compute gamma using average of values in patch
    gam_patch=nanmean(n2)*nanmean(chi)/2/nanmean(eps)/(nanmean(dtdz)^2);
   
    gam_all       = [gam_all(:)       ;gam(:)    ];
    gam_patch_all = [gam_patch_all(:) ; gam_patch];
    
end

% save the data
save(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/','eq14_INpatch.mat'),...
    'n2_patches','dtdz_patches','chi_patches','eps_patches')

%% Compare distributions from all points in patch with patch-average values

nbins=50

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(411)
h1=histogram(real(log10(cham.N2(:))),nbins,'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(n2_patches(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
grid on
xlim([-7 -2])
legend([h1 h2],'All data','patches','location','best')
xlabel('log_{10}N^2','fontsize',15)
ylabel('pdf')

subplot(412)
h1=histogram(real(log10(cham.DTDZ(:))),nbins,'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(dtdz_patches(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
grid on
xlim([-5 0])
legend([h1 h2],'All data','patches','location','best')
xlabel('log_{10}dT/dz','fontsize',15)
ylabel('pdf')

subplot(413)
h1=histogram(real(log10(cham.CHI(:))),nbins,'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(chi_patches(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
grid on
xlim([-13 -4])
legend([h1 h2],'All data','patches','location','best')
xlabel('log_{10}\chi','fontsize',15)
ylabel('pdf')

subplot(414)
h1=histogram(real(log10(cham.EPSILON(:))),nbins,'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(eps_patches(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
grid on
xlim([-11 -4])
legend([h1 h2],'All data','patches','location','best')
xlabel('log_{10}\epsilon','fontsize',15)
ylabel('pdf')

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
print(fullfile(figdir,'hists_allVspatch'),'-dpng')


%% Also compute gamma using ALL cham values (versus just those in patches)

%gam_cham=cham.N2 .* cham.CHI ./2 ./ cham.EPSILON ./ (cham.DTDZ.^2);
gam_cham=cham.N2 .* cham.CHI ./2 ./ cham.EPSILON ./ (cham.DTDZ_RHOORDER.^2);
%i3=find(gam_cham<10);
ig=find(cham.P>60 & cham.P<180);
gam_cham=gam_cham(ig);
%%
figure(1);clf
agutwocolumn(1);wysiwyg

ax1=subplot(311);
h1=histogram(log10(gam_all),50,'Normalization','pdf');
hold on
xlim([-5 4])
grid on
title(['all points in patches, median \Gamma=' num2str(roundx(nanmedian(gam_all),3))])
freqline(log10(0.2),'k--')
text(log10(0.2),0.5,'0.2')
freqline(nanmean(log10(gam_all)),'b')

ax2=subplot(312);
h2=histogram(log10(gam_patch_all),h1.BinEdges,'Normalization','pdf')
hold on
grid on
title(['patch average values, median \Gamma=' num2str(roundx(nanmedian(gam_patch_all),3))])
freqline(log10(0.2),'k--')
text(log10(0.2),0.5,'0.2')
freqline(nanmean(log10(gam_patch_all)),'b')


ax3=subplot(313);
h3=histogram(log10(gam_cham(:)),h1.BinEdges,'Normalization','pdf')
%legend([h1 h2 h3],'patch_all','patch_avg','all')
xlabel('log_{10}\Gamma')
grid on
ylabel('pdf')
title(['all points, median \Gamma=' num2str(roundx(nanmedian(gam_cham(:)),3))])
freqline(log10(0.2),'k--')
text(log10(0.2),0.5,'0.2')
freqline(nanmean(log10(gam_cham(i4))),'b')

linkaxes([ax1 ax2 ax3])

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
print(fullfile(figdir,'hists_gamma_allVspatch'),'-dpng')

%%
nanmean(gam_all)
nanmean(gam_patch_all)
nanmean(gam_cham(:))

nanmedian(gam_all)
nanmedian(gam_patch_all)
nanmedian(gam_cham(:))

%% Plot gamma vs patch size? 

figure(1);clf
agutwocolumn(0.5)
wysiwyg
scatter(patches(:,4),gam_patch_all,100,'filled','MarkerFaceAlpha',0.1)
grid on
ylim([0 1])
xlim([0 20])
xlabel('patch size [m]','fontsize',16)
ylabel('\Gamma','fontsize',16)

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches'
print(fullfile(figdir,'gamma_Vs_patchsize'),'-dpng')

%%