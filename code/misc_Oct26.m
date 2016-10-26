%% First let's look at gamma as computed from 1m binned chameleon data

clear ; close all

% load combined TIWE data
load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches','tiwe_comb_AP.mat') )
cham=tiwe;

% compute gamma using all data points
gam_cham=cham.N2 .* cham.CHI ./2 ./ cham.EPSILON ./ (cham.DTDZ.^2);

% find points not in mixed layer or near bottom of profiles
i4=find(cham.P>60 & cham.P<150);

%
numel(i4)/numel(cham.P)

%%
figure(1);clf
h1=histogram(real(log10(gam_cham(:))),'edgecolor','none','Normalization','pdf');
hold on
h2=histogram(log10(gam_cham(i4)),'edgecolor','none','Normalization','pdf');
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

nanmedian(gam_cham(:))
nanmedian(gam_cham(i4))


