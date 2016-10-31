%~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_tiwe_patches_jim_bill.m
%
% Compare Patch Data from FLUX91 from Jim and Bill
%
%
%--------------------
% 10/27/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Flux91/events_FLX91.mat')

load('/Users/Andy/Cruises_Research/ChiPod/Flux91/flx91_patch_out_nov94.mat')

numel(A.chi)
numel(flx91.chi)

% compute gamma from both
gam1= ComputeGamma(A.N2,A.tgrad,A.chi,A.eps);
gam2= ComputeGamma(flx91.bv.^2,flx91.dtdz,flx91.chi,flx91.eps);

%%

figure(1);clf
h1=histogram(log10(gam1(:)),'Normalization','pdf');
hold on
h2=histogram(log10(gam2(:)),'Normalization','pdf');
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')
freqline(log10(0.2),'k--')
xlabel('log_{10}\Gamma')
ylabel('pdf')
legend([h1 h2],'Bill patches','Jim Patches')
grid on

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures'
print( fullfile( figdir,'flx91_gam_hist_compare_bill_jim'), '-dpng')


%% Compare distributions of n2,dtdz,chi,eps

[h h1 h2]=Plot4hist_2(A.N2,A.tgrad,A.chi,A.eps,flx91.bv.^2,flx91.dtdz,flx91.chi,flx91.eps)
legend([h1 h2],'bill','jim')

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures'
print( fullfile( figdir,'flx91_hist_compare_bill_jim'), '-dpng')

%%