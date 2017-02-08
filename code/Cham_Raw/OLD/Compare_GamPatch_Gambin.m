%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_GamPatch_Gambin.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/eq14_cham_minOT_25_usetemp_1_patches_diffn2dtdzgamma.mat')

figure(1);clf
histogram2(log10(patches.gam_bulk(:)), log10(patches.gam_bin(:)) ,'DisplayStyle','Tile' )
xlim([-3 1])
ylim([-3 1])
xvec=linspace(-3,1,100);
hold on
plot(xvec,xvec,'k--','linewidth',2)
xlabel('\Gamma patch')
ylabel('\Gamma bin')

%%
figure(2);clf
to_plot = log10(patches.gam_bin ./ patches.gam_bulk) ;
histogram( to_plot , 'Normalization','pdf')
xlim([-3 3])
freqline(nanmean(to_plot) )
xlabel( 'log_{10}[\Gamma_{bin}/\Gamma_{patch}]')
ylabel('pdf')
title('ratio of \Gamma* to \Gamma in patches, EQ14')
grid on
text(1,0.6,['\mu = ' num2str(roundx(nanmean(to_plot),2)) ], 'fontsize',16)
text(1,0.55,['\sigma = ' num2str(roundx(nanstd(to_plot),2)) ], 'fontsize',16)
%%
10^(nanmean(to_plot))

%% Find binned gamma *outside* of patches



%%