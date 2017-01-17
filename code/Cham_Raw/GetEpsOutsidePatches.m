%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% GetEpsOutsidePatches.m
%
%
% Get chipod and Chameleon epsilon outside of patches
%
% 1/6/17 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

Params.gamma=0.03

% binned Chameleon epsilon
eps_bin_IN  = [] ;
eps_bin_OUT  = [] ;
% chipod estimates of epsilon using binned data
eps_chipod_binned_IN   = [] ;
eps_chipod_binned_OUT   = [] ;

% load patch info
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/EQ14_raw_patches_minOT_25_usetemp_1.mat')

for cnum=1:3100
    
    try
                        
        % regular chi-pod method on binned data
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax7Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        avg_chipod_binned=avg;clear avg   
        
        clear eps P
        eps = avg_chipod_binned.eps1 ;
        P = avg_chipod_binned.P ;
        
        % go through patches for this profile and nan out any values in
        % patch
        
        clear ipatch Np this_patch
        ipatch=find(new_patch_data(:,1)==cnum) ;
        Np=length(ipatch);
        this_patch = new_patch_data(ipatch,:);
        
                           
        % load binned chameleon profile
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/mat/eq14_' sprintf('%04d',cnum) '.mat'])
        
        clear epscham pcham epschamI
        epscham=avg.EPSILON;
        pcham=avg.P;
        epschamI=interp1(pcham,epscham,P);
        
        for ip=1:Np
           clear iz
           % find ids of data inside patches
           iz=isin(P,[this_patch(ip,2) this_patch(ip,3)]);
           
           % add IN-patch data
           eps_chipod_binned_IN = [eps_chipod_binned_IN ; eps(iz) ] ;
           eps_bin_IN = [ eps_bin_IN ; epschamI(iz) ] ;
           
           eps(iz)=nan;
           epschamI(iz)=nan;
        
        end
        
        % regular chi-pod method epsilons (10m smooth, gam=0.2 etc)
        eps_chipod_binned_OUT  =[eps_chipod_binned_OUT   ;     eps(:)    ];
        eps_bin_OUT = [eps_bin_OUT ; epschamI(:) ];
                       
    end % try
    
end % cnum

% ig=find(~isnan(eps_chipod_binned_OUT));
% eps_chipod_binned_OUT=eps_chipod_binned_OUT(ig);
% 
% ig=find(~isnan(eps_bin_OUT));
% eps_bin_OUT=eps_bin_OUT(ig);

% Save data

PatchInOut = struct('eps_bin_IN',eps_bin_IN,'eps_bin_OUT',eps_bin_OUT,'eps_chipod_binned_IN',...
    eps_chipod_binned_IN,'eps_chipod_binned_OUT',eps_chipod_binned_OUT);
PatchInOut.MakeInfo = ['Made ' datestr(now) ' w/ GetEpsOutsidePatches.m']
PatchInOut.Info='epsilon from chipod (applied to 1-m avg data) and Chaemelon, for depth ranges that are NOT in patches'

save(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data',['PatchesInOutEps_gamma' num2str(Params.gamma*100) '.mat']),'PatchInOut')


%%

saveplots=0

eps_bin(log10(eps_bin)<-8.5)=nan;
%eps_chipod_binned(log10(eps_chipod_binned)<-8.5)=nan;

figure(1);clf
histogram2(log10(eps_bin),log10(eps_chipod_binned),100,'DisplayStyle','Tile')
xlim([-8.5 -4])
ylim([-12 -4])
xvec=linspace(-8.5,-4,100);
hold on
plot(xvec,xvec,'k--','linewidth',2)
plot(xvec,xvec-1,'r--','linewidth',2)
plot(xvec,xvec+1,'r--','linewidth',2)
xlabel('log_{10}[\epsilon cham]')
ylabel('log_{10}[\epsilon \chi]')
title('OUTSIDE patches')


if saveplots==1
figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print(fullfile(figdir,'EpsBin_scatter_OUTSIDEpatches'),'-dpng')
end

figure(2);clf
to_plot= log10( eps_chipod_binned ./ eps_bin ) ;
histogram(to_plot,'Normalization','pdf')
xlim([-4.5 2])
ylabel('pdf')
freqline(nanmean(to_plot))
text(0.5,0.3,['\mu = ' num2str(roundx(nanmean(to_plot),2))])
text(0.5,0.25,['\mu = ' num2str(roundx(nanstd(to_plot),2))])
xlabel('log_{10}[\epsilon_{\chi}/\epsilon_{cham}]')
title('OUTSIDE patches')
grid on

if saveplots==1
figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print(fullfile(figdir,'EpsBin_hist_OUTSIDEpatches'),'-dpng')
end

%%