%~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Combine_ChipodMethodPatches.m
%
% This script combines (1) chi-pod estimates at patches using constant
% gamma, (2) chi-pod estimates at patches using actual gammas , (3) chi-pod
% estimates using binned data, and (4) chameleon estimates at patches, for
% all profiles and patches so we can compare them etc.
%
% - chi-pod estimates for patches are made in : ComputeChi_Chameleon_Eq14_PATCHES.m
%   - Chameleon chi/eps values at patch locations are also saved as
%   'eps_bin', 'eps_patch' etc
% - chi-pod estimates for binned profiles are made w/ ComputeChi_Chameleon_Eq14.m
%
% Produces a structure 'AllEps'
%
% Used to be part of CompareProfiles_bin_patch_cham.m
%
%--------------
% 1/6/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Combine results from all profiles (so we can make scatter plot etc) and save

clear ; close all

whN2dTdz='bulk'
%whN2dTdz='bulk2'
Params.gamma=0.2;

patch_size_min = 0.25
usetemp = 1

savedata=1;

% pre-allocate empty arrays
eps_patchN2dTdz_constGam = [] ;
eps_patchN2dTdzGam    = [] ;
eps_bin  = [] ;
eps_patch= [] ;
eps_chipod_binned    = [] ;

chi_patchN2dTdz_constGam = [] ;
chi_patchN2dTdzGam    = [] ;
chi_bin  = [] ;
chi_patch= [] ;
chi_chipod_binned    = [] ;

dir1=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches/')

for cnum=1:3100
    
    try
        
        % patch N^2,dTdz w/ constant gamma
        load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax7Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        avg_patchN2dTdz_constGam=avg;clear avg
        
        % patch N^2,dTdz w/ patch gamma
        load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax7Hz_respcorr0_fc_99hz_gammaPATCH_nfft_128_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        avg_patchN2dTdzGam=avg;clear avg
                
        % regular chi-pod method on binned data
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax7Hz_respcorr0_fc_99hz_gamma' num2str(Params.gamma*100) '_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        %load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma5_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        avg_chipod_binned=avg;clear avg   
        
        % get these (binned) values at patch locations
        eps3=nan*ones(size(avg_patchN2dTdz_constGam.eps1));
        chi3=nan*ones(size(avg_patchN2dTdz_constGam.eps1));
        for ip=1:length(avg_patchN2dTdz_constGam.P)
           [val,I]=nanmin(abs(avg_chipod_binned.P-avg_patchN2dTdz_constGam.P(ip)));
           eps3(ip)=avg_chipod_binned.eps1(I);
           chi3(ip)=avg_chipod_binned.chi1(I);
        end
        
        % regular chi-pod method epsilons (1m smooth, gam=0.2 etc)
        eps_chipod_binned  =[eps_chipod_binned   ; eps3(:)        ];
        chi_chipod_binned  =[chi_chipod_binned   ; chi3(:)        ];

        % chi-pod estimates at patches, using constant gamma
        eps_patchN2dTdz_constGam = [eps_patchN2dTdz_constGam ; avg_patchN2dTdz_constGam.eps1(:) ];
        chi_patchN2dTdz_constGam = [chi_patchN2dTdz_constGam ; avg_patchN2dTdz_constGam.chi1(:) ];
        
        % chi-pod estimates at patches, using actual patch gammas
        eps_patchN2dTdzGam       = [eps_patchN2dTdzGam       ; avg_patchN2dTdzGam.eps1(:)       ];
        chi_patchN2dTdzGam       = [chi_patchN2dTdzGam       ; avg_patchN2dTdzGam.chi1(:)       ];
        
        % chameleon epsilon binned data at patch locations
        eps_bin=[eps_bin ; avg_patchN2dTdz_constGam.eps_bin(:)];
        chi_bin=[chi_bin ; avg_patchN2dTdz_constGam.chi_bin(:)];
                
        % chameleon epsilon values computed over patches only
        eps_patch=[eps_patch ; avg_patchN2dTdzGam.eps_patch(:) ] ;
        chi_patch=[chi_patch ; avg_patchN2dTdzGam.chi_patch(:) ] ;
        
    end % try
    
end % cnum

%
ib=find(log10(eps_bin)<-8.5);
eps_bin(ib)=nan;

AllEps=struct('eps_chipod_binned',eps_chipod_binned,'eps_patchN2dTdz_constGam',...
    eps_patchN2dTdz_constGam,'eps_patchN2dTdzGam',eps_patchN2dTdzGam,...
    'eps_bin',eps_bin,'eps_patch',eps_patch,'chi_chipod_binned',chi_chipod_binned,'chi_patchN2dTdz_constGam',...
    chi_patchN2dTdz_constGam,'chi_patchN2dTdzGam',chi_patchN2dTdzGam,...
    'chi_bin',chi_bin,'chi_patch',chi_patch)
AllEps.MakeInfo=['Made ' datestr(now) ' w/ Combine_ChipodMethodPatches.m']

if savedata==1
% save data
sav_name = ['epsilons_N2dTdz_' num2str(whN2dTdz) '_otmin' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_chipodmethods.mat'] ; 
sav_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches'
save(fullfile(sav_dir,sav_name),'AllEps')
end

%%

figure(1);clf
histogram2( log10(AllEps.chi_bin(:)), log10(AllEps.chi_patch(:)), 100, 'DisplayStyle','Tile')
xlim([-12 -3])
ylim([-12 -3])

%%

figure(1);clf
histogram2( log10(AllEps.eps_bin(:)), log10(AllEps.eps_patch(:)), 100, 'DisplayStyle','Tile')
xlim([-8.5 -5])
ylim([-8.5 -5])


%%