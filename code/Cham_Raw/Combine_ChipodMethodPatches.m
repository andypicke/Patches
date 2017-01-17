%~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Combine_ChipodMethodPatches.m
%
%
%
%
% Used to be part of CompareProfiles_bin_patch_cham.m
%
%--------------
% 1/6/17 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Combine results from all profiles (so we can make scatter plot etc) and save

clear ; close all

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

whN2dTdz='bulk'
dir1=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches/')

for cnum=1:3100
    
    try
        
        % patch N^2,dTdz w/ constant gamma
        load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128'],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        avg_patchN2dTdz_constGam=avg;clear avg
        
        % patch N^2,dTdz w/ patch gamma
        load( fullfile( dir1, ['N2dTdz_' (whN2dTdz) '_fmax7Hz_respcorr0_fc_99hz_gammaPATCH_nfft_128'],['EQ14_' sprintf('%04d',cnum) 'avg.mat']))
        avg_patchN2dTdzGam=avg;clear avg
                
        % regular chi-pod method on binned data
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm1m_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        %load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/zsm10m_fmax7Hz_respcorr0_fc_99hz_gamma5_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat'])
        avg_chipod_binned=avg;clear avg   
        
        % get these values at patch locations
        eps3=nan*ones(size(avg_patchN2dTdz_constGam.eps1));
        chi3=nan*ones(size(avg_patchN2dTdz_constGam.eps1));
        for ip=1:length(avg_patchN2dTdz_constGam.P)
           [val,I]=nanmin(abs(avg_chipod_binned.P-avg_patchN2dTdz_constGam.P(ip)));
           eps3(ip)=avg_chipod_binned.eps1(I);
           chi3(ip)=avg_chipod_binned.chi1(I);
        end
        
        % regular chi-pod method epsilons (10m smooth, gam=0.2 etc)
        eps_chipod_binned  =[eps_chipod_binned   ; eps3(:)        ];
        chi_chipod_binned  =[chi_chipod_binned   ; chi3(:)        ];

        eps_patchN2dTdz_constGam = [eps_patchN2dTdz_constGam ; avg_patchN2dTdz_constGam.eps1(:) ];
        chi_patchN2dTdz_constGam = [chi_patchN2dTdz_constGam ; avg_patchN2dTdz_constGam.chi1(:) ];
        
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
% save data
sav_name = ['epsilons_N2dTdz_' num2str(whN2dTdz) '_chipodmethods.mat'] ; 
sav_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches'
save(fullfile(sav_dir,sav_name),'AllEps')

%%