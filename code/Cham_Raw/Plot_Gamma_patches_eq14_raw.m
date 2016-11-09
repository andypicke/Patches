%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_Gamma_patches_eq14_raw.m
%
%
%
%
%-----------------
% 11/09/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

min_patch_size=2

data_dir=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/',['minOT_' num2str(10*min_patch_size)]) 

n2_all=[];
dtdz_all=[];
chi_all=[];
eps_all=[];
P_all=[];

for cnum=1:3100
    
    try
        fname=['eq14_' sprintf('%04d',cnum) '.mat']
        load(fullfile(data_dir,fname))
        
        P_all = [P_all ; avg.P(:) ] ;
        n2_all = [n2_all ; avg.N2(:)] ;
        dtdz_all = [dtdz_all ; avg.dTdz(:) ] ;
        chi_all = [chi_all ; avg.CHI(:) ] ;
        eps_all = [eps_all ; avg.EPSILON(:) ];
        
    end
    
end % cnum

%%

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code

gam=ComputeGamma(n2_all,dtdz_all,chi_all,eps_all);

%ig=find(gam<2 );
%ig=find(P_all>60 & P_all<150 & log10(eps_all)>-8.5);
ig=find( log10(eps_all)>-8.5);


figure(1);clf
h1=histogram(gam(ig),0:0.1:nanmax(gam));
freqline(nanmedian(h1.Data) )
xlim([0 2])
xlabel('\Gamma','fontsize',16)
title(['min OT size=' num2str(min_patch_size) 'm'])

%%

nanmedian(gam)
nanmedian(gam(ig))

%% How does the 'bulk' N2 compare to N2 from the overturns code?


load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/EQ14_raw_patches_minOT_1_join_0_sep_50.mat')


%%

