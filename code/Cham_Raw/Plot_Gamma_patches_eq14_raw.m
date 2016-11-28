%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot_Gamma_patches_eq14_raw.m
%
% Plot eq14 Chameleon data computed for patches.
%
% Patches are computing in FindPatches_EQ14_Raw.m
% Data (chi,eps,N2,dt/dz etc) over patches are computed
% in run_eq14_for_PATCHES
%
%
%-----------------
% 11/09/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

min_patch_size=0.5
usetemp=1

data_dir=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/',['minOT_' num2str(10*min_patch_size) '_usetemp_' num2str(usetemp)]) 

n2_all=[];
n2_OT_all=[];
dtdz_all=[];
chi_all=[];
eps_all=[];
P_all=[];

patch_data=[];

hb=waitbar(0,'compiling patch data from all profiles');
for cnum=4:3100
    waitbar(cnum/3100,hb)
    clear avg
    try
        fname=['eq14_' sprintf('%04d',cnum) '.mat'];
        load(fullfile(data_dir,fname))
        
        if length(avg.N2)~=length(avg.n2_OT)
            disp(num2str(cnum))
        end
        
        P_all     = [P_all     ; avg.P(:)    ] ;
        n2_all    = [n2_all    ; avg.N2(:)   ] ;
        n2_OT_all = [n2_OT_all ; avg.n2_OT(:)] ;
        dtdz_all  = [dtdz_all  ; avg.dTdz(:) ] ;
        chi_all   = [chi_all   ; avg.CHI(:)  ] ;
        eps_all   = [eps_all   ; avg.EPSILON(:) ];
        
    end
    
end % cnum
delete(hb)

%% Compare N2 from OT code to N2 from Chameleon code

figure(1);clf
loglog(n2_all(:),n2_OT_all(:),'.')
grid on
xlabel('cham (''bulk'')')
ylabel('OT')
title(['eq14 N^2 over patches'])

%%

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code

gam1=ComputeGamma(n2_all,dtdz_all,chi_all,eps_all);
gam2=ComputeGamma(n2_OT_all,dtdz_all,chi_all,eps_all);

%%

ig1=find( log10(eps_all)>-8.5 & gam1<5);
ig2=find( log10(eps_all)>-8.5 & gam2<5);

figure(1);clf
h1=histogram(gam1(ig1),0:0.025:5);
hold on
h2=histogram(gam2(ig2),h1.BinEdges);
freqline(nanmedian(gam1) ,'b')
freqline(nanmedian(gam2), 'r' )
xlim([0 1])
xlabel('\Gamma','fontsize',16)
ylabel('count','fontsize',16)
title(['min OT size=' num2str(min_patch_size) 'm'])
legend([h1 h2],'cham N2','OT N2')

%%

figpath='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
figname=['eq14_gamma_patches_raw_minOT_' num2str(10*min_patch_size)]
print( fullfile( figpath, figname ), '-dpng')

%%
clc
nanmedian(gam1)
nanmedian(gam1(ig1))

nanmedian(gam2)
nanmedian(gam2(ig2))

%%

