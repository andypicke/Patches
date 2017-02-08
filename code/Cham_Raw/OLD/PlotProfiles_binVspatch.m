%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% PlotProfiles_binVspatch.m
%
% Compare binned chamelelon profiles to data computed just over patches.
%
% See also CompareChamRawtoNormal.m
%
% * save in separate folders by min patch size ?
%
%--------------------
% 11/29/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

saveplots=0

% load the patch data ( from Compute_N2_dTdz_ChamProfiles_V2.m)
%load(fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
%    'eq14_cham_patches_diffn2dtdzgamma.mat'))

patch_size_min = 0.25
usetemp = 1
saveplots = 1

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/profiles_binVspatch/'
figdir=['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/profiles_binVspatch/minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) ]
%figdir='/Users/Andy/Dropbox/AP_Share_With_JN/profiles_binVspatch'
ChkMkDir(figdir)


load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

% mean depth of each patch
patches.pmn=nanmean([patches.p1(:) patches.p2(:)],2);

for cnum=1:50:3100
    
    try
        
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/mat/eq14_' sprintf('%04d',cnum) '.mat'])
        
        clear ig
        ig=find(patches.cnum==cnum);
        
        %
        figure(1);clf
        set(gcf,'Name',['cast ' num2str(cnum)])
        agutwocolumn(1)
        wysiwyg
        set(gcf,'defaultaxesfontsize',15)
        
        yl=[0 200];
        
        ax1=subplot(221);
        xl=[-6 -2];
        hf=ShadePatchDepths(patches.p1(ig),patches.p2(ig),xl);
        hold on
        plot(log10(avg.N2),avg.P,'.-','color',0.5*[1 1 1])
        axis ij
        grid on
        plot( log10(patches.n2_bin(ig)),patches.pmn(ig),'bd','linewidth',2)
        %plot( log10(patches.nb(ig)),patches.pmn(ig),'kd','linewidth',2)
        plot( log10(patches.n3(ig)),patches.pmn(ig),'kd','linewidth',2)
        xlabel('log_{10}[N^2]')
        xlim(xl)
        ylim(yl)
        ylabel('P')
        %
        ax2=subplot(222);
        xl=[-4 0];
        hf=ShadePatchDepths(patches.p1(ig),patches.p2(ig),xl);
        hold on
        plot(log10(avg.DTDZ_RHOORDER),avg.P,'.-','color',0.5*[1 1 1])
        axis ij
        grid on
        %plot( log10(patches.dtdz2(ig)),patches.pmn(ig),'kd','linewidth',2)
        h1=plot( log10(patches.dtdz_bin(ig)),patches.pmn(ig),'bd','linewidth',2);
        
        h2=plot( log10(patches.dtdz3(ig)),patches.pmn(ig),'kd','linewidth',2);
        legend([h1 h2],'bin','patch','location','south','orientation','horizontal')
        xlabel('log_{10}[dT/dz]')
        xlim(xl)
        ylim(yl)
        %
        ax3=subplot(223);
        xl=[-11 -4];
        hf=ShadePatchDepths(patches.p1(ig),patches.p2(ig),xl);
        hold on
        plot( log10(avg.CHI),avg.P,'.-','color',0.5*[1 1 1])
        axis ij
        grid on
        plot( log10(patches.chi_bin(ig)),patches.pmn(ig),'bd','linewidth',2)
        plot( log10(patches.chi(ig)),patches.pmn(ig),'kd','linewidth',2)
        xlabel('log_{10}[\chi]')
        ylabel('P')
        xlim(xl)
        ylim(yl)
        %
        ax4=subplot(224);
        xl=[-11 -4];
        hf=ShadePatchDepths(patches.p1(ig),patches.p2(ig),xl);
        hold on
        plot( log10(avg.EPSILON),avg.P,'.-','color',0.5*[1 1 1])
        axis ij
        grid on
        plot( log10(patches.eps_bin(ig)),patches.pmn(ig),'bd','linewidth',2)
        plot( log10(patches.eps(ig)),patches.pmn(ig),'kd','linewidth',2)
        xlim(xl)
        xlabel('log_{10}[\epsilon]')
        ylim(yl)
        hh=vline(-8.5,'-');set(hh,'color',0.3*[1 1 1])
        linkaxes([ax1 ax2 ax3 ax4],'y')
        
        if saveplots==1
            print( fullfile(figdir,['cnum_' num2str(cnum) '_N2dtdzChiEps_profs']) , '-dpng' )
        end
        
        % Plot gammas
        
        figure(2);clf
        xl=[0 0.6]
        addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code
        clear gam
        gam=ComputeGamma(avg.N2,avg.DTDZ_RHOORDER,avg.CHI,avg.EPSILON);
        
        hf=ShadePatchDepths(patches.p1(ig),patches.p2(ig),xl);
        plot(gam,avg.P,'o-','color',0.5*[1 1 1])
        hold on
        % plot(patches.gam2(ig),patches.pmn(ig),'kd','linewidth',2)
        plot(patches.gam_bin(ig),patches.pmn(ig),'bd','linewidth',2)
        plot(patches.gam3(ig),patches.pmn(ig),'kd','linewidth',2)
        h1=freqline(nanmean(gam),'--');
        set(h1,'color',0.5*[1 1 1],'linewidth',2)
        %h2=freqline(nanmean(patches.gam2(ig)),'k--');
        %set(h2,'linewidth',2)
        h2=freqline(nanmean(patches.gam3(ig)),'k--');
        set(h2,'linewidth',2)
        axis ij
        grid on
        xlim(xl)
        ylim(yl)
        xlabel('\Gamma')
        ylabel('P')
        title(['cast # ' num2str(cnum)])
        %
        
        if saveplots==1
            print( fullfile(figdir,['cnum_' num2str(cnum) '_gamma_profs']) , '-dpng' )
        end
        
        %        pause
        
    end % try
    
end % cnum

%% Plot mean gamma versus cast #

mn_bin=[]
mn_p1=[]
mn_p2=[]
cnums=[]
for cnum=1:3100
    
    try
        
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/mat/eq14_' sprintf('%04d',cnum) '.mat'])
        
        clear ig
        ig=find(patches.cnum==cnum);
        
        mn_p1 = [mn_p1 ; nanmean(patches.gam2(ig)) ];
        
        clear gam
        gam=ComputeGamma(avg.N2,avg.DTDZ_RHOORDER,avg.CHI,avg.EPSILON);
        mn_bin = [ mn_bin ; nanmean(gam) ];
        
        cnums=[cnums cnum];
        
    end % try
    
end % cnum

%

figure(1);clf
h1=plot(cnums,mn_bin,'.','markersize',9);
hold on
h2=plot(cnums,mn_p1,'r.','markersize',9);
ylim([0 1])
grid on
xlabel('cast #')
ylabel('mean \Gamma')
legend([h1 h2],'binned','patch')
title('mean gamma for each cast')
%%

%%