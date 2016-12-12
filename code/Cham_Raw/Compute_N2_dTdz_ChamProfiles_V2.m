%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute_N2_dTdz_ChamProfiles_V2.m
%
% Modified to use patches found in Find_Patches_EQ14_Raw.m instead of
% re-finding them here. This way we will have exact same patches to lineup
% with chi and eps values..
%
% Compute N2 and dT/dz for overturns in EQ14 chameleon profiles using a few
% different methods to compare the results.
%
% See also LookAtProfilesOT.m
%
% - Run FindPatches_EQ14_Raw.m 1st to identify patches
% - then run run_eq14_for_PATCHES.m to process cham data over these patches
%
%--------------
% 11/22/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% load patch data (from FindPatches_EQ14_Raw.m)
patch_size_min=0.5  % min patch size
usetemp=1
datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc'
fname=['EQ14_raw_patches_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '.mat']
load(fullfile(datdir,fname))

patches=struct() ;
patches.cnum = new_patch_data(:,1) ;
patches.p1   = new_patch_data(:,2) ;
patches.p2   = new_patch_data(:,3) ;
patches.Lt   = new_patch_data(:,6) ;

% Make empty arrays for results
Npatches=length(patches.p1);
EmpVec=nan*ones(Npatches,1);
patches.n1=EmpVec;
patches.nb=EmpVec;
patches.n3=EmpVec;
patches.n4=EmpVec;

patches.dtdz1=EmpVec;
patches.dtdz2=EmpVec;
patches.dtdz3=EmpVec;

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/

cnum_loaded=-5;
warning off
hb=waitbar(0);

for ip=1:Npatches
    
    waitbar(ip/Npatches,hb)
    
    clear cnum
    cnum=patches.cnum(ip);
    
    % check if this cast already loaded to avoid unnecessary re-loading
    if cnum_loaded~=cnum
        % load chameleon cast
        clear cal cal2 head
        cham_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal';
        load( fullfile( cham_dir, ['eq14_' sprintf('%04d',cnum) '.mat'] ) )
        cnum_loaded=cnum;
    else
        
    end
    
    clear s t p lat
    s=cal.SAL(1:end-1); % (end-1) b/c last 2 values are same;
    t=cal.T1(1:end-1);
    p=cal.P(1:end-1);
    
    clear idot lat1 lat2
    idot=strfind(head.lat.start,'.');
    lat1=str2num(head.lat.start(1:idot-3));
    lat2=str2num(head.lat.start(idot-2:end))/60;
    lat=nanmean([lat1 lat2]);
    
    clear pstarts pstops
    clear iz t_ot s_ot p_ot
    
    iz=isin(p,[ patches.p1(ip) patches.p2(ip) ]);
    t_ot=t(iz);
    s_ot=s(iz);
    p_ot=p(iz);
    
    clear t_sort I
    [t_sort , I]=sort(t_ot,1,'descend');
    
    % compute potential temp
    clear t_pot t_pot_sort
    ptmp_ot=sw_ptmp(s_ot,t_ot,p_ot,0);
    [ptmp_sort , Iptmp]=sort(ptmp_ot,1,'descend');
    
    clear DT dz dTdz
    dT=nanmax(ptmp_ot)-nanmin(ptmp_ot);
    dz=nanmax(p_ot)-nanmin(p_ot);
    dTdz=-dT/dz;
    
    b=p_ot(end) - (1/dTdz)*ptmp_ot(end);
    tvec=linspace(nanmin(ptmp_ot),nanmax(ptmp_ot),100);
    pvec=linspace(nanmin(p_ot),nanmax(p_ot),100);
    yvec=dTdz*pvec - dTdz*b;
    
    % fit a line
    P=polyfit(p_ot,ptmp_sort,1);
    
    % save results
    patches.dtdz1(ip)=dTdz;
    patches.dtdz2(ip)=P(1);
    
    %~~ 'bulk gradient' method from Smyth et al 2001
    % essentially = rms T (btw sorted/unsorted) /  thorpe scale ?
    %    t_rms= sqrt( nanmean(( t_ot - t_sort ).^2) );
    t_rms= sqrt( nanmean(( ptmp_ot - ptmp_sort ).^2) );
    patches.dtdz3(ip)= t_rms / patches.Lt(ip) ;
    
    % Now do similar for density / N^2
    
    clear sgth sgth_ot
    sgth=sw_pden(s,t,p,0);
    sgth_ot=sw_pden(s_ot,t_ot,p_ot,0);
    
    clear sgth_sort I
    [sgth_sort , I]=sort(sgth_ot,1,'ascend');
    
    % try the range/dz method
    drho=nanmax(sgth_ot)-nanmin(sgth_ot);
    dz=nanmax(p_ot)-nanmin(p_ot);
    n2_1=9.81/nanmean(sgth_ot)*drho/dz;
    
    patches.n1(ip)=n2_1;
    
    % fit a line to sgth
    clear P1
    P1=polyfit(p_ot,sgth_sort,1);
    % calculate N^2 from this fit
    clear drhodz n2_2 drho dz n2_3
    drhodz=P1(1);
    n2_2=9.81/nanmean(sgth_ot)*drhodz;
    
    patches.nb(ip)=n2_2;
    
    % Smyth eta al says associated N^2=g*bulk gradient (assuming
    % density controlled by temperature??)
    % I think missing divide by rho_0 also?
    patches.n3(ip)= 9.81 / nanmean(sgth_ot) * t_rms / patches.Lt(ip)    ;
    
    % compute N^2 w/ sw_bfrq
    clear n2
    n2=sw_bfrq(s_ot(I),t_ot(I),p_ot,0.5);
    patches.n4(ip)=nanmean(n2);
    
    
end %ip
delete(hb)

warning on

%


%% Compute gamma using different N^2 and dT/dz values

% first, need to get chi and epsilon for each patch

data_dir=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/',['minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp)])

chi_all=[];
eps_all=[];
P_all=[];
cnum_all=[];

patches.eps=nan*ones(size(patches.p1));
patches.chi=nan*ones(size(patches.p1));

cnums=unique(patches.cnum);

hb=waitbar(0,'compiling patch data from all profiles');

for ic=1:length(cnums)
    waitbar(ic/length(cnums),hb)
    
    clear avg cnum
    
    cnum=cnums(ic);
    
    try
        % load the processed profile w/ chi and eps
        fname=['eq14_' sprintf('%04d',cnum) '.mat'];
        load(fullfile(data_dir,fname))
        
        % now loop over all patches
        clear ig
        ig=find(patches.cnum==cnum);
        for ip=1:length(ig)
            
            clear iz
            iz=isin( avg.P,[patches.p1(ig(ip)) patches.p2(ig(ip))]);
            
            if size(iz)==1
                patches.eps(ig(ip))=avg.EPSILON(iz);
                patches.chi(ig(ip))=avg.CHI(iz);
            end
            
        end %ip
        
    end
    
end % cnum
delete(hb)

%% Exclude values where epsilon is below noise floor
% note this makes significant difference in gamma median
ib=find(log10(patches.eps)<-8.5);
patches.eps(ib)=nan;

%%
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/
% use 'bulk gradient' for both
gam1=ComputeGamma(patches.n1,patches.dtdz1,patches.chi,patches.eps);
% use polyfit for both
gam2=ComputeGamma(patches.nb,patches.dtdz2,patches.chi,patches.eps);
% use range/dz for both
gam3=ComputeGamma(patches.n3,patches.dtdz3,patches.chi,patches.eps);
% plot histograms of each
gam4=ComputeGamma(patches.n4,patches.dtdz2,patches.chi,patches.eps);

figure(1);clf
h1=histogram(log10(gam1(:)));
hold on
histogram(log10(gam2(:)),h1.BinEdges);
histogram(real(log10(gam3(:))),h1.BinEdges);
histogram(real(log10(gam4(:))),h1.BinEdges);
%xlim([-3 1])
%histogram((gam1(:)));xlim([0 1])

iz1=find(gam1>0 & gam1<1);
iz2=find(gam2>0 & gam2<1);
iz3=find(gam3>0 & gam3<1 );%& log10(patches.eps)>-8);
iz4=find(gam4>0 & gam4<1);

figure(2);clf
h1=histogram( gam1(iz1), 40 );
hold on
h2=histogram( gam2(iz2) , h1.BinEdges );
h3=histogram( gam3(iz2) , h1.BinEdges );
h4=histogram( gam4(iz2) , h1.BinEdges );
freqline(nanmedian(gam1))
freqline(nanmedian(gam2))
freqline(nanmedian(gam3))
freqline(nanmedian(gam4))
text(nanmedian(gam1),3500,'\Gamma 1')
text(nanmedian(gam2),3600,'\Gamma 2')
text(nanmedian(gam3),3700,'\Gamma 3')
text(nanmedian(gam4),3800,'\Gamma 4')
legend([h1 h2 h3 h4],'\Gamma 1','\Gamma 2','\Gamma 3','\Gamma 4')
xlim([0 0.5])
grid on
xlabel('\Gamma')
ylabel('count')

print( fullfile( figdir, ['eq14_cham_patch_gammas'] ), '-dpng' )

%%

gams=[ nanmedian(gam1) nanmedian(gam2) nanmedian(gam3) nanmedian(gam4) ;
    nanmedian(gam1(iz1)) nanmedian(gam2(iz2)) nanmedian(gam3(iz3)) nanmedian(gam4(iz4)) ]

%% save results so we don't have to run everything again

patches.gam1=gam1;
patches.gam2=gam2;
patches.gam3=gam3;
patches.gam4=gam4;

save( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
    'eq14_cham_patches_diffn2dtdzgamma.mat'), 'patches' )

%

%%
% I want to compare gamma from binned data to the gamma I compute over
% patches. So for each patch I want to find the closest binned gamma and
% then compare those.

clear ; close all

% load binned chameleon data
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

% load the patch data ( from Compute_N2_dTdz_ChamProfiles_V2.m)
load(fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
    'eq14_cham_patches_diffn2dtdzgamma.mat'))

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/

Npatches=length(patches.cnum)
patches.gam_bin=nan*ones(size(patches.gam1));
patches.n2_bin=nan*ones(size(patches.gam1));
patches.dtdz_bin=nan*ones(size(patches.gam1));
patches.chi_bin=nan*ones(size(patches.gam1));
patches.eps_bin=nan*ones(size(patches.gam1));

for ip=1:Npatches
    clear cnum pbin pmn val I
    cnum=patches.cnum(ip);
    Icham=find(cham.castnumber==cnum);
    pbin=cham.P(:,Icham);
    pmn=nanmean([patches.p1(ip) patches.p2(ip)]);
    [val,I]=nanmin( abs(pbin-pmn));
    
    patches.n2_bin(ip)=cham.N2(I,Icham);
    patches.dtdz_bin(ip)=cham.DTDZ_RHOORDER(I,Icham);
    patches.chi_bin(ip)=cham.CHI(I,Icham);
    patches.eps_bin(ip)=cham.EPSILON(I,Icham);
    
    if log10(cham.EPSILON(I,Icham))>-8.5
        patches.gam_bin(ip)=ComputeGamma(cham.N2(I,Icham),cham.DTDZ_RHOORDER(I,Icham),cham.CHI(I,Icham),cham.EPSILON(I,Icham));
    end
    
end

% save again w/ gam_bin added
save(fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
    'eq14_cham_patches_diffn2dtdzgamma.mat'),'patches')


%% Plot 2D hist/scatter of different dt/dz methods

figure(1);clf
agutwocolumn(1)
wysiwyg
set(gcf,'defaultaxesfontsize',15)

subplot(211)
%loglog(abs(dtdz1(:)),abs(dtdz2(:)),'.')
histogram2(log10(abs(patches.dtdz1(:))),log10(abs(patches.dtdz2(:))),50,'displaystyle','tile')
xlabel('log_{10}[dtdz1]')
ylabel('log_{10}[dtdz2]')
%grid on
grid on

subplot(212)
%loglog(abs(dtdz1(:)),abs(dtdz2(:)),'.')
histogram2(log10(abs(patches.dtdz1(:))),log10(abs(patches.dtdz3(:))),50,'displaystyle','tile')
xlabel('log_{10}[dtdz1]')
ylabel('log_{10}[dtdz3]')
%grid on
grid on

figdir='/Users/Andy/Cruises_Research/ChiPod/Analyses/Patch_n2_dTdz'
print( fullfile( figdir, ['eq14_cham_patch_dTdzs_scatter'] ), '-dpng' )

%% Make histogram of different dt/dzs

figure(1);clf
agutwocolumn(0.7)
wysiwyg
set(gcf,'defaultaxesfontsize',14)

h1=histogram(log10(abs(patches.dtdz1(:))),'DisplayStyle','stair')%,'edgecolor','none');
hold on
h2=histogram(log10(abs(patches.dtdz2(:))),h1.BinEdges,'DisplayStyle','stair')%,'edgecolor','none');
h3=histogram(log10(abs(patches.dtdz3(:))),h1.BinEdges,'DisplayStyle','stair')%,'edgecolor','none');
grid on
xlim([-4 -0.5])
xlabel('log_{10}[dT/dz]')
ylabel('count')
legend([h1 h2 h3],'dtdz1','dtdz2','dtdz3')

figdir='/Users/Andy/Cruises_Research/ChiPod/Analyses/Patch_n2_dTdz'
print( fullfile( figdir, ['eq14_cham_patch_dTdzs'] ), '-dpng' )

%% Plot 2D histograms of different N2 methods

figure(2);clf
agutwocolumn(1)
wysiwyg
set(gcf,'defaultaxesfontsize',15)

xl=[-5.5 -3]
yl=[-6.75 -2.5]
m=3
n=1

subplot(m,n,1)
%loglog(n1(:),nb(:),'.')
histogram2(real(log10(patches.n1(:))),real(log10(patches.nb(:))),50,'displaystyle','tile')
grid on
xlim(xl)
ylim(yl)

xlabel('log_{10}[N^2_1]')
ylabel('log_{10}[N^2_2]')

subplot(m,n,2)
%loglog(n1(:),n3(:),'.')
histogram2(real(log10(patches.n1(:))),real(log10(patches.n3(:))),50,'displaystyle','tile')
grid on
xlabel('log_{10}[N^2_1]')
ylabel('log_{10}[N^2_3]')
xlim(xl)
ylim(yl)

subplot(m,n,3)
%loglog(n1(:),n3(:),'.')
histogram2(real(log10(patches.n1(:))),real(log10(patches.n4(:))),50,'displaystyle','tile')
grid on
xlabel('log_{10}[N^2_1]')
ylabel('log_{10}[N^2_4]')
xlim([xl])
ylim([yl])

print( fullfile( figdir, ['eq14_cham_patch_N2s_scatter'] ), '-dpng' )
%subplot(313)

%% Plot histogram of different N2s

figure(2);clf
agutwocolumn(0.7)
wysiwyg
set(gcf,'defaultaxesfontsize',14)

h1=histogram(real(log10(patches.n1(:)))) ;
hold on
h2=histogram(real(log10(patches.nb(:))),h1.BinEdges) ;
h3=histogram(real(log10(patches.n3(:))),h1.BinEdges) ;
h4=histogram(real(log10(patches.n4(:))),h1.BinEdges) ;
freqline(nanmedian(h1.Data))
freqline(nanmedian(h2.Data))
freqline(nanmedian(h3.Data))
freqline(nanmedian(h4.Data))
grid on
xlim([-6 -2.5])
xlabel('log_{10}[N^2]')
ylabel('count')
legend([h1 h2 h3 h4],'1  ','2 ','3 ','4 ','location','best')

print( fullfile( figdir, ['eq14_cham_patch_N2s'] ), '-dpng' )

%%

figure(1);clf
agutwocolumn(1)
wysiwyg

ax1=subplot(221);
h1=histogram(real(log10(patches.nb(:))),'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(patches.n2_bin(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
xlim([-6.5 -2.5])
grid on
xlabel('log_{10}[N^2]')
legend([h1 h2],'patch','bin','location','best')
freqline(nanmedian(h1.Data),'b--')
freqline(nanmedian(h2.Data),'r--')

ax2=subplot(222);
h1=histogram(real(log10(patches.dtdz2(:))),'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(patches.dtdz_bin(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
xlim([-5 0])
grid on
xlabel('log_{10}[dT/dz]')
legend([h1 h2],'patch','bin','location','best')
freqline(nanmedian(h1.Data),'b--')
freqline(nanmedian(h2.Data),'r--')

ax3=subplot(223);
h1=histogram(real(log10(patches.chi(:))),'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(patches.chi_bin(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
xlim([-12 -3])
grid on
xlabel('log_{10}[\chi]')
legend([h1 h2],'patch','bin','location','best')
freqline(nanmedian(h1.Data),'b--')
freqline(nanmedian(h2.Data),'r--')

ax4=subplot(224);
h1=histogram(real(log10(patches.eps(:))),'Normalization','pdf','edgecolor','none');
hold on
h2=histogram(real(log10(patches.eps_bin(:))),h1.BinEdges,'Normalization','pdf','edgecolor','none');
%xlim([-12 -3])
grid on
xlabel('log_{10}[\epsilon]')
legend([h1 h2],'patch','bin','location','best')
freqline(nanmedian(h1.Data),'b--')
freqline(nanmedian(h2.Data),'r--')

figdir='/Users/Andy/Cruises_Research/ChiPod/Analyses/Patch_n2_dTdz'
print( fullfile( figdir, ['eq14_cham_gamma_binVspatch_hists'] ), '-dpng' )

%%

ib=find(patches.gam_bin>1);

figure(1);clf
%scatter(log10(patches.gam_bin(:)),log10(patches.gam3(:)),'o','filled','MarkerFaceAlpha',0.05)
scatter(patches.gam_bin(:), patches.gam3(:),'o','filled','MarkerFaceAlpha',0.1) ; xlim([0 0.5]) ; ylim([0 0.5])
%xlim([-3 1]);ylim([-3 1])
xlabel('\Gamma bin')
ylabel('\Gamma patch')

% ig=find(~isnan(patches.gam_bin) & ~isnan(patches.gam2));
% P=polyfit(patches.gam_bin(ig),patches.gam3(ig),1);
% hold on
% plot(0:0.001:1,polyval(P,0:0.001:1),'r','linewidth',2)
grid on

figdir='/Users/Andy/Cruises_Research/ChiPod/Analyses/Patch_n2_dTdz'
print( fullfile( figdir, ['eq14_cham_gamma_binVspatch_scatter'] ), '-dpng' )


%%

gam_bin_all=ComputeGamma(cham.N2(:),cham.DTDZ_RHOORDER(:),cham.CHI(:),cham.EPSILON(:));

ig1=find(patches.gam_bin<=1);
ig2=find(patches.gam2<=1);

figure(2);clf
h1=histogram((patches.gam_bin(ig1)),'Normalization','pdf');
hold on
h2=histogram((patches.gam2(ig2)),'Normalization','pdf');
%histogram(log10(gam_bin_all(:)),'Normalization','pdf')
freqline(nanmedian(h1.Data))
freqline(nanmedian(h2.Data))
xlim([0 0.5])
grid on
xlabel('\Gamma','fontsize',15)
ylabel('pdf','fontsize',15)
legend([h1 h2],'bin','patch')
nanmedian(patches.gam_bin(:))
nanmedian(patches.gam3(:))

figdir='/Users/Andy/Cruises_Research/ChiPod/Analyses/Patch_n2_dTdz'
print( fullfile( figdir, ['eq14_cham_gamma_binVspatch'] ), '-dpng' )


%%

figure(2);clf
h1=histogram(log10(patches.gam_bin(:)),'Normalization','pdf');
hold on
h2=histogram(log10(patches.gam2(:)),'Normalization','pdf');
%histogram(log10(gam_bin_all(:)),'Normalization','pdf')
freqline(nanmedian(h1.Data),'b--')
freqline(nanmedian(h2.Data),'r--')
%freqline(nanmean(h1.Data),'b')
%freqline(nanmean(h2.Data),'r')
xlim([-4 2])
grid on
xlabel('log_{10}[\Gamma]')
ylabel('pdf')
legend([h1 h2],'bin','patch')
nanmedian(patches.gam_bin(:))
nanmedian(patches.gam3(:))
hf=freqline(log10(0.2),'k--')
set(hf,'linewidth',2)
print( fullfile( figdir, ['eq14_cham_LOGgamma_binVspatch'] ), '-dpng' )

%figdir='/Users/Andy/Cruises_Research/ChiPod/Analyses/Patch_n2_dTdz'
%print( fullfile( figdir, ['eq14_cham_gamma_binVspatch'] ), '-dpng' )

%%
