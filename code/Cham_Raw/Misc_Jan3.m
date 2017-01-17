%~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Misc_Jan3.m
%
%
%
%%

clear ; close all

cnum=18

BaseDir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/'

load(fullfile(BaseDir,'mfiles','Patches','data','ChipodPatches',['N2dTdz_bulk_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat']) )
avg_patch=avg; clear avg
load(fullfile(BaseDir,'Data','chameleon','processed','Cstar=0_032','mat',['eq14_' sprintf('%04d',cnum) '.mat'] ) )
avg_bin=avg;clear avg
load(fullfile(BaseDir,'Data','Cham_proc_AP','cal',['eq14_' sprintf('%04d',cnum) '.mat'] ) )

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/
n2=sw_bfrq(cal.SAL,cal.T1,cal.P,0.3);

figure(1);clf
agutwocolumn(1)
wysiwyg
semilogx(n2,cal.P(1:end-1),'.','color',0.75*[1 1 1])
hold on
hp=semilogx(avg_patch.N2,avg_patch.P,'ko')
hb=semilogx(avg_bin.N2,avg_bin.P)
legend(hp,'patch')
%semilogx(avg.N2,avg.P,'-')
axis ij
grid on

%semilogx(smooth(n2,25),cal.P(1:end-1))
%semilogx(medfilt1(n2),cal.P(1:end-1),'g.')
xlim([1e-6 1e-2])
ylabel('pressure [db]')
xlabel('log_{10}N^2')

%% 1st check if raw n2 (not smoothed, computed from cal data) is close to n2 from patches
% will find data points in each patch and average to compare them
% so need to loop over each patch

clear ; close all

cnum=68

BaseDir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/'

load(fullfile(BaseDir,'Data','Cham_proc_AP','cal',['eq14_' sprintf('%04d',cnum) '.mat'] ) )

load(fullfile(BaseDir,'mfiles','Patches','data','ChipodPatches',['N2dTdz_bulk_fmax7Hz_respcorr0_fc_99hz_gamma20_nfft_128/EQ14_' sprintf('%04d',cnum) 'avg.mat']) )
avg_patch=avg; clear avg

load(fullfile(BaseDir,'Data','chameleon','processed','Cstar=0_032','mat',['eq14_' sprintf('%04d',cnum) '.mat'] ) )
avg_bin=avg;clear avg

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/
n2=sw_bfrq(cal.SAL,cal.T1,cal.P,0.3);
n2=[n2 ; nan];
n2_med=medfilt1(n2);

% load patches data
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/EQ14_raw_patches_minOT_25_usetemp_1.mat')

% find patch data for this cast
igcast=find(new_patch_data(:,1)==cnum);
clear patches
patches=new_patch_data(igcast,:);

Np=length(igcast)
%
n2_p=nan*ones(Np,1);
n2_p_med=n2_p;

n2_bin=n2_p;
% loop over each patch
for ip=1:Np
    
    clear iz
    iz=isin(cal.P,[patches(ip,2),patches(ip,3)]) ;
    n2_p(ip)=nanmean(n2(iz)) ;
    n2_p_med(ip)=nanmean(n2_med(iz));
    
    % find binned values for patch locations
    clear iz
    iz=isin(avg_bin.P,[patches(ip,2),patches(ip,3)] ) ;
    n2_bin(ip)=nanmean(avg_bin.N2(iz)) ;
    
end

%
figure(1);clf
hp=semilogx(avg_patch.N2,avg_patch.P,'o')
hold on
h=semilogx(n2_p,avg_patch.P,'d')
h=semilogx(n2_p_med,avg_patch.P,'p')
axis ij
grid on

%
figure(1);clf
hraw=loglog(avg_patch.N2,n2_p,'o')
hold on
hb=loglog(avg_patch.N2,n2_bin,'d')
%hp=loglog(avg_patch.N2,n2_p_med,'o')
grid on
xlabel('patch values')
legend([hraw hb],'n2raw','n2bin','location','best')
xvec=linspace(1e-6,1e-2,100);
loglog(xvec,xvec,'k--')
xlim([1e-6 1e-2])
ylim([1e-6 1e-2])


%% Compile patch vs bin n2 for all profiles and see if there is a consistent bias? 

% Compared to the binned values, the patch values of N2 tend to be larger,
% and patch values of dT/dz tend to be smaller. This would increase the
% ratio N^2/dT/dz

clear ; close all
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/eq14_cham_minOT_25_usetemp_1_patches_diffn2dtdzgamma.mat')

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(311)
%loglog(patches.n2_bulk,patches.n2_bin,'.')
histogram2(real(log10(patches.n2_bin)),real(log10(patches.n2_bulk)),100,'DisplayStyle','Tile')
hold on
xvec=linspace(-6,-3,100);
loglog(xvec,xvec,'k--')
xlim([-6 -3])
ylim([-6 -2.5])
xlabel('N^2 bin')
ylabel('N^2 bulk')

subplot(312)
%loglog(patches.n2_bulk,patches.n2_bin,'.')
histogram2(real(log10(patches.dtdz_bin)),real(log10(patches.dtdz_bulk)),100,'DisplayStyle','Tile')
hold on
xvec=linspace(-5,0,100);
loglog(xvec,xvec,'k--')
xlim([-4.5 0])
ylim([-4 -0.5])
xlabel('T_z bin')
ylabel('T_z bulk')

% check if the ratio of N^2/dTdz^2 differs for binned vs patch computations ?
% yes, the ratio is higher for the patch values computed

%clear

%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/eq14_cham_minOT_25_usetemp_1_patches_diffn2dtdzgamma.mat')

rat_bin=patches.n2_bin ./ patches.dtdz_bin.^2 ;
rat_bulk=patches.n2_bulk ./ patches.dtdz_bulk.^2 ;
% 
% %figure(1);clf
% subplot(313)
% hbin=histogram(real(log10(rat_bin)),'Normalization','pdf','Edgecolor','none')
% hold on
% hbulk=histogram(real(log10(rat_bulk)),hbin.BinEdges,'Normalization','pdf','Edgecolor','none')
% xlim([-3 3])
% grid on
% legend([hbin hbulk],'bin','bulk')
% xlabel('log_{10}[N^2/<dTdz>^2]')
% ylabel('pdf')

%
%figure(2);clf
subplot(313)
histogram2( real(log10(rat_bin)), log10(rat_bulk) ,'DisplayStyle','Tile')
xvec=linspace(-4,4,100);
hold on
plot(xvec,xvec,'k--')
xlim([-4 4])
ylim([-2 2])
xlabel('log_{10}[N^2/<dTdz>^2] (bin)')
ylabel('log_{10}[N^2/<dTdz>^2] (patch)')

%%
figdir=

%% what percent of water column do patches cover?

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/eq14_cham_minOT_25_usetemp_1_patches_diffn2dtdzgamma.mat')

patches.dp=patches.p2-patches.p1;
cnums=unique(patches.cnum);

perc=nan*ones(1,length(cnums));

for ic=1:length(cnums)

    clear cnum ig ztot_patches
    cnum=cnums(ic);
    ig=find(patches.cnum==cnum) ;    
    
    ztot_patches=sum(patches.dp(ig));
    perc(ic)=ztot_patches/200*100;

end

%
figure(1);clf
histogram(perc)

%%