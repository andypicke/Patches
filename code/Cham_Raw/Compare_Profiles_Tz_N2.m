%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compare_Profiles_Tz_N2.m
%
% formerly named Misc_Jan27.m
%
% See also LookAtN2Profiles.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% (1)

% Using 1m binned profiles, compare my method of computing N^2 and T_z to
% the Chameleon method (fitting line over 1m segments)

clear ; close all

for cnum=1:50:3100
    
    try
        addpath /Users/Andy/Cruises_Research/seawater_ver3_2/
        
        % load raw chameleon profile (not binned)
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal/EQ14_' sprintf('%04d',cnum) '.mat'])
        clear cal
        cal=cal2;
        clear cal2
        
        % load binned chameleon profile ('avg') to compare with N2 and T_z computed from
        % linear fits over 1m bins
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/mat/eq14_' sprintf('%04d',cnum) '.mat'])
        
        % compute N^2 and T_z from binned profile
        N2_bin=sw_bfrq(avg.SAL,avg.T1,avg.P,0);
        N2_bin=[N2_bin(:) ; nan ] ;
        
        % negative sign so T_z>0 for temp increasing upward
        dtdz_bin= -diffs(avg.T1) ./ diffs(avg.P) ;
                
        figure(1);clf
        h1=plot( log10(avg.N2), avg.P);
        hold on
        h2=plot( log10(N2_bin), avg.P);
        axis ij
        grid on
        legend([h1 h2],'cham','my','location','best')
        
        figure(2);clf
        loglog(avg.N2(:),N2_bin(:),'d')
        xvec=linspace(1e-7,1e-2,100);
        hold on
        loglog(xvec,xvec,'k--')
        grid on
        axis tight
        xlim([1e-7 1e-2])
        ylim([1e-7 1e-2])
        
        pause(1)
    end
    
end

%%

% Ok, now how do N^2 and T_z from 1m profiles compare to those computed w/
% hi-res data?

clear ; close all

cnum=2549

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/

% load raw chameleon profile (not binned)
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal/EQ14_' sprintf('%04d',cnum) '.mat'])
clear cal
cal=cal2;
clear cal2

% load binned chameleon profile ('avg') to compare with N2 and T_z computed from
% linear fits over 1m bins
load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/mat/eq14_' sprintf('%04d',cnum) '.mat'])

% sorted temperature
[cal.T1sort,I] = sort(cal.T1,1,'descend');

% smooth salinity a little bit
cal.SAL=smooth(cal.SAL,20);

% compute pot. dens.
cal.sgth=sw_pden(cal.SAL,cal.T1,cal.P,0);

% sort pot. dens
[val,Irho] = sort(cal.sgth,1,'ascend');
%n2=sw_bfrq(cal.SAL,cal.T1,cal.P,0);
n2=sw_bfrq(cal.SAL(Irho),cal.T1(Irho),cal.P,0);
n2=[n2(:) ; nan ];

%dtdz=-diffs(cal.T1) ./ diffs(cal.P) ;
dtdz=-diffs(cal.T1sort) ./ diffs(cal.P) ;

% interpolate hi-res N2,T_z to 1m grid to compare

[C,IA,IC]=unique(cal.P);
n2i=interp1(cal.P(IA),n2(IA),avg.P);
dtdzi=interp1(cal.P(IA),dtdz(IA),avg.P);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Plotting

figure(1);clf
set(gcf,'defaultaxesfontsize',15)
set(gcf,'Name',['profile ' num2str(cnum)])
agutwocolumn(1)
wysiwyg

ax1=subplot(221);
hraw=plot(cal.T1,cal.P,'.','color',0.5*[1 1 1]);
hold on
hbin=plot(avg.T1,avg.P);
axis ij
axis tight
grid on
xlabel('T','fontsize',16)
ylabel('P [DB]','fontsize',16)
legend([hraw hbin],'raw','1mbin','location','best')

ax2=subplot(222);
plot(real(log10(dtdz)),cal.P,'.','color',0.7*[1 1 1])
hold on
plot(real(log10(avg.DTDZ_RHOORDER)),avg.P,'k')
plot(log10(dtdzi),avg.P,'r')
axis ij
axis tight
grid on
xlabel('log_{10}[T_z]','fontsize',16)

ax3=subplot(223);
plot(cal.sgth,cal.P,'.','color',0.5*[1 1 1])
hold on
plot(avg.SIGMA+1000,avg.P)
axis ij
axis tight
grid on
xlabel('pot. dens.','fontsize',16)
ylabel('P [DB]','fontsize',16)

ax4=subplot(224);
h1=plot(log10(n2),cal.P,'.','color',0.7*[1 1 1]);
hold on
h2=plot(log10(avg.N2),avg.P,'k');
h3=plot(log10(n2i),avg.P,'r');
axis ij
axis tight
grid on
xlabel('log_{10}[N^2]','fontsize',16)
legend([h1 h2 h3],'raw','1m','rawinterp','location','best')

linkaxes([ax1 ax2 ax3 ax4],'y')

% save plot
SetNotesFigDir
figname=['eq14_prof_' num2str(cnum) '_Tz_N2_raw_bin']
print( fullfile( NotesFigDir , figname ) , '-dpng' )

%% Scatterplot binned vs raw-interp

figure(2);clf
set(gcf,'defaultaxesfontsize',15)
set(gcf,'Name',['profile ' num2str(cnum)])
agutwocolumn(1)
wysiwyg

subplot(211)
loglog(avg.DTDZ(:),dtdzi(:),'d')
grid on
xvec=linspace(1e-4,1e0,100);
hold on
loglog(xvec,xvec,'k--')
grid on
axis tight
xlim([1e-3 1e0])
ylim([1e-3 1e0])
xlabel('T_z 1-m', 'fontsize', 16)
ylabel('T_z hi-res interp', 'fontsize', 16)

subplot(212)
loglog(avg.N2(:),n2i(:),'d')
grid on
xvec=linspace(1e-7,1e-2,100);
hold on
loglog(xvec,xvec,'k--')
grid on
axis tight
xlim([1e-7 1e-2])
ylim([1e-7 1e-2])
xlabel('N^2 1-m' , 'fontsize', 16)
ylabel('N^2 hi-res interp', 'fontsize', 16)

%% Plot histograms

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(211)
h1=histogram(real(log10(dtdz(:))),20,'Normalization','pdf','EdgeColor','none')
hold on
histogram(real(log10(avg.DTDZ(:))),h1.BinEdges,'Normalization','pdf','EdgeColor','none')
grid on
ylabel('pdf', 'fontsize', 16)
xlabel('log_{10}[T_z]', 'fontsize', 16)

subplot(212)
h1=histogram(real(log10(n2(:))),20,'Normalization','pdf','EdgeColor','none')
hold on
histogram(real(log10(avg.N2(:))),h1.BinEdges,'Normalization','pdf','EdgeColor','none')
grid on
ylabel('pdf', 'fontsize', 16)
xlabel('log_{10}[N^2]', 'fontsize', 16)

%%
figure(2);clf
h1=histogram(dtdz(:),'Normalization','pdf','EdgeColor','none')
hold on
histogram(avg.DTDZ(:),h1.BinEdges,'Normalization','pdf','EdgeColor','none')
grid on

%% Loop through all profiles to compare values

clear ; close all

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code

n2_raw=[];
n2_bin=[];
n2_raw_i=[];

dtdz_raw=[];
dtdz_bin=[];
dtdz_raw_i=[];

gam_raw_i=[];
gam_bin=[];

Pall=[];
warning off
for cnum=15:12:3100
    
    cnum
    
    clear cal avg cal2
    
    try
        
        % load binned chameleon profile ('avg') to compare with N2 and T_z computed from
        % linear fits over 1m bins
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_01366/mat/eq14_' sprintf('%04d',cnum) '.mat'])
        
        % load raw chameleon profile (not binned)
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal/EQ14_' sprintf('%04d',cnum) '.mat'])
        clear cal
        cal=cal2;
        clear cal2
        
        % sorted temperature
        [cal.T1sort,I] = sort(cal.T1,1,'descend');
        
        % smooth salinity a little bit
        cal.SAL=smooth(cal.SAL,20);
        
        % compute pot. dens.
        cal.sgth=sw_pden(cal.SAL,cal.T1,cal.P,0);
        
        % sort pot. dens
        [val,Irho] = sort(cal.sgth,1,'ascend');
        %n2=sw_bfrq(cal.SAL,cal.T1,cal.P,0);
        n2=sw_bfrq(cal.SAL(Irho),cal.T1(Irho),cal.P,0);
        n2=[n2(:) ; nan ];
        
        %dtdz=-diffs(cal.T1) ./ diffs(cal.P) ;
        dtdz=-diffs(cal.T1sort) ./ diffs(cal.P) ;
         
        % interpolate hi-res N2,T_z to 1m grid to compare
        clear IA n2i dtdzi
        [C,IA,IC]=unique(cal.P);
        n2i=interp1(cal.P(IA),n2(IA),avg.P);
        dtdzi=interp1(cal.P(IA),dtdz(IA),avg.P);
        
        n2_raw=[n2_raw(:) ; n2(:)];
        n2_bin=[n2_bin(:) ; avg.N2(:)];
        n2_raw_i=[n2_raw_i(:) ; n2i(:) ];
        
        dtdz_raw=[dtdz_raw(:) ; dtdz(:)];
        dtdz_bin=[dtdz_bin(:) ; avg.DTDZ(:)];
        dtdz_raw_i=[dtdz_raw_i(:) ; dtdzi(:) ];
        
        clear gam1 gam2
        gam1=ComputeGamma(n2i,dtdzi,avg.CHI,avg.EPSILON);
        gam2=ComputeGamma(avg.N2,avg.DTDZ,avg.CHI,avg.EPSILON);
        
        gam_raw_i=[gam_raw_i(:) ;  gam1(:) ] ;
        gam_bin = [gam_bin(:) ; gam2(:) ] ;
        
        Pall=[Pall(:) ; avg.P(:) ] ;
        
    end % try
    
end

warning on
%
dtdz_raw(find(dtdz_raw==Inf))=nan;
n2_raw(find(n2_raw==Inf))=nan;
n2_raw(find(n2_raw==-Inf))=nan;
%%

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(211)
hbin=histogram( dtdz_bin(:) , 80 , 'Normalization','pdf','Edgecolor','none');
hold on
hraw=histogram( dtdz_raw(:) , hbin.BinEdges, 'Normalization','pdf','Edgecolor','none');
legend([hbin hraw],'bin','raw','location','best')
grid on
xlim([-0.05 0.25])
xlabel('T_z', 'fontsize', 16)
ylabel('pdf', 'fontsize', 16)
freqline(nanmean(dtdz_bin(:)),'b--')
freqline(nanmean(dtdz_raw(:)),'r--')

subplot(212)
hraw=histogram( real(log10(n2_raw(:))) , 'Normalization','pdf','Edgecolor','none');
hold on
hbin=histogram( real(log10(n2_bin(:))), 'Normalization','pdf','Edgecolor','none');
legend([hbin hraw],'bin','raw','location','best')
grid on
xlim([-7 -1.5])
xlabel('log_{10}[N^2]', 'fontsize', 16)
ylabel('pdf', 'fontsize', 16)
freqline(nanmean( real(log10(n2_raw(:)))),'b')
freqline(nanmean( real(log10(n2_bin(:)))),'r')

% save plot
SetNotesFigDir
figname=['eq14_All_Tz_N2_raw_bin_hist']
print( fullfile( NotesFigDir , figname ) , '-dpng' )


%%

figure(1);clf

ax1=subplot(121);
boxplot(dtdz_bin(:))
ylim([-10 10])
grid on

ax2=subplot(122);
boxplot(dtdz_raw(:))
ylim([-10 10])
grid on

linkaxes([ax1 ax2],'y')

%% compare binned values to raw values interp to bins

figure(1);clf
agutwocolumn(1)
wysiwyg

to_plot=real(log10(dtdz_raw_i ./ dtdz_bin));
to_plot(find(to_plot==-Inf))=nan;
subplot(121)
histogram( to_plot , 'Normalization','pdf','Edgecolor','none')
xlim([-2 2])
grid on
xlabel('log_{10} [T_zi / T_zbin]')
freqline(nanmean(to_plot))

to_plot=real(log10(n2_raw_i ./ n2_bin));
to_plot(find(to_plot==-Inf))=nan;

subplot(122)
histogram( to_plot, 'Normalization','pdf','Edgecolor','none')
xlim([-3 3])
grid on
xlabel('log_{10} [N2_zi / N2_zbin]')
freqline(nanmean(to_plot))

rat_raw=n2_raw_i ./ (dtdz_raw_i.^2);
rat_raw(find(rat_raw==Inf))=nan;
rat_raw(find(rat_raw==-Inf))=nan;

rat_bin=n2_bin ./ (dtdz_bin.^2);
rat_bin(find(rat_bin==Inf))=nan;
rat_bin(find(rat_bin==-Inf))=nan;


figure(2);clf
hbin=histogram(real(log10(rat_bin)), 'Normalization','pdf','Edgecolor','none')
hold on
hraw=histogram(real(log10(rat_raw)),hbin.BinEdges, 'Normalization','pdf','Edgecolor','none')
xlim([-3.5 3.5])
legend([hbin hraw],'bin','raw')
xlabel('log_{10}[N^2/<T_{z}^{2}>]')
grid on
ylabel('pdf')
freqline(nanmean(real(log10(rat_bin))),'b--')
freqline(nanmean(real(log10(rat_raw))),'r--')

%%

rat1=avg.N2 ./ (avg.DTDZ.^2) ;
rat2=n2_raw ./ (dtdz_raw.^2) ;

figure(1);clf
histogram( real(log10(rat1(:))), 'Normalization','pdf','Edgecolor','none')
hold on
histogram( real(log10(rat2(:))), 'Normalization','pdf','Edgecolor','none')

%% compare gamma

ib=find(Pall<20);
gam_bin(ib)=nan;
gam_raw_i(ib)=nan;

figure(1);clf
hbin=histogram(real(log10(gam_bin(:))), 'Normalization','pdf','Edgecolor','none');
hold on
hrawi=histogram(real(log10(gam_raw_i(:))), 'Normalization','pdf','Edgecolor','none');
xlim([-4.5 2])
freqline(log10(0.2))
xlabel('log_{10}\Gamma')
grid on
legend([hbin hrawi],'bin','rawi','location','best')

%% compare the interpolated raw data in/out of patches

% load data on patches



%% compare the raw data in/out of patches

%%