%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Combine_TIWE.m
%
% Combine chameleon profiles from TIWE into a single structure for easier
% plotting and analysis etc.
%
% I'm using the processed mat files in /Tiwe91/mat_Greg_analysis/
%
% All the profiles are interpolated onto a common depth vector with 1 m
% spacing
%
%----------------------
% 10/21/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% combine all profiles into a single structure for plotting etc

p_com=0:210; % interp all profiles to this common pressure vector

%addpath /Users/Andy/Cruises_Research/seawater_ver3_2/
datdir='/Users/Andy/Dropbox/ap_share_with_jn/date_from_jim/Tiwe91/mat_Greg_analysis/'
file_list=dir( fullfile( datdir,['tw91*' '*.mat']) );

tiwe=struct()
tiwe.t=nan*ones(length(p_com),length(file_list));
tiwe.s=nan*ones(length(p_com),length(file_list));
tiwe.N2=nan*ones(length(p_com),length(file_list));
tiwe.DTDZ=nan*ones(length(p_com),length(file_list));
tiwe.CHI=nan*ones(length(p_com),length(file_list));
tiwe.EPSILON=nan*ones(length(p_com),length(file_list));

%warning off
for ip=1:length(file_list)
    clear avg
%    try
        load( fullfile('/Users/Andy/Dropbox/ap_share_with_jn/date_from_jim/Tiwe91/mat_Greg_analysis/',file_list(ip).name ) )
        
        % this is 1m binned data ; see if there are any overturns        
        tiwe.t(:,ip)=interp1(avg.P,avg.T1,p_com);
        tiwe.s(:,ip)=interp1(avg.P,avg.S,p_com);
        tiwe.N2(:,ip)=interp1(avg.P,avg.N2,p_com);
        tiwe.DTDZ(:,ip)=interp1(avg.P,avg.DTDZ,p_com);
        tiwe.CHI(:,ip)=interp1(avg.P,avg.CHI,p_com);
        tiwe.EPSILON(:,ip)=interp1(avg.P,avg.EPSILON,p_com);        
    
end % ip

%%
tiwe.P=p_com;
tiwe.ip=1:length(file_list);
%
tiwe.MakeInfo=['Made ' datestr(now) ' w/ Combine_TIWE.m, source is /Tiwe91/mat_Greg_analysis/']
%
% save combined structure
save(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches','tiwe_comb_AP.mat'),'tiwe')
%%

clear ; close all

load(fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches','tiwe_comb_AP.mat') )

%% Plot temperature and salinity

figure(1);clf
agutwocolumn(0.8)
wysiwyg

subplot(211)
ezpc(tiwe.ip,tiwe.P,tiwe.t)
colorbar
title('temperature')
xlabel('profile index')
ylabel('P')

subplot(212)
ezpc(tiwe.ip,tiwe.P,tiwe.s)
colorbar
caxis([32 36])
title('salinity')
xlabel('profile index')
ylabel('P')

%%

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print( fullfile( figdir,'tiwe_comb_t_ts'), '-dpng')

%% Plot N2,dtdz,chi,eps

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(411)
ezpc(tiwe.ip,tiwe.P,real(log10(tiwe.N2)))
colorbar
%caxis([-11 -5])
xlabel('profile index')
ylabel('P')
title('N^2')

cmap=colormap;
colormap([0.5*[1 1 1] ; cmap ])

subplot(412)
ezpc(tiwe.ip,tiwe.P,real(log10(tiwe.DTDZ)))
colorbar
%caxis([-11 -5])
xlabel('profile index')
ylabel('P')
title('dT/dz')

subplot(413)
ezpc(tiwe.ip,tiwe.P,log10(tiwe.CHI))
colorbar
caxis([-11 -5])
xlabel('profile index')
ylabel('P')
title('\chi')

subplot(414)
ezpc(tiwe.ip,tiwe.P,log10(tiwe.EPSILON))
colorbar
caxis([-11 -5])
xlabel('profile index')
ylabel('P')
title('\epsilon')

%%

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print( fullfile( figdir,'tiwe_comb_n2_dtdz_chi_eps'), '-dpng')


%%
% compare to EQ distributions

dir_base='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed';
clear cham
load( fullfile( dir_base, '/Cstar=0_01366/sum/eq14_sum.mat') )

n2=tiwe.N2;
dtdz=tiwe.DTDZ;
chi=tiwe.CHI;
eps=tiwe.EPSILON;

Nm='pdf';

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
h1=histogram(real(log10(n2(:))),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(real(log10(cham.N2(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}N^2')
grid on


subplot(222)
h1=histogram(real(log10(dtdz(:))),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(real(log10(cham.DTDZ(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
%h3=histogram(real(log10(dtdz_2(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}dT/dz')
grid on

subplot(223)
h1=histogram(log10(chi(:)),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(real(log10(cham.CHI(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}\chi')
grid on
freqline(nanmedian(log10(chi(:))),'b')
freqline(nanmedian(log10(cham.CHI(:))),'r')

subplot(224)
h1=histogram(log10(eps(:)),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(real(log10(cham.EPSILON(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}\epsilon')
grid on
legend([h1 h2],'TIWE','EQ14')

%%

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print( fullfile( figdir,'tiwe_eq14_hist_compare'), '-dpng')

%% get rid of those weird dtdz values

% try removing the weird spike in values near 1e-4
%dtdz_2=dtdz;
ib2=find(log10(dtdz)<-3.7);
%dtdz_2(ib2)=nan;

n2_2=n2;n2_2(ib2)=nan;
dtdz_2=dtdz;dtdz_2(ib2)=nan;
chi_2=chi;chi_2(ib2)=nan;
eps_2=eps;eps_2(ib2)=nan;

% plot distributions again
Nm='pdf';

figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
h1=histogram(real(log10(n2_2(:))),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(real(log10(cham.N2(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}N^2')
grid on


subplot(222)
h1=histogram(real(log10(dtdz_2(:))),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(real(log10(cham.DTDZ(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
%h3=histogram(real(log10(dtdz_2(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}dT/dz')
grid on

subplot(223)
h1=histogram(log10(chi_2(:)),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(real(log10(cham.CHI(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}\chi')
grid on
freqline(nanmedian(log10(chi_2(:))),'b')
freqline(nanmedian(log10(cham.CHI(:))),'r')

subplot(224)
h1=histogram(real(log10(eps_2(:))),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(real(log10(cham.EPSILON(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}\epsilon')
grid on
legend([h1 h2],'TIWE','EQ14')

%%

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print( fullfile( figdir,'tiwe_eq14_hist_compare_fixdtdz'), '-dpng')


%% Compute gamma

%gam_tiwe=n2 .* chi ./2 ./ eps ./ (dtdz.^2);
gam_tiwe=n2 .* chi ./2 ./ eps ./ (dtdz_2.^2);
gam_eq14=cham.N2 .* cham.CHI ./2 ./cham.EPSILON ./(cham.DTDZ.^2);

% plot distributions
figure(1);clf
agutwocolumn(0.6)
wysiwyg
h1=histogram(real(log10(gam_tiwe(:))),'edgecolor','none','Normalization','pdf');
hold on
h2=histogram(log10(gam_eq14(:)),h1.BinEdges,'edgecolor','none','Normalization','pdf');
%hold on
%h3=histogram(real(log10(gam_tiwe_2(:))),'edgecolor','none','Normalization','pdf');
grid on
xlim([-5 4])
legend([h1 h2],'tiwe','eq14')
xlabel('log_{10}\Gamma','fontsize',16)
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')

%%

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures/'
print( fullfile( figdir,'tiwe_eq14_gamma_compare'), '-dpng')

%%
