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


%%
