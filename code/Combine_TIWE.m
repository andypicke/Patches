%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Combine_TIWE.m
%
% Combine chameleon profiles from TIWE into a single structure for easier
% plotting and analysis etc.
%
% I'm using the processed mat files in /Tiwe91/mat_Greg_analysis/
%
%-------------
% 10/21/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% combine all profiles into a single structure for plotting etc

p_com=0:210; % interp all profiles to this common pressure vector

%addpath /Users/Andy/Cruises_Research/seawater_ver3_2/
datdir='/Users/Andy/Dropbox/ap_share_with_jn/date_from_jim/Tiwe91/mat_Greg_analysis/'
file_list=dir( fullfile( datdir,['tw91*' '*.mat']) );
%dat_sum=[]

tiwe=struct()
tiwe.t=nan*ones(length(p_com),length(file_list));
tiwe.s=nan*ones(length(p_com),length(file_list));
tiwe.n2=nan*ones(length(p_com),length(file_list));
tiwe.dtdz=nan*ones(length(p_com),length(file_list));
tiwe.chi=nan*ones(length(p_com),length(file_list));
tiwe.eps=nan*ones(length(p_com),length(file_list));

%warning off
for ip=1:length(file_list)
    clear avg
%    try
        load( fullfile('/Users/Andy/Dropbox/ap_share_with_jn/date_from_jim/Tiwe91/mat_Greg_analysis/',file_list(ip).name ) )
        
        % this is 1m binned data ; see if there are any overturns
        
        tiwe.t(:,ip)=interp1(avg.P,avg.T1,p_com);
        tiwe.s(:,ip)=interp1(avg.P,avg.S,p_com);
        tiwe.n2(:,ip)=interp1(avg.P,avg.N2,p_com);
        tiwe.dtdz(:,ip)=interp1(avg.P,avg.DTDZ,p_com);
        tiwe.chi(:,ip)=interp1(avg.P,avg.CHI,p_com);
        tiwe.eps(:,ip)=interp1(avg.P,avg.EPSILON,p_com);
        
    
end % ip

%
tiwe.P=p_com;
tiwe.ip=1:length(file_list);

% save combined structure

%%

figure(1);clf
%ezpc(tiwe.ip,tiwe.P,tiwe.t)
ezpc(tiwe.ip,tiwe.P,log10(tiwe.eps))
colorbar
caxis([-11 -5])
xlabel('profile index')
ylabel('P')

%%

%gam=tiwe.n2 .* tiwe.chi .

