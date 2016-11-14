function [pstarts pstops]=IdentifyPatches(s,t,p,lat)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% [pstarts pstops]=IdentifyPatches(s,t,p,lat)
%
% Identify overturn patches from profile of potential density.
% Modified from compute_overturns_discrete_AP.m **
%
% INPUT
% s - profile of salinity
% t - profile of temperature
% p - profile of pressure
% lat - latitude of profile
%
% OUTPUT
%
% pstarts - Pressures of top of patches
% pstops  - Pressures of bottom of patches
%
% DEPENDENCIES
% sw_pden.m
% sw_ptmp.m
% sw_bfrq.m
%
%-------------------------
% 10/7/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear pden ptmp
refd=0;
plots=0;

% compute potential density and temperature
pden = sw_pden(s(:),t(:),p(:),refd);
ptmp = sw_ptmp(s(:),t(:),p(:),refd);

%
if plots==1
    figure(1);clf
    subplot(121)
    plot(t,p);axis ij;grid on
    xlabel('t')
    subplot(122)
    plot(s,p);axis ij;grid on
    xlabel('s')
end

% choose temperature or density to use
%if usetemp
%    V=ptmp;
%else
V=pden;
%end

clear xx isort n2 p0 ig
% sort density profile
[xx,isort]=sort(pden);

% compute N^2 with sorted (stable) density profile
[n2,q,p_ave] = sw_bfrq(s(isort),t(isort),p,lat);

p0=p(:);  % full depth

% find good (not NaN density values)
ig=find(~isnan(pden));

%    if numel(ig)>1 % only do if we have good data

% make vectors with only good values
clear pg
pg=p(ig);
ptmp=ptmp(ig);
pden=pden(ig);
V=V(ig);

% make sure presure is a column vector
pg=pg(:);

% distinguish between up/down casts?
sig = sign(nanmedian(diff(V)));

% sort profile
[Vsort,ind]=sort(sig*V);
% tsort=sig*V; % AP 13 Feb
psort = pg(ind);
dz = pg-psort;

csdz = cumsum(-dz);
thresh = 0.0000001;

% Find start and stop indices for overturn region(s)
start = find(csdz(1:end-1)<thresh & csdz(2:end)>=thresh)+1;
if dz(1)<0
    start = [1;start];
end;
stops = find(csdz(1:end-1)>=thresh & csdz(2:end)<thresh)+1;

%
if plots==1
    figure(1);clf
    
    ax1=subplot(131);
    plot(n2,p_ave)
    xlim([0 nanmax(n2)])
    xlabel('N^2')
    axis ij
    grid on
    ylabel('p [db]')
    
    ax2=subplot(132);
    plot(sig*V,pg,Vsort,pg)
    hold on
    plot(V(start),pg(start),'bo')
    plot(V(stops),pg(stops),'rd')
    axis ij
    xlabel('sgth')
    grid on
    ytloff
    
    ax3=subplot(133);
    plot(dz,pg,'o-')
    xlabel('dz')
    axis ij
    grid on
    ytloff
    
    linkaxes([ax1 ax2 ax3],'y')
end

%% require overturn size > min resolvable overturn based on N2 and density noise level

clear ipass
ipass=[];
sigma=5e-4;

for ii=1:length(start)
    clear ind drho n2avg indp
    ind=start(ii):stops(ii);
    indp=find( p_ave>min(pg(ind)) & p_ave<max(pg(ind)) );
    drho=(max(pden(ind))-min(pden(ind)));
    n2avg=nanmean(n2(indp));
    delz=abs(max(pg(ind))-min(pg(ind)));
    
    if   delz>(2*9.8/n2avg*sigma/1027)
        ipass=[ipass ; ii];
    else
    end
    
end

pstarts=p0(start(ipass));
pstops=p0(stops(ipass));


%%