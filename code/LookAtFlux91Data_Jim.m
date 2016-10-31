%~~~~~~~~~~~~~~~~~~~~~~~~~
%
% LookAtFlux91Data_Jim.m
%
% Look at flux 91 data that Jim shared. Will compute gamma etc. similar to
% EQ14 and compare results..
%
% I downloaded the data to my laptop at /Chipod/Flux91/
%
%-------------
% 10/17/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Flux91/flx91_patch_out_nov94.mat')

% compute gamma from this data

n2=flx91.bv.^2;% flx91.bv = N?
dtdz=flx91.dtdz;
chi=flx91.chi;
eps=flx91.eps;

gam =  n2 .* chi ./2 ./ eps ./ (dtdz.^2);

%%

h=Plot4hist(n2,dtdz,chi,eps)

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures'
print( fullfile( figdir,'flx91_patches_jim_hist'), '-dpng')

%%

figure(2);clf
agutwocolumn(0.5)
wysiwyg

histogram(gam(:),30)
freqline(nanmedian(gam))
title(['median=' num2str(roundx(nanmedian(gam),2))])
grid on

figdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/figures'
print( fullfile( figdir,'flx91_patches_jim_hist_gam'), '-dpng')

%%
