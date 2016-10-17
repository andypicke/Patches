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

figure(1);clf
histogram(gam(:),30)

nanmedian(gam)

%%
