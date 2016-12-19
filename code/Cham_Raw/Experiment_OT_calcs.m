%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Experiment_OT_calcs.m
%
% Experimenting w/ how different parameters affect OT calculations
%
%---------------
% 12/14/16 - AP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% Add all the paths we need from mixing software
mixpath='/Users/Andy/Cruises_Research/mixingsoftware/'
addpath(fullfile(mixpath,'seawater'))


datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal'

% patch options
save_data=1         % save data at the end
patch_size_min=0.15  % min patch size
usetemp=1           % 1=use pot. temp, 0= use density

%patch_data=[];

% only do profiles that are done in chameleon processing
cnums_to_do=[4:12 14:46 48:87 374:519 550:597 599:904 906:909 911:1070 ...
    1075:1128 1130:1737 1739:2550 2552:2996 2998:3092];
%for cnum=1:3100%
%for cnum= cnums_to_do;
cnum=801


%    try

close all
clear cal head
clear tpspec kspec kkspec fspec kks ks

% Load the data for this cast
load(fullfile(datdir,['eq14_' sprintf('%04d',cnum) '.mat']))

clear cal
cal=cal2;

clear s t p lat
s=cal.SAL(1:end-1); % (end-1) b/c last 2 values are same;
t=cal.T1(1:end-1);
p=cal.P(1:end-1);

clear idot lat1 lat2
idot=strfind(head.lat.start,'.');
lat1=str2num(head.lat.start(1:idot-3));
lat2=str2num(head.lat.start(idot-2:end))/60;
lat=nanmean([lat1 lat2]);

% find overturns
%        clear pstarts pstops
%       [pstarts pstops]=IdentifyPatches(s,t,p,lat);

clear Params
Params.lat=lat;
Params.plotit=1;
Params.sigma=1e-7;
Params.runlmin=0;
Params.minotsize=0.15 %patch_size_min;
Params.usetemp=usetemp;
addpath /Users/Andy/Standard-Mixing-Routines/ThorpeScales/
clear OT
OT=compute_overturns_discrete_AP(p,t,s,Params);

% end

%%