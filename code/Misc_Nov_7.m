%
%
% Misc_Nov_7.m
%
%
%

clear ; close all

% Add all the paths we need from mixing software
mixpath='/Users/Andy/Cruises_Research/mixingsoftware/'
addpath(fullfile(mixpath,'seawater'))

% for IdentifyPatches.m
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/

datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal'

% patch options
save_data=1
patch_size_min=0.5  % min patch size
join_patches=0    % join nearby patches
patch_sep_min=0.5 %

patch_data=[];
% loop through each cast
warning off
%hb=waitbar(0,'working on profiles');

for cnum=800%:3100%length(Flist)
    
    %   waitbar(cnum/3100,hb)
    
    %    try
    
    close all
    clear cal head
    clear tpspec kspec kkspec fspec kks ks
    
    % Load the data for this cast
    load(fullfile(datdir,['eq14_' sprintf('%04d',cnum) '.mat']))
    
    clear cal
    cal=cal2;
    
    clear s t p lat
    s=cal.SAL(1:end-1); % (end-1) b/c last 2 values are same
    t=cal.T1(1:end-1);
    p=cal.P(1:end-1);
    
    clear idot lat1 lat2
    idot=strfind(head.lat.start,'.');
    lat1=str2num(head.lat.start(1:idot-3));
    lat2=str2num(head.lat.start(idot-2:end))/60;
    lat=nanmean([lat1 lat2]);
    
    Params.lat=lat
    Params.plotit=1
    Params.sigma=1e-4
    Params.runlmin=0
    Params.minotsize=0.5;
    Params.usetemp=1
    addpath /Users/Andy/Standard-Mixing-Routines/ThorpeScales/
    OT=compute_overturns_discrete_AP(p,t,s,Params)
    
end

%%