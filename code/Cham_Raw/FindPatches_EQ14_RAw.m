%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FindPatches_EQ14_Raw.m
%
% Find patches (overturns) in EQ14 chameleon profiles, using raw (not
% binned/averaged) data.
%
% Ideally i'd like to do this for tiwe data and compare to what Bill got
% for his patches. But i'm doing it on this better-documented data set
% first to figure it out. I'll compare the patches from this to what I get
% using the 1m binned/avg data.
%
% The raw mat files for each chameleon cast are made w/
% ProcessEq14Cham_AP.m, which was modified from Sally's code so I could
% make files to apply chipod method to. See also ComputeChi_Chameleon_Eq14.m
%
% Dependencies:
%   - compute_overturns_discrete_AP.m
%
%-----------------
% 10/27/16 - A.Pickering - andypicke@gmail.com
% 11/01/16 - AP - Use specific cast #s so we can match up with other data
% later more easily
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% Add all the paths we need from mixing software
mixpath='/Users/Andy/Cruises_Research/mixingsoftware/'
addpath(fullfile(mixpath,'seawater'))

datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal'

% patch options
save_data = 1 ;        % save data at the end
patch_size_min = 0.25 ;% min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density

patch_data=[];

% loop through each cast
warning off
hb=waitbar(0,'working on profiles');

% only do profiles that are done in chameleon processing
cnums_to_do=[4:12 14:46 48:87 374:519 550:597 599:904 906:909 911:1070 ...
    1075:1128 1130:1737 1739:2550 2552:2996 2998:3092];

for cnum= cnums_to_do;
    
    waitbar(cnum/length(cnums_to_do),hb)
    
    try
        
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
        
        clear Params
        Params.lat=lat;
        Params.plotit=0;
        Params.sigma=1e-5;
        Params.runlmin=0;
        Params.minotsize=patch_size_min;
        Params.usetemp=usetemp;
        addpath /Users/Andy/Standard-Mixing-Routines/ThorpeScales/
        clear OT
        OT=compute_overturns_discrete_AP(p,t,s,Params);
        
        pstarts=OT.pstarts_each;
        pstops=OT.pstops_each;
        
        for i=1:length(pstarts)
            patch_data=[patch_data ; cnum pstarts(i) pstops(i)  ( pstops(i) - pstarts(i) ) OT.Otnsq_each(i) OT.Lt_each(i) ];
        end
        
    end % try
    
end % cnum

delete(hb)
warning on


new_patch_data=patch_data;

if save_data==1
    savedir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/'
    fname=['EQ14_raw_patches_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '.mat']
    save( fullfile( savedir,fname), 'new_patch_data')
end

%%
