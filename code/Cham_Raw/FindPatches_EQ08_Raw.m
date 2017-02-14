%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% FindPatches_EQ08_Raw.m
%
% Modified from FindPatches_EQ14_Raw.m
%
% Find patches (overturns) in EQ08 chameleon profiles, using raw (not
% binned/averaged) data.
%
%
% The raw mat files for each chameleon cast are made w/
% ProcessEq14Cham_AP.m, which was modified from Sally's code so I could
% make files to apply chipod method to. See also ComputeChi_Chameleon_Eq14.m
% They are in the folder /ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal/
%
% Dependencies:
%   - compute_overturns_discrete_AP.m
%
%-----------------
% 2/14/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% Add all the paths we need from mixing software
mixpath='/Users/Andy/Cruises_Research/mixingsoftware/'
addpath(fullfile(mixpath,'seawater'))
addpath /Users/Andy/Standard-Mixing-Routines/ThorpeScales/

datdir='/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/casts'

% patch options
save_data = 1 ;         % save data at the end
patch_size_min = 0.25 ; % min patch size
usetemp   = 1 ;         % 1=use pot. temp, 0= use density

patch_data=[];

% loop through each cast
warning off
hb=waitbar(0,'working on profiles');

%% only do profiles that are done in chameleon processing

Flist=dir(fullfile(datdir,'*.mat'))
%%
for ic= 1:length(Flist)
    
    waitbar(ic/length(Flist),hb)
    
    try
        
        close all
        clear cal cal2 head
        
        clear cnum
        cnum=str2double(Flist(ic).name(5:8));%
        
        % Load the data for this cast
        load(fullfile(datdir,Flist(ic).name))
        
        cal=cal2; clear cal2
        
        clear s1 t1 p1 s t p lat
        %s=smooth( cal.S(1:end-1), 20 ); % (end-1) b/c last 2 values are same;
        s1=cal.S(1:end-1);
        t1=cal.T(1:end-1);
        p1=cal.P(1:end-1) ;
        
        % get only unqiue values of p (otherwise interp doesnt work in
        % sw_bfreq)
        [p,IA,IC] = unique(p1);
        s=s1(IA);
        t=t1(IA);
        
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
        
        clear OT
        OT=compute_overturns_discrete_AP(p,t,s,Params);
        
        clear pstarts pstops
        pstarts=OT.pstarts_each;
        pstops=OT.pstops_each;
        
        for i=1:length(pstarts)
            % don't keep patches shallower than 10m depth
            if pstarts(i)>10
                patch_data=[patch_data ; cnum pstarts(i) pstops(i)  ( pstops(i) - pstarts(i) ) OT.Otnsq_each(i) OT.Lt_each(i) ];
            end
        end
        
    end % try
    
end % cnum

delete(hb)
warning on

if save_data==1
    savedir = '/Users/Andy/Cruises_Research/ChiPod/Eq08/data/patches/'
    fname   = ['EQ08_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '.mat']
    save( fullfile( savedir,fname), 'patch_data')
end

%%

figure(1);clf
plot(t,p)
axis ij
grid on
hold on
plot(20,patch_data(:,2),'o')


%%
