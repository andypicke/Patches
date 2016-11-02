%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute_Patch_Data_EQ14_Raw.m
%
% See also FindPatches_EQ14_Raw.m ; in this code I also try to compute
% n2,dtdz, chi, and eps over each patch.
%
%
% 10/31/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% Make a list of all the chameleon casts we have (made w/ ProcessEq14Cham_AP.m)
datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal'
Flist=dir(fullfile(datdir,'*EQ14*.mat'))

% load patch data (from FindPatches_EQ14_Raw.m )
% data contains: [ icast ; pstart ; pstop ; patch size ]
patch_size_min=1  % min patch size
join_patches=0     % join nearby patches
patch_sep_min=0.15 %
savedir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/'
fname=['EQ14_raw_patches_minOT_' num2str(patch_size_min) '_join_' num2str(join_patches) '_sep_' num2str(100*patch_sep_min) '.mat']
load(fullfile(savedir,fname))

icast=1%:length(Flist)

%    waitbar(icast/length(Flist),hb)
%    disp(icast)
close all
clear cal head
clear tpspec kspec kkspec fspec kks ks

% Load the data for this cast
load(fullfile(datdir,Flist(icast).name))

clear cal
cal=cal2;

clear s t p lat
s=cal.SAL;
t=cal.T1;
p=cal.P;

clear idot lat1 lat2
idot=strfind(head.lat.start,'.');
lat1=str2num(head.lat.start(1:idot-3));
lat2=str2num(head.lat.start(idot-2:end))/60;
lat=nanmean([lat1 lat2]);

sgth=sw_pden(s,t,p,0);

% now LOOP through each patch
clear ig patch_prof
ig=find(new_patch_data(:,1)==icast)
patch_prof=new_patch_data(ig,:);

if length(ig)>1
    
    % make empty arrays
    n2=nan*ones(length(ig),1);
    dtdz=nan*ones(length(ig),1);
    chi=nan*ones(length(ig),1);
    eps=nan*ones(length(ig),1);
            
    for ii=1%:length(ig)
        
        % find indices of data in this patch
        clear iz dz del_rho
        iz=isin(cal.P,[patch_prof(ii,2) patch_prof(ii,3) ] );
        
        dz=nanmax(cal.P(iz)) - nanmin(cal.P(iz)) ;
        dT=nanmax(cal.T1(iz)) - nanmin(cal.T1(iz));
        del_rho= nanmax(sgth(iz)) - nanmin(sgth(iz)) ;
        
        % compute N2 and dT/dz (just do for entire patch)
        n2(ii)= 9.81 / nanmean(sgth(iz)) * del_rho / dz ;
        dtdz(ii) = dT/dz ;
        
        % compute chi
        
        
        % compute epsilon
        
        %
    end % ii (each patch)
    
end % if we have patches for this profile

%%