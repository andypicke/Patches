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
%
%
%-----------------
% 10/27/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% Add all the paths we need from mixing software
mixpath='/Users/Andy/Cruises_Research/mixingsoftware/'
addpath(fullfile(mixpath,'seawater'))

% Make a list of all the chameleon casts we have (made w/ ProcessEq14Cham_AP.m)
datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal'
Flist=dir(fullfile(datdir,'*EQ14*.mat'))

% options
save_data=1
patch_size_min=0.15   % min patch size
join_patches=0     % join nearby patches
patch_sep_min=0.15 %

patch_data=[];
% loop through each cast
warning off
hb=waitbar(0,'working on profiles');
for icast=1:length(Flist)
    waitbar(icast/length(Flist),hb)
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
    
    % find overturns
    clear pstarts pstops
    [pstarts pstops]=IdentifyPatches(s,t,p,lat);
    
    % should also join patches separated by less than 15cm?
    
    for i=1:length(pstarts)
        
        if ( pstops(i) - pstarts(i) ) > patch_size_min % reject patches thinner than 15cm
            patch_data=[patch_data ; icast pstarts(i) pstops(i)  ( pstops(i) - pstarts(i) ) ];
        else
        end
        
    end
    
end

delete(hb)
warning on

if join_patches==1
    % Join patches if they are separated by less than patch_sep_min
    
    prof_nums=unique(patch_data(:,1));
    new_patch_data=[]
    for ip=1:length(prof_nums)
        
        clear ig patch_small new_patch
        ig=find(patch_data(:,1)==prof_nums(ip));
        
        if length(ig)>1
            patch_small=patch_data(ig,1:3);
            ii=1;
            while ii<size(patch_small,1)
                clear new_patch
                % check if separated by <15cm
                if ( patch_small(ii+1,2) - patch_small(ii,3) ) < patch_sep_min
                    % join patch
                    new_patch=[prof_nums(ip) patch_small(ii,2) patch_small(ii+1,3)];
                    patch_small=[patch_small(1:ii-1,:) ; new_patch ; patch_small(ii+2:end,:)];
                else
                end
                ii=ii+1;
                
            end % patches for 1 profile
            
            % add patch size column
            patch_small=[patch_small patch_small(:,3)-patch_small(:,2)];
            
            %
            new_patch_data=[new_patch_data ; patch_small];
        else
            new_patch_data=[new_patch_data ; patch_data(ig,:)];
        end % have more than 1 patch
        
    end %ip (which profile
 
else % don't join patches
    new_patch_data=patch_data;
end % if join_patches==1

% save data
if save_data==1
savedir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/'
fname=['EQ14_raw_patches_minOT_' num2str(patch_size_min) '_join_' num2str(join_patches) '_sep_' num2str(100*patch_sep_min) '.mat']
save( fullfile( savedir,fname), 'new_patch_data')
end
%%

figure(1);clf
h1=histogram(new_patch_data(:,4));
freqline(nanmedian(h1.Data))
xlim([0 20])
grid on
xlabel('patch size')

%%

figure(2);clf
h1=histogram(new_patch_data(:,4));
freqline(nanmedian(h1.Data))
xlim([0 10])

%% plot patch locations over pcolor and compare to patches from 1m binned data

figure(1);clf
plot(new_patch_data(:,1),new_patch_data(:,2),'go')
hold on
plot(new_patch_data(:,1),new_patch_data(:,3),'ro')
axis ij

%%
figure(1);clf
plot(patch_data(:,1),patch_data(:,2),'go')
hold on
plot(patch_data(:,1),patch_data(:,3),'ro')
axis ij

%%

