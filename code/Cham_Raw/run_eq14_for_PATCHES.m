%~~~~~~~~~~~~~~~
%
% run_eq14_for_PATCHES.m
%
% Modified version of run_eq14.m . Only difference is that it calls
% average_data_PATCH_AP.m instead of average_data_gen1.m . This computes
% all data for specified patches, instead of the normal evenly-spaced bins.
% Data paths also modified from Sally's version.
%
% output saved to:
% '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/';
%
%--------
% 11/01/16 - A.Pickering
%~~~~~~~~~~~~~~~
%
% run_eq14_AP.m
%
% Modified from run_eq14.m on 06/28/16 by A.Pickering.
%
% Same as run_eq14, but for chi calc use (1) no response correction and (2)
% differnt fmax in order to directly compare w/ chipod method.
%
% ** note I tried to pass in fmax to chi_calc_AP but couldn't get it to
% work, so have to manually change it for each run :( .
%
% Changes:
%
% 1) average_data_gen1 > average_data_gen1_AP
%
% 2) In average_dat_gen1_AP:
% calc_chi > calc_chi_AP
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Adapting Sasha's processing from Dynamo11 for the eq14 cruise. The
% realtime processing doesn't go throuhg the cleaning and filtering steps
% that these codes do.
%
% By running this code it should process, combine, and clean all raw cham
% files. It calls functions: raw_load_cham2, cali_eq14, modify_header_eq14,
% fix_COND_eq14, average_data_gen1, set_filters_eq14, sum_eq14,
% tag_file_eq14, fix_summary_salt_eq14, make_mld, along with a lot of
% functions saved in mixingsoftware.
%
% all notes and comments added by sally
%
% Sally Warner, March 2015
%%

clear all

% option to test w/ 1m patches to make sure we get same results as original
% binned processing
run_test=0;

% load patch data (from FindPatches_EQ14_Raw.m)
patch_size_min=0.15  % min patch size
usetemp=1

savedir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc'
fname=['EQ14_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '.mat']
load(fullfile(savedir,fname))
clear savedir fname

addpath /Users/Andy/Cruises_Research/mixingsoftware/Chameleon2/Version2015/
addpath /Users/Andy/Cruises_Research/mixingsoftware/seawater/
addpath /Users/Andy/Cruises_Research/mixingsoftware/general/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/
addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/calibrate/
addpath /Users/Andy/Cruises_Research/mixingsoftware/calibrate/
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/mfiles/

% define important data and paths

cruise_id = 'eq14';
depth_max=250;

%path_cham = '~/GDrive/data/eq14/chameleon/';
path_cham='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/'
%path_raw = [path_cham 'raw/'];
path_raw='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/raw/'

% make directories for the single drop .mat files and for the summary file
% turn warnings off so they don't say that a directory already exists
warning off

if run_test==1
    path_save = ['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/test1m/'];
else
    path_save = ['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '/'];
end
ChkMkDir(path_save)
warning on

% define imporant variables
% define global variables which are called by many of the scripts
global data head cal q

% make structure q.script which will be used eventually to call the correct
% file by raw_load
q.script.pathname =  path_raw;
q.script.prefix = cruise_id;

% define the variables to be processed
% (note Pavan changed a number of the variable names, so they are different
% from the dynamo versions of the code)
% (note: neet to have T1, T2, SAL before epsilon and chi)
% (note: no need to list dT/dz or N2 here. They are calculated from 1m bin
% averaged data and therefore are not needed in this list)
series = {'fallspd','t1','t2','tp','cond','sal','theta','sigma',...
    'epsilon1','epsilon2',...
    'chi','az2','scat','ax_tilt','ay_tilt',...
    'varaz','varlt'};

% old versions from Sasha's code (first was NaNed out, second was used)
% series={'fallspd','t','tp','c','t2','mhc','s','theta','sigma','sigma_order','n2','epsilon1','epsilon2',...
%         'chi','az2','dtdz','drhodz','varaz','varlt','scat','ax_tilt','ay_tilt'};
% series={'fallspd','t','tp','c','t2','mhc','s','theta','sigma','epsilon1','epsilon2',...
%         'chi','az2','varaz','varlt','scat','ax_tilt','ay_tilt'};
warning off


% PROCESS EACH CAST

% manually define the cast numbers to be processed
% (note: in the realtime code, the "bad" flag is used for casts that are no
% good for whatever reason. Here, you're just supposed to have a matrix so
% the bad files are skipped over. In other words, "bad" isn't implemented.)
for cast = [4:12 14:46 48:87 374:519 550:597 599:904 906:909 911:1070 ...
        1075:1128 1130:1737 1739:2550 2552:2996 2998:3092];
    
    try
        
        disp(['processing cast number: ' num2str(cast)]);
        q.script.num=cast;
        q.series=series;
        temp1=q;
        
        % weird cludgy fix for redefining q for each cast (need to delete the
        % other variables that are saved in q for previous casts)
        clear global head data cal q
        global data head cal q
        q=temp1;
        
        %%%%%%%% LOAD %%%%%%%%
        % changing this to raw_load_cham2 to load the new chameleon files
        dummy = num2str(cast + 10000);
        load_file = [path_raw cruise_id '_' dummy(2) '.' dummy(3:5)];
        [data head]=raw_load_cham2(load_file);
        
        %%%%%%%% CALIBRATE %%%%%%%%
        % calibrate the raw voltages into useful data using the calibration
        % coefficients saved in the header OR new calibration coefficients
        % which are now defined in MODIFY_HEADER_**
        cali_eq14;
        
        % The variable "bad" comes from cali_realtime. If NOT bad, continue
        % with code, otherwise do not save a .mat file for this cast.
        if bad ~= 1
            
            %%%%%%%% 1m AVERAGE AND CALCULATE EPSILON AND CHI %%%%%%%%
            % Now, average the data into 1m bins and calculate epsilon and chi
            % from the 1m wavenumber spectra
            % Need to create SET_FILTERS_** (called by average_data_gen1)
            % which sets how the shear spectra should be filtered.
            nfft=128;
            warning off
            
            % *** replace with average_data_PATCH_AP.m, need to load patch data also
            % and specify in inputs ***
            nfft=128;
            
            clear igp pstarts pstops this_patch
            
            if run_test==1
                pstarts=1:200;
                pstops=2:201;
            else
                igp=find(patch_data(:,1)==cast);
                this_patch=patch_data(igp,:);
                pstarts=patch_data(igp,2);
                pstops=patch_data(igp,3);
            end
            
            avg=average_data_PATCH_AP(q.series,nfft,pstarts,pstops);
            
            if length(avg.P)~=length(pstarts)
                disp('uh oh')
            end
            
            %~~~ also add N2 from the overturns code for each patch, to
            %compare with the N2 I calc in average_data_PATCH_AP
            avg.n2_OT=nan*ones(length(avg.P),1);
            if run_test~=1
                for ip=1:length(avg.P)
                    clear ig
                    ig=find(avg.P(ip)>pstarts & avg.P(ip)<pstops) ;
                    avg.n2_OT(ip) = this_patch(ig,5);
                end
                
            end
            
            %%%%%%%% calculate N2, dSAL/dz, dSAL/dz_rhoorder, dT/dz, dT/dz_rhoorder %%%%%%%%
            % done by using polyfit over 1m sections of density (or reordered
            % density) to get the slope (drho/dz, dsal/dz, dT/dz, etc)
            % (note: the following all used positive pressure and no negative sign
            % before g/rho*drhodz. So N2 had the right sign, but all of the rest
            % had the wrong sign. Added appropriate negative signs, because z
            % should be positive upward. [SJW, March 2015])
            %         g=9.81;
            %         rhoav=nanmean(avg.SIGMA)+1000;
            %         for ii=1:length(avg.P)
            %             % find data points (in cal) for each 1m section
            %             in=find(cal.P>=avg.P(ii)-0.5 & cal.P<avg.P(ii)+0.5);
            %
            %             % linear fit to get slope: pp(1) = drho/dz, then calc Nsq (N2)
            %             pp=polyfit(-cal.P(in),cal.SIGMA_ORDER(in),1);
            %             avg.N2(ii)=-g/rhoav*pp(1);
            %
            %             % linear fit to get slope of salinity: pp(1) = dSAL/dz
            %             pp=polyfit(-cal.P(in),cal.SAL(in),1);
            %             avg.DSALDZ(ii)=pp(1);
            %
            %             % linear fit to salinity (that has been reordered by density)
            %             pp=polyfit(-cal.P(in),cal.SAL_RHOORDER(in),1);
            %             avg.DSALDZ_RHOORDER(ii)=pp(1);
            %
            %             % linear fit to get slope of theta: pp(1) = dT/dz
            %             pp=polyfit(-cal.P(in),cal.THETA(in),1);
            %             avg.DTDZ(ii)=pp(1);
            %
            %             % linear fit to theta (that has been reordered by density)
            %             pp=polyfit(-cal.P(in),cal.THETA_RHOORDER(in),1);
            %             avg.DTDZ_RHOORDER(ii)=pp(1);
            %         end
            %
            %%%%%%%% NAN OUT GLITCHES %%%%%%%%
            % VARAZ is the variance of the vertical acceleration. When glitches
            % occur as chameleon is falling, VARAZ will be MUCH LARGER than it
            % would be while in freefall
            % find times when VARAZ is > 1e-2 and NaN out those bad values
            idaz=find(avg.VARAZ>1.e-02);
            avg.EPSILON1(idaz)=NaN;
            avg.EPSILON2(idaz)=NaN;
            avg.EPSILON(idaz)=NaN;
            
            % find times when the squared vertical acceletion is greater than a
            % threshold and nan out bad values
            temp2=find(log10(avg.AZ2)>-5);
            avg.EPSILON1(temp2)=NaN;
            avg.EPSILON2(temp2)=NaN;
            avg.EPSILON(temp2)=NaN;
            
            % find values of EPSILON that are too large and NaN them out
            bad1=find(avg.EPSILON1>1e-4);
            avg.EPSILON1(bad1)=NaN;
            bad2=find(avg.EPSILON2>1e-4);
            avg.EPSILON2(bad2)=NaN;
            bad=find(avg.EPSILON>1e-4);
            avg.EPSILON(bad)=NaN;
            
            % Flag out any data collected shallower than 5m depth
            idsur=find(avg.P<=5);
            avg.EPSILON1(idsur)=NaN;
            avg.EPSILON2(idsur)=NaN;
            avg.EPSILON(idsur)=NaN;
            
            % if there are specific casts that need to use epsilon1 rather than
            % epsilon 2 or vice versa, state these now.
            %     if cast<=77 || (cast>=1866 && cast<=2232) || (cast>=2917 && cast<=2941)
            %         avg.EPSILON=avg.EPSILON1;
            %     elseif cast>=2913 && cast<=2916
            %         avg.EPSILON=avg.EPSILON2;
            %     end
            
            warning backtrace
            
            %%%%%%%% DYNAMIC HEIGHT %%%%%%%%
            head=calc_dynamic_z(avg,head);
            
            %%%%%%%% ADD MAXIMUM DEPTH TO HEAD %%%%%%%%
            head.p_max=max(cal.P);
            
            %%%%%%%% SAVE INDIVIDUAL CASTS %%%%%%%%
            temp=num2str(q.script.num+10000);
            fn=[q.script.prefix '_' temp(2:5)];
            
            avg.readme = strvcat('Processed Chameleon Data from individual cast',...
                ['cruise = ' q.script.prefix],...
                'Data has been averaged into 1m bins',...
                'Epsilon calculated from wavenumber spectra of 1m bins',...
                'Epsilon spectra are filtered to removed narrowband strumming',...
                'Further flags on epsilon and other variables to be implemented',...
                '   in the sum file',...
                'See header for more details on this file',...
                ['saved on ' datestr(now) ' by run_' q.script.prefix '.m']);
            
            
            save(fullfile(path_save,fn),'avg','head')
            
            
            % if the data is not good, (i.e. if bad == 1), do not save anything
        else
            
        end
        
    catch
    end % try
    
end % cast

%%
%%%%%%%% COMBINE ALL CASTS AND SAVE SUM FILE %%%%%%%%
%sum_eq14(path_save,path_sum,depth_max,cruise_id);


