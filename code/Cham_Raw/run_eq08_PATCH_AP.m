%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% run_eq08_PATCH_AP.m
%
% Modified version to compute data only for patches.
%
% Calls average_data_PATCH_AP.m instead of average_data_gen1.m
%
% 11/02/16 - A.Pickering
%
%---------------------
%
%
% Script to run processing of EQ08 Chameleon data.
%
% Modified from run_eq08.m, by Sasha.
%
% Dependencies:
% tag_file_eq08.m
% raw_load.m
% cali_eq08.m
%
%------------------
% 03/15/16 - A.Pickering - apickering@coas.oregonstate.edu
% 03/18/16 - AP - Modifying to save N2 and dTdz also for use in chi calcs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear all

addpath /Users/Andy/Cruises_Research/mixingsoftware/marlcham/calibrate/

% folder with raw Chameleon data
path_raw='/Volumes/SP PHD U3/NonBackup/EQ08/raw/'

% path to save processed files to
path_save='/Users/Andy/Cruises_Research/ChiPod/EQ08/Data/cham_proc/casts/'

% check if path_save exits, make new directory if not
ChkMkDir(path_save)

% run script to tag bad casts and depth ranges before processing
tag_file_eq08

global data head cal q
q.script.pathname =  path_raw;
q.script.prefix = 'eq08';
series={'fallspd','t','tp','c','mht','s','theta','sigma','sigma_order','n2','epsilon1','epsilon2',...
    'chi','az2','dtdz','drhodz','varaz','varlt','varlt_theta','scat','ax_tilt','ay_tilt'};
warning off
ic=0;
hb=waitbar(0)
% restart at 992
 for cast=[16:78,80:543,545:547,549,551:685,687:719,721:789,791,...
        793:798,800,802:988,990:1083,1085:1189,1191:1199,1201:1414,...
        1416:1420,1422:1500,1502,1504:1624,1631:1861,1863:1978,...
        1980:2124,2126:2140,2142:2298,2300:2668]
%     
    
    try
        %1524
        
        % bad, nonfixable files:
        % 79,792,799,801,1415,1421,1501,1503,1625,1862,2141,...
        %
        % up or surface profiles:
        % 544,548,550,686,720,989,1084,1190,1200,1979,2125,2299
        ic=ic+1;
        
        waitbar(ic/2624,hb)
        
        disp(['working on cast ' num2str(cast)]);
        q.script.num=cast;
        q.series=series;
        temp1=q;
        clear global head data cal q
        global data head cal q
        q=temp1;
        
        % load raw Chameleon data
        [data head]=raw_load(q);
        
        % apply calibrations to data
        cali_eq08;
        warning off
        
        % * 3/15/16 - AP - skip this part, just want calibrated T' *
        %     average calibrated data into 1m bins
        %     %     nfft=128;
        %     nfft=256;
        %     avg=average_data_gen1(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
        %     %            avg=average_data(q.series,'binsize',1,'nfft',nfft,'whole_bins',1);
        %     %
        %     % remove glitches
        %     % flag AZ vibrations
        %     idaz=find(avg.VARAZ>1.e-02);
        %     avg.EPSILON1(idaz)=NaN; avg.EPSILON2(idaz)=NaN; avg.EPSILON(idaz)=NaN;
        %     %     avg.WD(idaz)=NaN; avg.WD2(idaz)=NaN;
        %     temp2=find(log10(avg.AZ2)>-4.5);
        %     avg.EPSILON1(temp2)=NaN; avg.EPSILON2(temp2)=NaN; avg.EPSILON(temp2)=NaN;
        %     bad1=find(avg.EPSILON1>1e-4);
        %     avg.EPSILON1(bad1)=NaN;
        %     bad2=find(avg.EPSILON2>1e-4);
        %     avg.EPSILON2(bad2)=NaN;
        %     bad=find(avg.EPSILON>1e-4);
        %     avg.EPSILON(bad)=NaN;
        %     %         avg.WD(temp2)=NaN; avg.WD2(temp2)=NaN;
        %     % flag surface
        %     idsur=find(avg.P<=5);
        %     avg.EPSILON1(idsur)=NaN; avg.EPSILON2(idsur)=NaN; avg.EPSILON(idsur)=NaN;
        %     if avg.EPSILON(end)>20*avg.EPSILON(end-1);
        %         avg.EPSILON(end)=NaN;avg.EPSILON(end)=NaN;avg.EPSILON(end)=NaN;
        %     end
        %     if q.script.num>=544 && q.script.num<=549
        %         avg.EPSILON=avg.EPSILON1;
        %     end
        %     %         avg.WD(idsur)=NaN; avg.WD2(idsur)=NaN;
        %
        %     %             epsilon_glitch_factor=6;
        %     %             avg.EPSILON=(avg.EPSILON1+avg.EPSILON2)/2;
        %     %             % determine if EPSILON1>>>EPSILON2
        %     %             a=find(avg.EPSILON1>epsilon_glitch_factor*avg.EPSILON2 | isnan(avg.EPSILON1));
        %     %             avg.EPSILON(a)=avg.EPSILON2(a);
        %     %             % determine if EPSILON2>>>EPSILON1
        %     %             a=find(avg.EPSILON2>epsilon_glitch_factor*avg.EPSILON1 | isnan(avg.EPSILON2));
        %     %             avg.EPSILON(a)=avg.EPSILON1(a);
        
        %
        %    warning backtrace
        %
        % calc dynamic height
        %
        %    head=calc_dynamic_z(avg,head);
        %
        % create a seperate .mat data file containing 1m binned data and header
        %
        %     temp=num2str(q.script.num+10000);
        %     fn=[q.script.prefix temp(2:5)];
        %     head.p_max=max(cal.P);
        %     eval(['save ' path_save fn ' avg head']);
        
        % * 3/15/16 - AP - need to make all variables same length
        % make a new modified structure that I will use for chipod method
        cal2=struct();
        
        %~~ check if TP,p,c,t sampled at same rate; if not, interpolate to same
        % length; I want all to be same length as TP
        %%
        irP=head.irep.P;
        irTP=head.irep.TP;
        irT=head.irep.T;
        irS=head.irep.S;
        irFspd=head.irep.FALLSPD;
        %%
        if head.irep.P~=irTP
            clear P2
            % make new P vector length of TP (want to preserve TP higher
            % sampling rate)
            P2=nan*ones(length(cal.TP),1);
            P2( 1 : irTP/irP : end )=cal.P;
            P2=NANinterp(P2,1);
            cal2.P=P2;
            cal2.TP=cal.TP;
        end
        
        
        if irT~=irTP
            clear P2
            P2=nan*ones(length(cal.T),1);
            P2( 1 : irT/irP : end )=cal.P;
            P2=NANinterp(P2,1);
            ig=find(diffs(P2)>0);
            cal2.T=interp1(P2(ig),cal.T(ig),cal2.P);
        end
        
        if head.irep.S~=irTP
            clear P2
            P2=nan*ones(length(cal.S),1);
            P2( 1 : irS/irP : end )=cal.P;
            P2=NANinterp(P2,1);
            ig=find(diffs(P2)>0);
            cal2.S=interp1(P2(ig),cal.S(ig),cal2.P);
            
        end
        
        if head.irep.FALLSPD~=irTP
            clear P2
            P2=nan*ones(length(cal.FALLSPD),1);
            P2( 1 : irFspd/irP : end )=cal.P;
            P2=NANinterp(P2,1);
            ig=find(diffs(P2)>0);
            cal2.fspd=interp1(P2(ig),cal.FALLSPD(ig),cal2.P);
            
        end
        %~~
        
        % 03/18/16 - AP - save N2 and dTdz also
        cal2.N2=interp1(cal.P,cal.N2,cal2.P);
        cal2.dTdz=interp1(cal.P,cal.DTDZ,cal2.P);
        cal2.Note='N2 and dTdz computed in cali_eq08.m'
        cal.MakeInfo=['Made by AP ' datestr(now) ' w/ run_eq08_AP.m'];
        cal2.MakeInfo=['Made by AP ' datestr(now) ' w/ run_eq08_AP.m'];
        
        % * 3/15/16 - AP - need to save calibrated cast files here
        clear temp fn
        temp=num2str(q.script.num+10000);
        fn=[q.script.prefix temp(2:5)];
        save(fullfile(path_save,fn),'cal2','head')
        
    end % try
    
end % cast #
delete(hb)

% sum_eq08_uppend
%sum_eq08

%%