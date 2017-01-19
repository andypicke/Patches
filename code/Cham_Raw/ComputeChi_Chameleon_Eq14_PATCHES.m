%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ComputeChi_Chameleon_Eq14_PATCHES.m
%
% ** Modified version to do chi-pod method for just patches in chameleon
% profiles **
%
% Modified from ComputeChi_Chameleon_Eq14.m
%
% Chameleon mat files are made w/ ProcessEq14Cham_AP.m, which must be run
% first.
%
% After this, combine all profiles w/ Combine_EQ14_chi.m
%
% Modified from ComputeChi_Chameleon_Flx91.m
%
%
%------------------------
% 12/14/16 - A.Pickering - apickering@coas.oregostate.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% Add all the paths we need from mixing software
mixpath='/Users/Andy/Cruises_Research/mixingsoftware/'
addpath(fullfile(mixpath,'seawater'))
addpath(fullfile(mixpath,'general')) % kraichnan.m
addpath(fullfile(mixpath,'marlcham')) % integrate.m
addpath(fullfile(mixpath,'Chameleon2','Version2015')) % fast_psd
addpath(fullfile(mixpath,'CTD_Chipod','mfiles')) ;
addpath(fullfile(mixpath,'chipod','compute_chi')); % get_chipod_chi.md

makeplots=0
savespec =0  % save wavenumber spectra

patch_size_min = 0.25
usetemp = 1

load( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

% directory for chameleon casts we have (made w/ ProcessEq14Cham_AP.m)
datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal'

% Params for chipod calculations
Params.z_smooth = 0   ;   % distance to smooth N^2 and dTdz over
Params.nfft     = 128 ;   %
Params.fmax     = 7   ;   %
Params.TPthresh = 1e-6;   %
Params.resp_corr= 0   ;   % correct TP spectra for freq response of thermistor
Params.fc       = 99  ;   % cutoff frequency for response correction
Params.gamma    = 0.2 ;   % mixing efficiency

% option to use gamma computed in patches, instead of a constant value
use_patch_gam=0;

if Params.resp_corr==0
    Params.fc=99;
end

whN2dTdz='bulk'
%whN2dTdz='bulk2'

% Make directory to save processed casts in (name based on Params)
if use_patch_gam==1
    datdirsave=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches',...
        ['N2dTdz_' num2str(whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gammaPATCH_nfft_' num2str(Params.nfft)]);
else
    datdirsave=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChipodPatches',...
        ['N2dTdz_' num2str(whN2dTdz) '_fmax' num2str(Params.fmax) 'Hz_respcorr' num2str(Params.resp_corr) '_fc_' num2str(Params.fc) 'hz_gamma' num2str(Params.gamma*100) '_nfft_' num2str(Params.nfft)]);
end

disp(['Data will be saved to ' datdirsave])

% check if directory exists, make new if not
ChkMkDir(datdirsave)

%%
tstart=tic;
% loop through each cast
for cnum=[4:12 14:46 48:87 374:519 550:597 599:904 906:909 911:1070 ...
        1075:1128 1130:1737 1739:2550 2552:2996 2998:3092]
        
    disp(['Working on cnum ' num2str(cnum) ])
    
    close all
    clear cal head
    clear tpspec kspec kkspec fspec kks ks
    
    try
        %%
        % Load the data for this cast
        load( fullfile(datdir,['EQ14_' sprintf('%04d',cnum) '.mat']) )
        
        clear cal
        cal=cal2;
        clear cal2
        
        %%
        clear ctd z_smooth
        ctd=struct();
        ctd.t1=cal.T1;
        ctd.s1=cal.SAL;
        ctd.p=cal.P;
        
        % add in lat and lon (from sum_eq14.m)
        clear idot lat1 lat2
        idot=strfind(head.lon.start,'.');
        lon1=str2num(head.lon.start(1:idot-3));
        lon2=str2num(head.lon.start(idot-2:end))/60;
        ctd.lon=nanmean([lon1 lon2]);
        
        clear idot lat1 lat2
        idot=strfind(head.lat.start,'.');
        lat1=str2num(head.lat.start(1:idot-3));
        lat2=str2num(head.lat.start(idot-2:end))/60;
        ctd.lat=nanmean([lat1 lat2]);
        
        ctd.lat=nanmean([str2num(head.lat.start) str2num(head.lat.end)]);
        
        %% find patches for this cnum
        clear igc
        igc=find(patches.cnum==cnum);
        
        %% get start/stop indices for the patch depths
        
        Np=length(igc);
        todo_inds=nan*ones(Np,2);
        for ip=1:Np
            clear iz
            iz=isin(cal.P,[patches.p1(igc(ip)) patches.p2(igc(ip)) ] );
            todo_inds(ip,:)=[iz(1) iz(end)] ;
        end
        
        Nwindows=Np;
        
        %% need to replace this, returns todo_inds, the indices of data for each window (start/stop)
        % Make windows for chi calcs
        clear TP
        TP=cal.TP;
        %[todo_inds,Nwindows]=MakeCtdChiWindows(TP,Params.nfft);
        
        %%
        
        %~ make 'avg' structure for the processed data
        clear avg
        avg=struct();
        avg.Params=Params;
        tfields={'datenum','P','N2','dTdz','fspd','T','S','P','theta','sigma',...
            'chi1','eps1','KT1','TP1var','kstart','kstop'};
        for n=1:length(tfields)
            avg.(tfields{n})=NaN*ones(Nwindows,1);
        end
        
        avg.FspdStd=avg.fspd;
        
        % Get TP samplerate (Hz)
        avg.samplerate=head.slow_samp_rate * head.irep.TP;
        
        % Get average time, pressure, and fallspeed in each window
        for iwind=1:Nwindows
            clear inds
            inds=todo_inds(iwind,1) : todo_inds(iwind,2);
            avg.P(iwind)=nanmean(cal.P(inds));
            avg.fspd(iwind)=nanmean(cal.FALLSPD(inds))/100;
            avg.FspdStd(iwind)=nanstd(cal.FALLSPD(inds)/100);
        end
        
        %%
        %~ save spectra to plot after
        % make empty arrays for spectra
        if savespec==1
            % observed wavenumber
            ks=nan*ones(Nwindows,Params.nfft /2);
            % observed wave# spectra
            kspec=nan*ones(Nwindows,Params.nfft /2);
            % fit wavenumber
            kks=nan*ones(Nwindows,56);
            % fit wave# spectra
            kkspec=nan*ones(Nwindows,56);
            % observed frequency spectra
            tpspec=nan*ones(Nwindows,Params.nfft /2);
            %~
        end
        
        %%
        
        if makeplots==1
            %~~ plot histogram of avg.P to see how many good windows we have in
            %each 10m bin
            if nanmax(avg.P>10)
                figure
                hi=histogram(avg.P,0:10:nanmax(avg.P));
                hi.Orientation='Horizontal';axis ij;
                ylabel('P [db]')
                xlabel('# good data windows')
                % title([whSN ' cast ' cast_suffix ' - ' C.castdir 'cast'],'interpreter','none')
                % print('-dpng',fullfile(chi_fig_path_specific,[whSN '_' cast_suffix '_Fig' num2str(whfig) '_' C.castdir 'cast_chi_' whsens '_avgPhist']))
                % whfig=whfig+1;
            end
        end
        
        %%
        
        % get N2, dTdz for each window
        good_inds=find(~isnan(ctd.p));
        % interpolate ctd data to same pressures as chipod
        %    avg.N2  =interp1(ctd.p(good_inds),ctd.N2(good_inds),avg.P);
        %   avg.dTdz=interp1(ctd.p(good_inds),ctd.dTdz(good_inds),avg.P);
        
        % Get N2 and dTdz from patches structure instead
        
        if whN2dTdz==2
            avg.N2   = patches.nb(igc)   ;
            avg.dTdz = patches.dtdz2(igc);
            avg.gamma= patches.gam2(igc) ;
        elseif strcmp(whN2dTdz,'bulk')
            avg.N2   = patches.n2_bulk(igc)   ;
            avg.dTdz = patches.dtdz_bulk(igc);
            avg.gamma= patches.gam_bulk(igc) ;
        elseif strcmp(whN2dTdz,'bulk2')
            avg.N2   = patches.n2_bulk_2(igc)   ;
            avg.dTdz = patches.dtdz_bulk(igc);
            avg.gamma= patches.gam_bulk(igc) ;
        end
        
        % add binned eps and chi so we can compare after
        avg.chi_bin   = patches.chi_bin(igc) ;
        avg.chi_patch = patches.chi(igc) ;
        
        avg.eps_bin  =patches.eps_bin(igc);
        avg.eps_patch=patches.eps(igc)    ;
        
        % ** find unique ctd.p values
        [p2,IA,IC] = unique(ctd.p(good_inds));
        avg.T   =interp1(ctd.p(good_inds(IA)),ctd.t1(good_inds(IA)),avg.P);
        avg.S   =interp1(ctd.p(good_inds(IA)),ctd.s1(good_inds(IA)),avg.P);
                
        
        %%
        % note sw_visc not included in newer versions of sw?
        avg.nu=sw_visc_ctdchi(avg.S,avg.T,avg.P);
        avg.tdif=sw_tdif_ctdchi(avg.S,avg.T,avg.P);
        
        avg.n_iter=nan*ones(size(avg.P));
        %%
        % loop through each window and do the chi computation
        for iwind=1:Nwindows
            
            clear inds
            % get inds for this window
            inds=todo_inds(iwind,1) : todo_inds(iwind,2);
            
            % integrate dT/dt spectrum
            clear tp_power freq
            [tp_power,freq]=fast_psd(TP(inds),Params.nfft,avg.samplerate);
            avg.TP1var(iwind)=sum(tp_power)*nanmean(diff(freq));
            
            if avg.TP1var(iwind)>Params.TPthresh
                
                % apply filter correction for sensor response?
                
                if Params.resp_corr==1
                    thermistor_filter_order=2; % 2=double pole?
                    thermistor_cutoff_frequency=Params.fc;%32;
                    analog_filter_order=4;
                    analog_filter_freq=50;
                    tp_power=invert_filt(freq,invert_filt(freq,tp_power,thermistor_filter_order, ...
                        thermistor_cutoff_frequency),analog_filter_order,analog_filter_freq);
                end
                
                clear gamma
                if use_patch_gam==1
                    gamma=avg.gamma(iwind);
                else
                    gamma=Params.gamma;
                end
                
                % Compute chi using iterative procedure
                [chi1,epsil1,k,spec,kk,speck,stats]=get_chipod_chi(freq,tp_power,abs(avg.fspd(iwind)),avg.nu(iwind),...
                    avg.tdif(iwind),avg.dTdz(iwind),'nsqr',avg.N2(iwind),'fmax',Params.fmax,'gamma',gamma);%,'doplots',1);
                %            pause
                %            'doplots',1 for plots
                avg.chi1(iwind)=chi1(1);
                avg.eps1(iwind)=epsil1(1);
                avg.KT1(iwind)=0.5*chi1(1)/avg.dTdz(iwind)^2;
                avg.kstart(iwind)=stats.k_start;
                avg.kstop(iwind)=stats.k_stop;
                
                avg.n_iter(iwind)=length(chi1);
                
                if savespec==1
                    % 02/17/16 - AP - save spectra
                    
                    fspec=freq;
                    tpspec(iwind,:)=tp_power;
                    
                    % observed wavenumber
                    ks(iwind,:)=k;
                    
                    % observed spectra
                    kspec(iwind,:)=spec;
                    
                    % best-fit theoreticl spectra
                    kkspec(iwind,:)=speck;
                    
                    if ~isnan(kk)
                        % theoretical fit wavenumber
                        kks(iwind,:)=kk;
                    end
                end
                
            end % if T1Pvar>threshold
            
        end % windows
        %%
        
        if makeplots==1
            %~ plot summary figure
            ax=CTD_chipod_profile_summary(avg,cal,TP);
            axes(ax(1))
            title(Flist(icast).name,'interpreter','none')
            %~~~
        end
        
        if savespec==1
            %~
            avg.tpspec=tpspec;
            avg.kspec=kspec;
            avg.kkspec=kkspec;
            avg.ks=ks;
            avg.kks=kks;
            avg.fspec=fspec;
        end
        
        avg.lat=ctd.lat;
        avg.z=sw_dpth(avg.P,avg.lat);
        avg.MakeInfo=['Made ' datestr(now) ' w/ ComputeChi_Chameleon_Eq14_PATCHES.m'];
        
        % save results
        save(fullfile(datdirsave,['EQ14_' sprintf('%04d',cnum) 'avg.mat']),'avg')
        
        
    end % try
    
end % icast

%delete(hb)
%%