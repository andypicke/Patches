%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute_N2_dTdz_ChamProfiles_V2.m
%
% Compute N2 and dT/dz for overturns in EQ14 chameleon profiles using a few
% different methods.
%
% Uses patches found in Find_Patches_EQ14_Raw.m
%
% OUTPUT:
% 'patches' structure
%
% See also LookAtProfilesOT.m
%
% - Run FindPatches_EQ14_Raw.m 1st to identify patches
% - then run run_eq14_for_PATCHES.m to process cham data over these patches
%
%--------------
% 11/22/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% patch options
patch_size_min = 0.25  % min patch size
usetemp = 1

% load patch data (from FindPatches_EQ14_Raw.m)
datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/'
fname=['EQ14_raw_patches_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '.mat']
load(fullfile(datdir,fname))

patches = struct() ;
patches.cnum = patch_data(:,1) ;
patches.p1   = patch_data(:,2) ;
patches.p2   = patch_data(:,3) ;
patches.n2_ot= patch_data(:,5) ;
patches.Lt   = patch_data(:,6) ;

% Make empty arrays for results
Npatches=length(patches.p1);
EmpVec=nan*ones(Npatches,1);

% Different methods of computing N^2
patches.n2_range=EmpVec;
patches.n2_line=EmpVec;
patches.n2_bulk=EmpVec;
patches.n4=EmpVec;
patches.n2_bulk_2=EmpVec;

patches.drhodz_bulk=EmpVec;
patches.drhodz_line=EmpVec;

% Different methods of computing dT/dz
patches.dtdz_range=EmpVec;
patches.dtdz_line=EmpVec;
patches.dtdz_bulk=EmpVec;

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/

cnum_loaded=-5;
warning off
hb=waitbar(0);

for ip=1:Npatches
    
    waitbar(ip/Npatches,hb)
    
    try
        
        clear cnum
        cnum=patches.cnum(ip);
        
        % check if this cast already loaded to avoid unnecessary re-loading
        if cnum_loaded~=cnum
            % load chameleon cast
            clear cal cal2 head
            cham_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal';
            load( fullfile( cham_dir, ['eq14_' sprintf('%04d',cnum) '.mat'] ) )
            cnum_loaded=cnum;
            
            % compute pot. temp, pot. density etc.
            clear s t p lat ptmp sgth
            s = cal.SAL(1:end-1); % (end-1) b/c last 2 values are same;
            s = smooth(s,20);
            t = cal.T1 (1:end-1);
            p = cal.P  (1:end-1);
            ptmp=sw_ptmp(s,t,p,0);
            sgth=sw_pden(s,t,p,0);
            
            % get latitude for profile
            clear idot lat1 lat2
            idot=strfind(head.lat.start,'.');
            lat1=str2num(head.lat.start(1:idot-3));
            lat2=str2num(head.lat.start(idot-2:end))/60;
            lat=nanmean([lat1 lat2]);
            
        else        
            % this profile already loaded
        end
        
        
        clear pstarts pstops
        clear iz t_ot s_ot p_ot ptmp_ot sgth_ot
        
        % get indices of data in patch
        iz=isin(p,[ patches.p1(ip) patches.p2(ip) ]);
        
        if length(iz)>10
            
            t_ot=t(iz);
            s_ot=s(iz);
            p_ot=p(iz);
            ptmp_ot=ptmp(iz);
            sgth_ot=sgth(iz);
            
            clear t_sort I
            [t_sort , I]=sort(t_ot,1,'descend');
            
            % sort potential temp
            clear t_pot t_pot_sort
            [ptmp_sort , Iptmp]=sort(ptmp_ot,1,'descend');
            
            clear DT dz dTdz
            dT=nanmax(ptmp_ot)-nanmin(ptmp_ot);
            dz=nanmax(p_ot)-nanmin(p_ot);
            dTdz=-dT/dz;
                         
            % fit a line
            P=polyfit(p_ot,ptmp_sort,1);
            
            % save results
            patches.dtdz_range(ip)=dTdz;
            
            patches.dtdz_line(ip)=-P(1);
            
            %~~ 'bulk gradient' method from Smyth et al 2001
            % essentially = rms T (btw sorted/unsorted) /  thorpe scale ?
            %    t_rms= sqrt( nanmean(( t_ot - t_sort ).^2) );
            t_rms= sqrt( nanmean(( ptmp_ot - ptmp_sort ).^2) );
            patches.dtdz_bulk(ip)=  t_rms / patches.Lt(ip) ;
            
            %~~ Now do similar for density / N^2
            
            [sgth_sort , I]=sort(sgth_ot,1,'ascend');
            
            % try the range/dz method
            drho=nanmax(sgth_ot)-nanmin(sgth_ot);
            dz=nanmax(p_ot)-nanmin(p_ot);
            n2_1=9.81/nanmean(sgth_ot)*drho/dz;
            
            patches.n2_range(ip)=n2_1;
            
            % fit a line to sgth
            clear P1
            P1=polyfit(p_ot,sgth_sort,1);
            
            % calculate N^2 from this fit
            clear drhodz n2_2 drho dz n2_3
            drhodz=-P1(1);
            patches.drhodz_line(ip)=drhodz;
            
            n2_2=-9.81/nanmean(sgth)*drhodz;
            patches.n2_line(ip)=n2_2;
            
            % Smyth eta al says associated N^2=g*bulk gradient (assuming
            % density controlled by temperature??)
            % I think missing divide by rho_0 also?
            %patches.n2_bulk(ip)= 9.81 / nanmean(sgth) * t_rms / patches.Lt(ip)    ;
            alpha = -0.2641; % from fit of sgth vs theta
            patches.n2_bulk(ip)= -9.81 / nanmean(sgth) * alpha * t_rms / patches.Lt(ip)    ;
            
            % try computing drho/dz using 'bulk' method (instead of
            % assuming rho only depends on T
            clear rho_rms
            rho_rms = sqrt( nanmean(( sgth_ot - sgth_sort ).^2) );
            patches.drhodz_bulk(ip) = - rho_rms / patches.Lt(ip) ;
            patches.n2_bulk_2(ip) = - 9.81 / nanmean(sgth) * patches.drhodz_bulk(ip);
            
            % compute N^2 w/ sw_bfrq
            clear n2
            n2=sw_bfrq(s_ot(I),t_ot(I),p_ot,0.5);
            patches.n4(ip)=nanmean(n2);
            
        end
        
    end % try
    
end %ip
delete(hb)

warning on

%% Add (patch-wise) Chameleon chi and epsilon for each patch

% these are calulated in run_eq14_for_PATCHES.m

data_dir=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc/',['minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp)])

chi_all = [] ;
eps_all = [] ;
P_all = [] ;
cnum_all = [] ;

patches.eps=nan*ones(size(patches.p1));
patches.chi=nan*ones(size(patches.p1));

cnums=unique(patches.cnum);

hb=waitbar(0,'compiling patch data from all profiles');

for ic=1:length(cnums)
    waitbar(ic/length(cnums),hb)
    
    clear avg cnum
    
    cnum=cnums(ic);
    
    try
        % load the processed profile w/ chi and eps
        fname=['eq14_' sprintf('%04d',cnum) '.mat'];
        load(fullfile(data_dir,fname))
        
        % now loop over all patches and get average Eps and chi from those
        % depths
        clear ig
        ig=find(patches.cnum==cnum);
        for ip=1:length(ig)
            
            clear iz
            iz=isin( avg.P,[patches.p1(ig(ip)) patches.p2(ig(ip))]);
            
            if size(iz)==1
                patches.eps(ig(ip)) = nanmean(avg.EPSILON(iz)) ;
                patches.chi(ig(ip)) = nanmean(avg.CHI(iz))     ;
            end
            
        end %ip
        
    end % try
    
end % cnum
delete(hb)

%% Exclude values where epsilon is below noise floor
% note this makes significant difference in gamma median
ib=find(log10(patches.eps)<-8.5);
patches.eps(ib)=nan;

%%
addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/
% use 'bulk gradient' for both
gam_range = ComputeGamma(patches.n2_range, patches.dtdz_range, patches.chi, patches.eps);
% use polyfit for both
gam_line = ComputeGamma(patches.n2_line, patches.dtdz_line, patches.chi, patches.eps);
% use bulk for both
gam_bulk = ComputeGamma(patches.n2_bulk, patches.dtdz_bulk, patches.chi, patches.eps);
%
gam_bulk_2 = ComputeGamma(patches.n2_bulk_2, patches.dtdz_bulk, patches.chi, patches.eps);
%
gam4 = ComputeGamma(patches.n4, patches.dtdz_line, patches.chi , patches.eps );

% save results so we don't have to run everything again
patches.gam_range = gam_range;
patches.gam_line  = gam_line;
patches.gam_bulk  = gam_bulk;
patches.gam_bulk_2 = gam_bulk_2;
patches.gam4 = gam4;

save( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )


%% Get binned data (N^2,dTdz,chi,eps, gamma) at patch locations
% I want to compare gamma from binned data to the gamma I compute over
% patches. So for each patch I want to find the closest binned gamma and
% then compare those.

% load binned chameleon data
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/

Npatches = length(patches.cnum) ;
patches.gam_bin = nan*ones(size(patches.gam_bulk)) ;
patches.n2_bin  = nan*ones(size(patches.gam_bulk)) ;
patches.dtdz_bin= nan*ones(size(patches.gam_bulk)) ;
patches.chi_bin = nan*ones(size(patches.gam_bulk)) ;
patches.eps_bin = nan*ones(size(patches.gam_bulk)) ;
patches.drhodz_bin=nan*ones(size(patches.gam_bulk)) ;

for ip=1:Npatches
    
    clear cnum pbin pmn val I Icham
    cnum = patches.cnum(ip) ;
    Icham= find(cham.castnumber==cnum) ;
    pbin = cham.P(:,Icham) ;
    pmn  = nanmean([patches.p1(ip) patches.p2(ip)]) ;
    
    clear ig
    ig=find(~isnan(pbin)) ;
    patches.n2_bin(ip)    = interp1(pbin(ig) , cham.N2(ig,Icham) , pmn);
    patches.dtdz_bin(ip)  = interp1(pbin(ig) , cham.DTDZ_RHOORDER(ig,Icham) , pmn);
    patches.chi_bin(ip)   = interp1(pbin(ig) , cham.CHI(ig,Icham) , pmn);
    patches.eps_bin(ip)   = interp1(pbin(ig) , cham.EPSILON(ig,Icham) , pmn);
    patches.drhodz_bin(ip)= patches.n2_bin(ip) * (nanmean(cham.SIGMA(:,Icham))+1000) / -9.81 ;
    
    if log10(patches.eps_bin(ip))>-8.5
        patches.gam_bin(ip)=ComputeGamma(patches.n2_bin(ip),patches.dtdz_bin(ip),patches.chi_bin(ip),patches.eps_bin(ip));
    end
    
    %  [val,I]=nanmin( abs(pbin-pmn)) ;
    
    %    patches.n2_bin(ip)  = cham.N2(I,Icham) ;
    %   patches.dtdz_bin(ip)= cham.DTDZ_RHOORDER(I,Icham) ;
    %patches.chi_bin(ip) = cham.CHI(I,Icham) ;
    %patches.eps_bin(ip) = cham.EPSILON(I,Icham) ;
    %    patches.drhodz_bin(ip)= cham.N2(I,Icham) * (nanmean(cham.SIGMA(:,Icham))+1000) / -9.81;
    
    %     if log10(cham.EPSILON(I,Icham))>-8.5
    %         patches.gam_bin(ip)=ComputeGamma(cham.N2(I,Icham),cham.DTDZ_RHOORDER(I,Icham),cham.CHI(I,Icham),cham.EPSILON(I,Icham));
    %     end
    
    
end

%% save again

patches.MakeInfo = ['Made ' datestr(now) ' w/ Compute_N2_dTdz_ChamProfiles_V2.m'] 

save( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/data/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(100*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

%%

