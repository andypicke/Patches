%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute_N2_dTdz_ChamProfiles_V2.m
%
% Modified to use patches found in Find_Patches_EQ14_Raw.m instead of
% re-finding them here. This way we will have exact same patches to lineup
% with chi and eps values..
%
% Compute N2 and dT/dz for overturns in EQ14 chameleon profiles using a few
% different methods to compare the results.
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

% options
patch_size_min=0.25  % min patch size
usetemp=1

% load patch data (from FindPatches_EQ14_Raw.m)
datdir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc'
fname=['EQ14_raw_patches_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '.mat']
load(fullfile(datdir,fname))

patches=struct() ;
patches.cnum = new_patch_data(:,1) ;
patches.p1   = new_patch_data(:,2) ;
patches.p2   = new_patch_data(:,3) ;
patches.Lt   = new_patch_data(:,6) ;

% Make empty arrays for results
Npatches=length(patches.p1);
EmpVec=nan*ones(Npatches,1);

% Different methods of computing N^2
patches.n2_range=EmpVec;
patches.n2_line=EmpVec;
patches.n2_bulk=EmpVec;
patches.n4=EmpVec;

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
        else
            
        end
        %
        
        clear s t p lat
        s=cal.SAL(1:end-1); % (end-1) b/c last 2 values are same;
        t=cal.T1(1:end-1);
        p=cal.P(1:end-1);
        
        clear idot lat1 lat2
        idot=strfind(head.lat.start,'.');
        lat1=str2num(head.lat.start(1:idot-3));
        lat2=str2num(head.lat.start(idot-2:end))/60;
        lat=nanmean([lat1 lat2]);
        
        clear pstarts pstops
        clear iz t_ot s_ot p_ot
        
        iz=isin(p,[ patches.p1(ip) patches.p2(ip) ]);
        
        if length(iz)>10
            
            t_ot=t(iz);
            s_ot=s(iz);
            p_ot=p(iz);
            
            clear t_sort I
            [t_sort , I]=sort(t_ot,1,'descend');
            
            % compute potential temp
            clear t_pot t_pot_sort
            ptmp_ot=sw_ptmp(s_ot,t_ot,p_ot,0);
            [ptmp_sort , Iptmp]=sort(ptmp_ot,1,'descend');
            
            clear DT dz dTdz
            dT=nanmax(ptmp_ot)-nanmin(ptmp_ot);
            dz=nanmax(p_ot)-nanmin(p_ot);
            dTdz=-dT/dz;
            
            b=p_ot(end) - (1/dTdz)*ptmp_ot(end);
            tvec=linspace(nanmin(ptmp_ot),nanmax(ptmp_ot),100);
            pvec=linspace(nanmin(p_ot),nanmax(p_ot),100);
            yvec=dTdz*pvec - dTdz*b;
            
            % fit a line
            P=polyfit(p_ot,ptmp_sort,1);
            
            % save results
            patches.dtdz_range(ip)=dTdz;
            
            patches.dtdz_line(ip)=P(1);
            
            %~~ 'bulk gradient' method from Smyth et al 2001
            % essentially = rms T (btw sorted/unsorted) /  thorpe scale ?
            %    t_rms= sqrt( nanmean(( t_ot - t_sort ).^2) );
            t_rms= sqrt( nanmean(( ptmp_ot - ptmp_sort ).^2) );
            patches.dtdz_bulk(ip)= t_rms / patches.Lt(ip) ;
            
            %~~ Now do similar for density / N^2
            
            clear sgth sgth_ot
            sgth=sw_pden(s,t,p,0);
            sgth_ot=sw_pden(s_ot,t_ot,p_ot,0);
            
            clear sgth_sort I
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
            drhodz=P1(1);
            n2_2=9.81/nanmean(sgth_ot)*drhodz;
            
            patches.n2_line(ip)=n2_2;
            
            % Smyth eta al says associated N^2=g*bulk gradient (assuming
            % density controlled by temperature??)
            % I think missing divide by rho_0 also?
            patches.n2_bulk(ip)= 9.81 / nanmean(sgth_ot) * t_rms / patches.Lt(ip)    ;
            
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

data_dir=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/',['minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp)])

chi_all=[];
eps_all=[];
P_all=[];
cnum_all=[];

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
            
            if size(iz)>0
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
gam1=ComputeGamma(patches.n2_range,patches.dtdz_range,patches.chi,patches.eps);
% use polyfit for both
gam_line=ComputeGamma(patches.n2_line,patches.dtdz_line,patches.chi,patches.eps);
% use range/dz for both
gam_bulk=ComputeGamma(patches.n2_bulk,patches.dtdz_bulk,patches.chi,patches.eps);
% plot histograms of each
gam4=ComputeGamma(patches.n4,patches.dtdz_line,patches.chi,patches.eps);

%% save results so we don't have to run everything again

patches.gam_range=gam_range;
patches.gam_line=gam_line;
patches.gam_bulk=gam_bulk;
patches.gam4=gam4;

save( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )


%% Get binned data (N^2,dTdz,chi,eps, gamma) at patch locations
% I want to compare gamma from binned data to the gamma I compute over
% patches. So for each patch I want to find the closest binned gamma and
% then compare those.

%clear ; close all

% load binned chameleon data
%load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')
load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed_AP_7hz/sum/eq14_sum_clean.mat')

% load the patch data ( from Compute_N2_dTdz_ChamProfiles_V2.m)
%load(fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
%    'eq14_cham_patches_diffn2dtdzgamma.mat'))

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/

Npatches=length(patches.cnum)
patches.gam_bin = nan*ones(size(patches.gam1)) ;
patches.n2_bin  = nan*ones(size(patches.gam1)) ;
patches.dtdz_bin= nan*ones(size(patches.gam1)) ;
patches.chi_bin = nan*ones(size(patches.gam1)) ;
patches.eps_bin = nan*ones(size(patches.gam1)) ;

for ip=1:Npatches
    clear cnum pbin pmn val I
    cnum = patches.cnum(ip) ;
    Icham= find(cham.castnumber==cnum) ;
    pbin = cham.P(:,Icham) ;
    pmn  = nanmean([patches.p1(ip) patches.p2(ip)]) ;
    [val,I]=nanmin( abs(pbin-pmn)) ;
    
    patches.n2_bin(ip)  = cham.N2(I,Icham) ;
    patches.dtdz_bin(ip)= cham.DTDZ_RHOORDER(I,Icham) ;
    patches.chi_bin(ip) = cham.CHI(I,Icham) ;
    patches.eps_bin(ip) = cham.EPSILON(I,Icham) ;
    
    if log10(cham.EPSILON(I,Icham))>-8.5
        patches.gam_bin(ip)=ComputeGamma(cham.N2(I,Icham),cham.DTDZ_RHOORDER(I,Icham),cham.CHI(I,Icham),cham.EPSILON(I,Icham));
    end
    
end

% save again
save( fullfile( '/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc',...
    ['eq14_cham_minOT_' num2str(10*patch_size_min) '_usetemp_' num2str(usetemp) '_patches_diffn2dtdzgamma.mat']), 'patches' )

%%

