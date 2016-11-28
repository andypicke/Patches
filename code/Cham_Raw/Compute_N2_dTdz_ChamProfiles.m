%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Compute_N2_dTdz_ChamProfiles.m
%
% Compute N2 and dT/dz for overturns in EQ14 chameleon profiles using a few
% different methods to compare the results.
%
% See also LookAtProfilesOT.m
%
%--------------
% 11/22/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

% Make empty arrays for results
n1=[];
nb=[];
n3=[];
n4=[];

dtdz1=[];
dtdz2=[];
dtdz3=[];

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/

for cnum=1:3100 % choose a cast # to look at
    
    cnum
    
    try
        
        % load chameleon cast
        clear cal cal2 head
        cham_dir='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/Cham_proc_AP/cal';
        load( fullfile( cham_dir, ['eq14_' sprintf('%04d',cnum) '.mat'] ) )
        
        clear s t p lat
        s=cal.SAL(1:end-1); % (end-1) b/c last 2 values are same;
        t=cal.T1(1:end-1);
        p=cal.P(1:end-1);
        
        clear idot lat1 lat2
        idot=strfind(head.lat.start,'.');
        lat1=str2num(head.lat.start(1:idot-3));
        lat2=str2num(head.lat.start(idot-2:end))/60;
        lat=nanmean([lat1 lat2]);
        
        % Find overturns
        clear Params
        Params.lat=lat;
        Params.plotit=0;
        Params.sigma=1e-4;
        Params.runlmin=0;
        Params.minotsize=1;
        Params.usetemp=1;
        addpath /Users/Andy/Standard-Mixing-Routines/ThorpeScales/
        clear OT
        OT=compute_overturns_discrete_AP(p,t,s,Params);
        
        clear pstarts pstops
        pstarts=OT.pstarts_each;
        pstops=OT.pstops_each;
        
        for iot=1:length(OT.pstarts_each)
            
            clear iz t_ot s_ot p_ot
            iz=isin(p,[ OT.pstarts_each(iot) OT.pstops_each(iot) ]);
            t_ot=t(iz);
            s_ot=s(iz);
            p_ot=p(iz);
            
            clear t_sort I
            [t_sort I]=sort(t_ot,1,'descend');
            
            clear DT dz dTdz
            dT=nanmax(t_ot)-nanmin(t_ot);
            dz=nanmax(p_ot)-nanmin(p_ot);
            dTdz=-dT/dz;
            
            b=p_ot(end) - (1/dTdz)*t_ot(end);
            tvec=linspace(nanmin(t_ot),nanmax(t_ot),100);
            pvec=linspace(nanmin(p_ot),nanmax(p_ot),100);
            
            yvec=dTdz*pvec - dTdz*b;
            
            % fit a line
            P=polyfit(p_ot,t_sort,1);
            
            % save results
            dtdz1=[dtdz1 ; dTdz];
            dtdz2=[dtdz2 ; P(1) ];
            
            %~~ 'bulk gradient' method from Smyth et al 2001
            % essentially rms T (btw sorted/unsorted) /  thorpe scale
            t_rms= sqrt( nanmean(( t_ot - t_sort ).^2) );
            dtdz3=[ dtdz3 ; t_rms / OT.Lt_each(iot) ];
            
            % Now do similar for density / N^2
            
            clear sgth sgth_ot
            sgth=sw_pden(s,t,p,0);
            sgth_ot=sw_pden(s_ot,t_ot,p_ot,0);
            
            clear sgth_sort I
            [sgth_sort I]=sort(sgth_ot,1,'ascend');
            
            % try the range/dz method
            drho=nanmax(sgth_ot)-nanmin(sgth_ot);
            dz=nanmax(p_ot)-nanmin(p_ot);
            n2_1=9.81/nanmean(sgth_ot)*drho/dz;
            
            n1=[ n1 ; n2_1];
            
            
            % fit a line to sgth
            clear P1
            P1=polyfit(p_ot,sgth_sort,1);
            
            % calculate N^2 from this fit and compare to other N2 estimates
            clear drhodz n2_2 drho dz n2_3
            drhodz=P1(1);
            n2_2=9.81/nanmean(sgth_ot)*drhodz;
            
            nb=[nb ; n2_2];
            
            % compute N^2
            clear n2
            %           n2=sw_bfrq(s_ot,t_ot,p_ot,1);
            n2=sw_bfrq(s_ot(I),t_ot(I),p_ot(I),1);
            
            n3=[n3 ; nanmean(n2)];
            
            % Smyth eta al says associated N^2=g*bulk gradient (assuming
            % density controlled by temperature??)
            % I think missing divide by rho_0 also?
            n4=[n4 ; 9.81/nanmean(sgth_ot)*t_rms];            

        end %iot
        
    end % try
    
end % cnum


%%

figure(1);clf
agutwocolumn(1)
wysiwyg
set(gcf,'defaultaxesfontsize',15)

subplot(211)
%loglog(abs(dtdz1(:)),abs(dtdz2(:)),'.')
histogram2(log10(abs(dtdz1(:))),log10(abs(dtdz2(:))),50,'displaystyle','tile')
xlabel('log_{10}[dtdz1]')
ylabel('log_{10}[dtdz2]')
%grid on
grid on

subplot(212)
%loglog(abs(dtdz1(:)),abs(dtdz2(:)),'.')
histogram2(log10(abs(dtdz1(:))),log10(abs(dtdz3(:))),50,'displaystyle','tile')
xlabel('log_{10}[dtdz1]')
ylabel('log_{10}[dtdz3]')
%grid on
grid on

figdir='/Users/Andy/Cruises_Research/ChiPod/Analyses/Patch_n2_dTdz'
print( fullfile( figdir, ['eq14_cham_patch_dTdzs_scatter'] ), '-dpng' )

%%

figure(1);clf
agutwocolumn(0.7)
wysiwyg
set(gcf,'defaultaxesfontsize',14)

h1=histogram(log10(abs(dtdz1(:))),'DisplayStyle','stair')%,'edgecolor','none');
hold on
h2=histogram(log10(abs(dtdz2(:))),h1.BinEdges,'DisplayStyle','stair')%,'edgecolor','none');
h3=histogram(log10(abs(dtdz3(:))),h1.BinEdges,'DisplayStyle','stair')%,'edgecolor','none');
grid on
xlim([-4 -0.5])
xlabel('log_{10}[dT/dz]')
ylabel('count')
legend([h1 h2 h3],'dtdz1','dtdz2','dtdz3')

figdir='/Users/Andy/Cruises_Research/ChiPod/Analyses/Patch_n2_dTdz'
print( fullfile( figdir, ['eq14_cham_patch_dTdzs'] ), '-dpng' )

%%

figure(2);clf
agutwocolumn(1)
wysiwyg
set(gcf,'defaultaxesfontsize',15)

xl=[-5.5 -3]
yl=[-6.75 -2.5]
m=3
n=1

subplot(m,n,1)
%loglog(n1(:),nb(:),'.')
histogram2(real(log10(n1(:))),real(log10(nb(:))),50,'displaystyle','tile')
grid on
xlim(xl)
ylim(yl)

xlabel('log_{10}[N^2_1]')
ylabel('log_{10}[N^2_2]')

subplot(m,n,2)
%loglog(n1(:),n3(:),'.')
histogram2(real(log10(n1(:))),real(log10(n3(:))),50,'displaystyle','tile')
grid on
xlabel('log_{10}[N^2_1]')
ylabel('log_{10}[N^2_3]')
xlim(xl)
ylim(yl)

subplot(m,n,3)
%loglog(n1(:),n3(:),'.')
histogram2(real(log10(n1(:))),real(log10(n4(:))),50,'displaystyle','tile')
grid on
xlabel('log_{10}[N^2_1]')
ylabel('log_{10}[N^2_4]')
xlim([xl])
ylim([yl])

print( fullfile( figdir, ['eq14_cham_patch_N2s_scatter'] ), '-dpng' )
%subplot(313)

%%

figure(2);clf
agutwocolumn(0.7)
wysiwyg
set(gcf,'defaultaxesfontsize',14)

h1=histogram(real(log10(n1(:)))) ;
hold on
h2=histogram(real(log10(nb(:))),h1.BinEdges) ;
h3=histogram(real(log10(n3(:))),h1.BinEdges) ;
h4=histogram(real(log10(n4(:))),h1.BinEdges) ;
grid on
xlim([-6 -2.5])
xlabel('log_{10}[N^2]')
ylabel('count')
legend([h1 h2 h3 h4],'1  ','2 ','3 ','4 ','location','best')

print( fullfile( figdir, ['eq14_cham_patch_N2s'] ), '-dpng' )

%% Compute gamma using different N^2 and dT/dz values

% first, need to get chi and epsilon for each patch

%clear ; close all

min_patch_size=1
usetemp=1

data_dir=fullfile('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/',['minOT_' num2str(10*min_patch_size) '_usetemp_' num2str(usetemp)]) 

%n2_all=[];
%n2_OT_all=[];
%dtdz_all=[];
chi_all=[];
eps_all=[];
P_all=[];

hb=waitbar(0,'compiling patch data from all profiles');
for cnum=1:3100
    waitbar(cnum/3100,hb)
    clear avg
    try
        fname=['eq14_' sprintf('%04d',cnum) '.mat'];
        load(fullfile(data_dir,fname))
        
        if length(avg.N2)~=length(avg.n2_OT)
            disp(num2str(cnum))
        end
        
        P_all     = [P_all     ; avg.P(:)    ] ;
%        n2_all    = [n2_all    ; avg.N2(:)   ] ;
 %       n2_OT_all = [n2_OT_all ; avg.n2_OT(:)] ;
  %      dtdz_all  = [dtdz_all  ; avg.dTdz(:) ] ;
        chi_all   = [chi_all   ; avg.CHI(:)  ] ;
        eps_all   = [eps_all   ; avg.EPSILON(:) ];
        
    end
    
end % cnum
delete(hb)
%%
% use 'bulk gradient' for both

% use polyfit for both

% use range/dz for both

% plot histograms of each

%%