%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ExamineRawPatchDataEq14.m
%
%
%
% 11/01/16 - A.Pickering
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

clear ; close all

n2=[];
dtdz=[]
chi=[];
eps=[];

for cnum=1:3100
    clear avg
    try
        load(['/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/ChamRawProc/eq14_' sprintf('%04d',cnum) '.mat'])
        %        avg1=avg;clear avg
        n2=[n2 ; avg.N2(:)];
        dtdz=[dtdz ; avg.dTdz(:) ];
        chi=[chi ; avg.CHI(:) ];
        eps=[eps ; avg.EPSILON(:) ];
        
    end
    
end

%%

h=Plot4hist(n2,dtdz,chi,eps)
%%

gam=ComputeGamma(n2,dtdz,chi,eps);

ig=find(gam<1);
figure(1);clf
histogram(gam(ig))
freqline(nanmedian(gam))

%%