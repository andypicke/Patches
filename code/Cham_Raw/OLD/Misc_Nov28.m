

%%

clear ; close all

load('/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed/Cstar=0_032/sum/eq14_sum_clean.mat')

addpath /Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/mfiles/Patches/code/

% compute gamma from binned data
gam_bin = ComputeGamma(cham.N2,cham.DTDZ,cham.CHI,cham.EPSILON);

ib=find(log10(cham.EPSILON)<-8.5);
gam_bin(ib)=nan;

figure(1);clf
histogram( log10(gam_bin(:)) , 'edgecolor','none')
xlim([-5 2])

nanmedian(gam_bin(:))

%%
%%

figure(1);clf
loglog(patches.n2_bin(:),patches.n3(:),'.')
loglog(patches.dtdz_bin(:),patches.dtdz3(:),'.')
loglog(patches.chi_bin(:),patches.chi(:),'.')
loglog(patches.eps_bin(:),patches.eps(:),'.')
grid on


%%

figure(2);clf
%h1=histogram(patches.n3 - patches.n2_bin);xlim(2e-4*[-1 1])
%h1=histogram(patches.dtdz3 - patches.dtdz_bin);xlim(0.03*[-1 1])
h1=histogram(patches.chi ./ patches.chi_bin);%xlim(0.03*[-1 1])
h1=histogram(patches.eps ./ patches.eps_bin);%xlim(0.03*[-1 1])
freqline(nanmedian(h1.Data))


%%

%%

addpath /Users/Andy/Cruises_Research/seawater_ver3_2/

figure(1);clf
plot(cham.N2(:,100),cham.P(:,100))
axis ij
n2=sw_bfrq(cham.SAL(:,100),cham.T1(:,100),cham.P(:,100),nanmean(cham.lat))
hold on
plot(n2,cham.P(1:end-1,100))
%%