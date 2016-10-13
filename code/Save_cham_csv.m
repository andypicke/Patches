%%

% save a csv file to look at in R

%%

clear ; close all

dir_base='/Users/Andy/Cruises_Research/ChiPod/Cham_Eq14_Compare/Data/chameleon/processed';

load( fullfile( dir_base, '/Cstar=0_01366/sum/eq14_sum.mat') )

chi=cham.CHI(:);
eps=cham.EPSILON(:);
n2=cham.N2(:);
dtdz=cham.DTDZ(:);

M=[n2(:) dtdz(:) chi(:) eps(:)];

csvwrite('test.csv',M)

%%