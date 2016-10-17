%%

gam_cham=cham.N2 .* cham.CHI ./2 ./ cham.EPSILON ./ (cham.DTDZ.^2);
gam_cham1=cham.N2 .* cham.CHI ./2 ./ cham.EPSILON ./ (cham.DTDZ_RHOORDER.^2);

%%

figure(1);clf
h1=histogram(log10(gam_cham),'edgecolor','none')
hold on
h2=histogram(log10(gam_cham1),h1.BinEdges,'edgecolor','none')

%%
nanmedian(gam_cham(:))
nanmedian(gam_cham1(:))
%%

