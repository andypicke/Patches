%
%
% Comparing patch calculations of N2 and gamma for different methods.
% Specifically, I was computing the bulk N2 using just temperature (from
% Bill's paper). But I added a 2nd version where I compute the bulk
% gradient using sgth (so including temp and salinity).
%
%
%
%%

figure(1);clf

h1=histogram(log10(patches.gam_bulk(:)));
hold on
histogram(log10(patches.gam_bulk_2(:)),h1.BinEdges,'EdgeColor','none')

%%

figure(2);clf
histogram2( log10(patches.gam_bulk(:)), log10(patches.gam_bulk_2(:)), 'DisplayStyle','Tile')

%%

figure(2);clf
histogram2( log10(patches.n2_bulk(:)), log10(patches.n2_bulk_2(:)), 'DisplayStyle','Tile')
xlim([-6 -2])
ylim([-6 -2])

%%