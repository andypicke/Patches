function [h h1 h2]=Plot4hist_2(n2,dtdz,chi,eps,n2_2,dtdz_2,chi_2,eps_2)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Plot4hist_2.m
%
% Same as Plot4hist.m, but plot 2 different distributions to compare.
%
% Make a 2X2 figure w/ histograms of N^2,dT/dz,chi, and eps
%
% Was doing this a lot for this analysis so wrote function to simplify
% thigns
%
%---------------------
% 10/25/16 - A.Pickering - andypicke@gmail.com
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
Nm='pdf';

h=figure(1);clf
agutwocolumn(1)
wysiwyg

subplot(221)
h1=histogram(real(log10(n2(:))),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(real(log10(n2_2(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}N^2')
grid on
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')
ylabel(Nm)

subplot(222)
h1=histogram(real(log10(dtdz(:))),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(real(log10(dtdz_2(:))),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}dT/dz')
grid on
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')
ylabel(Nm)

subplot(223)
h1=histogram(log10(chi(:)),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(log10(chi_2(:)),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}\chi')
grid on
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')
ylabel(Nm)

subplot(224)
h1=histogram(log10(eps(:)),'edgecolor','none','Normalization',Nm);
hold on
h2=histogram(log10(eps_2(:)),h1.BinEdges,'edgecolor','none','Normalization',Nm);
xlabel('log_{10}\epsilon')
freqline(nanmedian(h1.Data),'b')
freqline(nanmedian(h2.Data),'r')
grid on
ylabel(Nm)

%%