%
%
%
%
%
%%

figure(1);clf

h1=histogram(log10(patches.gam_bin(:)), 'Normalization','pdf');
hold on
h2=histogram(log10(patches.gam_line(:)), 'Normalization','pdf');
h3=histogram(log10(patches.gam_bulk(:)), 'Normalization','pdf');
freqline(log10(0.2))
xlim([-3.5 1])
grid on
legend([h1 h2 h3],'bin','line','bulk')
xlabel('log_{10}[\Gamma]')

SetNotesFigDir
print( fullfile(NotesFigDir,'eq14_gamma_hist_line_bulk'), '-dpng')


%%

figure(1);clf
h1=histogram( real(log10(patches.dtdz_bin(:))), 'Normalization','pdf');
hold on
h2=histogram(log10(patches.dtdz_line(:)), 'Normalization','pdf');
h3=histogram(log10(patches.dtdz_bulk(:)), 'Normalization','pdf');
freqline(log10(0.2))
xlim([-3.5 1])
grid on
legend([h1 h2 h3],'bin','line','bulk')
xlabel('log_{10}[T_z]')


%%

figure(1);clf
h1=histogram( real(log10(patches.n2_bin(:))), 'Normalization','pdf');
hold on
h2=histogram(log10(patches.n2_line(:)), 'Normalization','pdf');
h3=histogram(log10(patches.n2_bulk(:)), 'Normalization','pdf');
freqline(log10(0.2))
xlim([-6.5 -2.5])
grid on
legend([h1 h2 h3],'bin','line','bulk')
xlabel('log_{10}[N^2]')


%%