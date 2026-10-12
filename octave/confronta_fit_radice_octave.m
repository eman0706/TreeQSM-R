% Passo 2 della prova di unita' di cylinder_fitting: legge i file scritti da esporta_fit_radice (R) ed esegue il cylinder_fitting ORIGINALE.
% Uso (in Octave): cd('<cartella prova_fit>'); addpath(genpath('<TreeQSM-master>/src')); addpath('<cartella con cylinder_fitting_orig.m>'); confronta_fit_radice_octave
1;
warning('off', 'all'); pkg load statistics;
P = csvread('fit_P.csv', 1, 0); Points = csvread('fit_Points.csv', 1, 0); Ind = csvread('fit_Ind.csv', 1, 0); nl = size(Ind, 1);
fprintf('P %d punti, Points %d, strati %d\n', size(P, 1), numel(Points), nl);
[cyl, Reg] = cylinder_fitting_orig(P, Points, Ind, nl, 1);
M = double([cyl.radius cyl.length cyl.start cyl.axis cyl.SurfCov cyl.mad double(cyl.conv) double(cyl.rel)]);
dlmwrite('fit_Oct_out.csv', M, 'precision', 12);
fprintf('cilindri: %d\n  z_start   lunghezza   D(cm)   SurfCov\n', size(M, 1)); disp([M(:, 5) M(:, 2) 200 * M(:, 1) M(:, 9)](1:min(25, end), :));
