% DIAGNOSI_STADI_OCTAVE: percorre i passi di treeqsm.m (originale) e stampa, per ciascuno degli 8 modelli di un albero, grandezze di ogni stadio.
% Serve a vedere in quale stadio l'originale e il porting in R cominciano a divergere (stesse colonne di diagnosi_stadi_R.R).
% Uso (in Octave, dalla cartella di lavoro, dopo aver fatto addpath(genpath('...src'))):
%   albero_diag = 35; diagnosi_stadi_octave
% (facoltativo: cart_diag = 'alberi_corretti_v2'; e' il valore predefinito). Dura pochi minuti. Non scrive nulla tranne diagnosi_stadi_octave_<N>.csv.
1;
warning('off', 'all'); pkg load statistics;
if ~exist('albero_diag', 'var'), albero_diag = 35; end
if ~exist('cart_diag', 'var'), cart_diag = 'alberi_corretti_v2'; end
create_input;
inputs.plot = 0; inputs.disp = 0; inputs.savemat = 0; inputs.savetxt = 0; inputs.name = 'albero'; inputs.tree = 1; inputs.model = 1; inputs.filter.plot = 0;
P = csvread(fullfile(cart_diag, sprintf('albero_%03d.csv', albero_diag)), 1, 0); P = double(P(:, 1:3)); P = P(P(:, 3) >= 0.15, :);
PatchDiam1 = inputs.PatchDiam1; PatchDiam2Min = inputs.PatchDiam2Min; PatchDiam2Max = inputs.PatchDiam2Max; BallRad1 = inputs.BallRad1; BallRad2 = inputs.BallRad2;
fprintf('albero %d: %d punti\n', albero_diag, size(P, 1));
fprintf('%2s %2s %2s | %6s %5s %4s | %6s %6s | %6s %5s %4s | %5s %5s | %8s %8s %8s\n', 'h', 'i', 'j', 'ns1', 'seg1', 'lay1', 'RSmed', 'RSmean', 'ns2', 'seg2', 'lay2', 'ncil', 'ncilT', 'r3base', 'D130cm', 'Lfusto');
out = zeros(0, 16);
for h = 1:length(PatchDiam1)
  Inputs = inputs; Inputs.PatchDiam1 = PatchDiam1(h); Inputs.BallRad1 = BallRad1(h);
  cover1 = cover_sets(P, Inputs);
  [cover1, Base, Forb] = tree_sets(P, cover1, Inputs);
  segment1 = segments(cover1, Base, Forb);
  segment1 = correct_segments(P, cover1, segment1, Inputs, 0, 1, 1);
  for i = 1:length(PatchDiam2Max)
    Inputs.PatchDiam2Max = PatchDiam2Max(i); Inputs.BallRad2 = BallRad2(i);
    for j = 1:length(PatchDiam2Min)
      Inputs.PatchDiam2Min = PatchDiam2Min(j);
      RS = relative_size(P, cover1, segment1);
      cover2 = cover_sets(P, Inputs, RS);
      [cover2, Base, Forb] = tree_sets(P, cover2, Inputs, segment1);
      segment2 = segments(cover2, Base, Forb);
      segment2 = correct_segments(P, cover2, segment2, Inputs, 1, 1, 0);
      cylinder = cylinders(P, cover2, segment2, Inputs);
      T = find(cylinder.branch == 1); r = double(cylinder.radius(T)); L = double(cylinder.length(T)); cs = cumsum(L);
      k = find(cs >= 1.3, 1); if isempty(k), d130 = NaN; else d130 = 200 * r(k); end
      row = [h i j numel(cover1.center) numel(segment1.segments) numel(segment1.segments{1}) median(double(RS)) mean(double(RS)) numel(cover2.center) numel(segment2.segments) numel(segment2.segments{1}) numel(cylinder.radius) numel(T) 100 * mean(r(1:min(3, numel(r)))) d130 sum(L)];
      out(end+1, :) = row;
      fprintf('%2d %2d %2d | %6d %5d %4d | %6.1f %6.1f | %6d %5d %4d | %5d %5d | %8.2f %8.2f %8.2f\n', row(1:3), row(4:6), row(7:8), row(9:11), row(12:13), row(14:16));
    end
  end
end
dlmwrite(sprintf('diagnosi_stadi_octave_%d.csv', albero_diag), out);
fprintf('\nMediane sugli 8 modelli: ns1 %.0f seg1 %.0f lay1 %.0f RSmean %.1f ns2 %.0f seg2 %.0f lay2 %.0f ncil %.0f ncilT %.0f r3base %.2f D130 %.2f Lfusto %.2f\n', median(out(:, 4:16)));
