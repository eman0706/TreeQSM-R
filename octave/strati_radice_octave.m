% STRATI_RADICE_OCTAVE: struttura degli strati del segmento radice (originale). Uso:  albero_diag = 35; strati_radice_octave
1;
warning('off', 'all'); pkg load statistics;
if ~exist('albero_diag', 'var'), albero_diag = 35; end
if ~exist('cart_diag', 'var'), cart_diag = 'alberi_corretti_v2'; end
create_input; inputs.plot = 0; inputs.disp = 0; inputs.savemat = 0; inputs.savetxt = 0; inputs.name = 'albero'; inputs.tree = 1; inputs.model = 1; inputs.filter.plot = 0;
P = csvread(fullfile(cart_diag, sprintf('albero_%03d.csv', albero_diag)), 1, 0); P = double(P(:, 1:3)); P = P(P(:, 3) >= 0.15, :);
zq = [0.25 0.5 1 1.5 2]; res = zeros(0, 5); np10 = zeros(0, 10); ns10 = zeros(0, 10); info = zeros(0, 5);
for h = 1:length(inputs.PatchDiam1)
  Inputs = inputs; Inputs.PatchDiam1 = inputs.PatchDiam1(h); Inputs.BallRad1 = inputs.BallRad1(h);
  cover1 = cover_sets(P, Inputs); [cover1, Base, Forb] = tree_sets(P, cover1, Inputs);
  segment1 = segments(cover1, Base, Forb); segment1 = correct_segments(P, cover1, segment1, Inputs, 0, 1, 1);
  for i = 1:length(inputs.PatchDiam2Max)
    Inputs.PatchDiam2Max = inputs.PatchDiam2Max(i); Inputs.BallRad2 = inputs.BallRad2(i);
    for j = 1:length(inputs.PatchDiam2Min)
      Inputs.PatchDiam2Min = inputs.PatchDiam2Min(j);
      RS = relative_size(P, cover1, segment1); cover2 = cover_sets(P, Inputs, RS);
      [cover2, Base, Forb] = tree_sets(P, cover2, Inputs, segment1);
      segment2 = segments(cover2, Base, Forb); segment2 = correct_segments(P, cover2, segment2, Inputs, 1, 1, 0);
      SPar = segment2.ParentSegment; bases = find(SPar(:, 1) == 0); root = bases(1);
      Seg = segment2.segments{root}; nl = numel(Seg); cen = double(cover2.center);
      nset = zeros(nl, 1); zc = zeros(nl, 1); npt = zeros(nl, 1);
      for l = 1:nl
        S = double(Seg{l}(:)); nset(l) = numel(S); zc(l) = mean(P(cen(S), 3));
        for q = 1:numel(S), npt(l) = npt(l) + numel(cover2.ball{S(q)}); end
      end
      r = NaN(1, 5); for q = 1:5, k = find(zc >= zq(q), 1); if ~isempty(k), r(q) = k; end, end
      res(end+1, :) = r; np10(end+1, :) = npt(1:10)'; ns10(end+1, :) = nset(1:10)'; info(end+1, :) = [nl sum(nset) sum(npt) zc(1) max(zc)];
    end
  end
end
fprintf('\nStrato (indice) in cui la quota media dei centri raggiunge z = %s, per ciascun modello:\n', mat2str(zq)); disp(res);
fprintf('Strati, insiemi, punti del segmento radice; quota del primo e dell''ultimo strato:\n'); disp(round(100 * info) / 100);
fprintf('\nMEDIANE: strati a z: %s | punti nei primi 10 strati: %s | insiemi nei primi 10 strati: %s\n', mat2str(median(res, 1)), mat2str(median(np10, 1)), mat2str(median(ns10, 1)));
