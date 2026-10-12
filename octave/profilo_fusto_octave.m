% PROFILO_FUSTO_OCTAVE: diametro (cm) del fusto a varie lunghezze della catena di cilindri, dopo e PRIMA di adjustments (UnmodRadius), per gli 8 modelli.
% Uso (in Octave):  albero_diag = 35; profilo_fusto_octave        (cart_diag = 'alberi_corretti_v2' e' il predefinito)
1;
warning('off', 'all'); pkg load statistics;
if ~exist('albero_diag', 'var'), albero_diag = 35; end
if ~exist('cart_diag', 'var'), cart_diag = 'alberi_corretti_v2'; end
create_input;
inputs.plot = 0; inputs.disp = 0; inputs.savemat = 0; inputs.savetxt = 0; inputs.name = 'albero'; inputs.tree = 1; inputs.model = 1; inputs.filter.plot = 0;
P = csvread(fullfile(cart_diag, sprintf('albero_%03d.csv', albero_diag)), 1, 0); P = double(P(:, 1:3)); P = P(P(:, 3) >= 0.15, :);
quote = [0.25 0.5 0.75 1 1.25 1.5 2 3 4 6];
PatchDiam1 = inputs.PatchDiam1; PatchDiam2Min = inputs.PatchDiam2Min; PatchDiam2Max = inputs.PatchDiam2Max; BallRad1 = inputs.BallRad1; BallRad2 = inputs.BallRad2;
fin = zeros(0, numel(quote)); unm = zeros(0, numel(quote)); sc = zeros(0, numel(quote)); ncl = zeros(0, 3);
for h = 1:length(PatchDiam1)
  Inputs = inputs; Inputs.PatchDiam1 = PatchDiam1(h); Inputs.BallRad1 = BallRad1(h);
  cover1 = cover_sets(P, Inputs); [cover1, Base, Forb] = tree_sets(P, cover1, Inputs);
  segment1 = segments(cover1, Base, Forb); segment1 = correct_segments(P, cover1, segment1, Inputs, 0, 1, 1);
  for i = 1:length(PatchDiam2Max)
    Inputs.PatchDiam2Max = PatchDiam2Max(i); Inputs.BallRad2 = BallRad2(i);
    for j = 1:length(PatchDiam2Min)
      Inputs.PatchDiam2Min = PatchDiam2Min(j);
      RS = relative_size(P, cover1, segment1); cover2 = cover_sets(P, Inputs, RS);
      [cover2, Base, Forb] = tree_sets(P, cover2, Inputs, segment1);
      segment2 = segments(cover2, Base, Forb); segment2 = correct_segments(P, cover2, segment2, Inputs, 1, 1, 0);
      cylinder = cylinders(P, cover2, segment2, Inputs);
      T = find(cylinder.branch == 1); L = double(cylinder.length(T)); cs = cumsum(L); r = double(cylinder.radius(T)); u = double(cylinder.UnmodRadius(T)); s = double(cylinder.SurfCov(T));
      f = NaN(1, numel(quote)); g = f; w = f;
      for q = 1:numel(quote), k = find(cs >= quote(q), 1); if ~isempty(k), f(q) = 200 * r(k); g(q) = 200 * u(k); w(q) = s(k); end, end
      fin(end+1, :) = f; unm(end+1, :) = g; sc(end+1, :) = w; ncl(end+1, :) = [numel(T) sum(abs(r - u) > 1e-6) mean(L)];
    end
  end
end
fprintf('\nDiametro DOPO adjustments (cm), colonne = lunghezza catena %s m\n', mat2str(quote)); disp(round(10 * fin) / 10);
fprintf('Diametro PRIMA di adjustments (UnmodRadius, cm)\n'); disp(round(10 * unm) / 10);
fprintf('SurfCov del cilindro\n'); disp(round(100 * sc) / 100);
fprintf('Cilindri del fusto | cilindri modificati da adjustments | lunghezza media (m)\n'); disp(round(100 * ncl) / 100);
fprintf('\nMEDIANE: dopo:   %s\n', mat2str(round(10 * median(fin)) / 10)); fprintf('MEDIANE: prima:  %s\n', mat2str(round(10 * median(unm)) / 10)); fprintf('MEDIANE: SurfCov: %s\n', mat2str(round(100 * median(sc)) / 100));
