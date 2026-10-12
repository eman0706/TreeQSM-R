% ESEGUI_TREEQSM_OCTAVE  (versione per GNU Octave, gratuito: https://octave.org)
% Esegue il codice ORIGINALE di TreeQSM (sorgenti MATLAB non modificati) con l'interprete GNU Octave sulle stesse nuvole usate per il porting in R,
% con i parametri di default di create_input.m. Equivale a esegui_treeqsm_matlab.m, ma usa solo funzioni di base (niente tabelle).
% TESTATO con GNU Octave 10.1.0 (Windows) e il pacchetto statistics di Octave Forge sulle nuvole 15 e 11 (risultati nel manoscritto, Sezione 4.10). Non testato con MATLAB. Funzioni non supportate: solo alphaShape (sostituita da uno stub, vedi sotto) e pdist2 (pacchetto statistics).
% Preparazione (una volta, in Octave):   pkg install -forge statistics
% Uso: addpath(genpath('percorso/TreeQSM/src')); poi, dalla cartella di lavoro:  esegui_treeqsm_octave
% Prova minima (impostazioni predefinite): solo l'albero 16 della nuvola 15. Tre alberi: alberi_prova = [4 16 35]. Corsa completa: alberi_prova = []; nuvole_da_fare = [1 2] (Octave e' piu' lento di MATLAB: non si sa quanto, la prova minima serve a misurarlo).
% USCITE: risultati_octave_<cartella>.csv (stesse colonne di quello di MATLAB, leggibile da confronta_matlab_r.R) e un riepilogo a schermo.
1;   % il file e' uno script che contiene anche una funzione
function shp = alphaShape(varargin)
  % SOSTITUTO MINIMO di alphaShape, che Octave non ha (https://www.octave.org/missing.html). TreeQSM la usa solo in tree_data.m (crown_measures) per
  % CrownAreaAlpha e CrownVolumeAlpha: queste due misure della chioma NON sono valide in questa esecuzione (NaN); diametri, cilindri, altezza e volumi del fusto non ne dipendono.
  shp.area = NaN; shp.volume = NaN;
end
warning('off', 'Octave:language-extension'); warning('off', 'Octave:mixed-string-concat');   % avvisi di sintassi inoffensivi di dlmwrite
pkg load statistics;
cartelle = {'alberi_corretti_v2', 'alberi_corretti_las11'}; rif_file = {'rif_num_v2.csv', 'rif_num_las11.csv'};
porting_F = [3.50 3.50]; procedura_greco = [3.00 3.07];
% Le impostazioni si danno dalla finestra dei comandi PRIMA di lanciare lo script (non serve modificare il file). Se non sono definite, valgono queste di PROVA, le piu' piccole:
if ~exist('alberi_prova', 'var'), alberi_prova = 16; end   % un solo albero (il piu' piccolo); [] = tutti gli alberi di campo; per esempio [4 16 35] per tre alberi
if ~exist('nrip', 'var'), nrip = 1; end                    % ripetizioni indipendenti (TreeQSM usa numeri casuali)
if ~exist('nuvole_da_fare', 'var'), nuvole_da_fare = 1; end % 1 = nuvola 15 (alberi_corretti_v2), 2 = nuvola 11 (alberi_corretti_las11), [1 2] = entrambe
create_input;             % script di TreeQSM: definisce 'inputs' con i valori di default
inputs.filter.plot = 0; inputs.plot = 0; inputs.disp = 0; inputs.savemat = 0; inputs.savetxt = 0; inputs.name = 'albero'; inputs.tree = 1; inputs.model = 1;
for c = nuvole_da_fare
  cart = cartelle{c}; file_out = ['risultati_octave_' cart '.csv'];
  rif = csvread(rif_file{c}, 1, 0);          % colonne: N, ref_cm
  alberi = rif(:, 1)';
  if ~isempty(alberi_prova), alberi = intersect(alberi, alberi_prova); end
  righe = zeros(0, 11); t0 = tic;
  for a = alberi
    P = csvread(fullfile(cart, sprintf('albero_%03d.csv', a)), 1, 0); P = P(:, 1:3); P = P(P(:, 3) >= 0.15, :);
    for r = 1:nrip
      try
        QSM = treeqsm(P, inputs);
      catch ME
        fprintf('%s albero %d, ripetizione %d: errore: %s\n   in %s, linea %d\n', cart, a, r, ME.message, ME.stack(1).name, ME.stack(1).line); continue
      end
      for m = 1:numel(QSM)
        td = QSM(m).treedata; inp = QSM(m).rundata.inputs; pmd = NaN;
        if isfield(QSM(m), 'pmdistance') && ~isempty(QSM(m).pmdistance) && isfield(QSM(m).pmdistance, 'mean'), pmd = QSM(m).pmdistance.mean; end
        righe(end+1, :) = [a, r, m, inp.PatchDiam1(1), inp.PatchDiam2Min(1), inp.PatchDiam2Max(1), 100*td.DBHcyl, 100*td.DBHqsm, td.TreeHeight, pmd, numel(QSM(m).cylinder.radius)]; %#ok<AGROW>
      end
    end
    fprintf('%s: albero %d fatto (%.1f minuti dall''inizio)\n', cart, a, toc(t0)/60); fflush(stdout);
  end
  if isempty(righe), fprintf('%s: nessun modello prodotto (vedere i messaggi di errore sopra)\n', cart); continue; end
  fid = fopen(file_out, 'w'); fprintf(fid, 'N,rip,modello,PD1,PD2min,PD2max,dbh_cyl_cm,dbh_qsm_cm,altezza_m,pmd_mean,n_cil\n'); fclose(fid);
  dlmwrite(file_out, righe, '-append', 'delimiter', ',', 'precision', 8);
  % riepilogo con la stessa regola del porting: mediana dei modelli validi (>= 5 cm) per ripetizione, poi mediana delle ripetizioni
  err = nan(numel(alberi), 1);
  fprintf('\n=== %s ===\n%4s %8s %10s %10s %10s\n', cart, 'N', 'rif_cm', 'DBH_cm', 'err_%', 'altezza_m');
  for k = 1:numel(alberi)
    a = alberi(k); v = rif(rif(:, 1) == a, 2); meds = [];
    for r = 1:nrip
      sel = righe(:, 1) == a & righe(:, 2) == r & righe(:, 7) >= 5;
      if any(sel), meds(end+1) = median(righe(sel, 7)); end %#ok<SAGROW>
    end
    if ~isempty(meds)
      est = median(meds); err(k) = 100 * (est - v) / v; fprintf('%4d %8.2f %10.1f %10.1f %10.1f\n', a, v, est, err(k), median(righe(righe(:, 1) == a, 9)));
    end
  end
  ok = ~isnan(err);
  fprintf('\nTreeQSM originale in Octave (%s): mediana |errore| %.2f%% | media %.2f%% | segno mediano %+.1f%% | oltre 15%%: %d\n', cart, median(abs(err(ok))), mean(abs(err(ok))), median(err(ok)), sum(abs(err(ok)) > 15));
  fprintf('  porting R (DBHcyl fedele, con riparazione): %.2f%% | procedura di Greco et al. (alberi isolati dagli autori): %.2f%%\n\n', porting_F(c), procedura_greco(c));
end
