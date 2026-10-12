import json, csv, collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
J = json.load(open('/home/claude/albero_chiamate.json')); E = J['archi']; eff = J['effettiva']; defs = J['definizioni']; reach = set(J['raggiunte'])
# stato di audit (dal lavoro svolto): V verificata riga per riga e fedele; C difetto trovato e corretto; D deviazione documentata; N non letta riga per riga
STATO = {}
for f in ('cover_sets', 'generate_cover_uniform', 'generate_cover_variable', 'segments', 'define_cut', 'cut_components', 'study_components', 'component_classification',
          'correct_segments', 'modify_topology', 'remove_small', 'modify_parent', 'search_stem_top', 'search_branch_top', 'segment_direction', 'cylinder_fitting', 'adjustments', 'build_cover_structure', 'create_prioritized_permutation'): STATO[f] = 'V'
for f in ('relative_size', 'define_main_branches', 'surface_coverage'): STATO[f] = 'C'
for f in ('define_trunk', 'define_base_forb', 'least_squares_cylinder', 'surface_coverage_filtering', 'cylinders', 'tree_data', 'filter_layer_outliers'): STATO[f] = 'D'
for f in ('parent_cylinder', 'branches', 'point_model_distance', 'make_tree_connected', 'cubical_partition', 'cubical_partition_with_cubes', 'least_squares_circle_centre', 'surface_coverage2',
          'func_grad_cylinder', 'growth_volume_correction', 'treeqsm', 'tree_sets'): STATO[f] = 'N' if f not in ('tree_sets', 'treeqsm') else 'V'
ERR = {'cylinders', 'cylinder_fitting', 'adjustments'}      # funzioni che catturano errori e proseguono (misurati: zero errori su 3 alberi)
COL = {'V': '#7bc47f', 'C': '#6fa8dc', 'D': '#f6b26b', 'N': '#d9d9d9'}
ETI = {'V': 'letta riga per riga, fedele', 'C': 'difetto trovato e corretto', 'D': 'deviazione documentata', 'N': 'non letta riga per riga'}
ETI_EN = {'V': 'read line by line, faithful', 'C': 'defect found and corrected', 'D': 'documented deviation', 'N': 'not read line by line'}
# verifica che ogni arco mostrato esista davvero nel grafo estratto dal codice
fasi = [('cover_sets', ['generate_cover_uniform', 'generate_cover_variable', 'cubical_partition']), ('relative_size', []),
        ('tree_sets', ['define_trunk', 'define_base_forb', 'define_main_branches', 'make_tree_connected']),
        ('segments', ['define_cut', 'cut_components', 'study_components', 'component_classification']),
        ('correct_segments', ['modify_topology', 'remove_small', 'modify_parent']),
        ('cylinders', ['cylinder_fitting', 'parent_cylinder', 'adjustments', 'filter_layer_outliers']),
        ('branches', []), ('tree_data', ['robust_trunk_reference', 'diameter_at_height']), ('point_model_distance', [])]
mancanti = [(a, b) for a, bs in fasi for b in bs if b not in E.get(a, {})]
secondo = {'modify_topology': ['search_stem_top', 'search_branch_top'], 'cylinder_fitting': ['surface_coverage_filtering', 'least_squares_cylinder'], 'least_squares_cylinder': ['surface_coverage'], 'adjustments': ['least_squares_circle_centre', 'surface_coverage2']}
mancanti += [(a, b) for a, bs in secondo.items() for b in bs if b not in E.get(a, {})]
print("archi mostrati ma assenti dal grafo estratto dal codice:", mancanti)
fig, ax = plt.subplots(figsize=(13.2, 7.6)); ax.set_xlim(0, 13.2); ax.set_ylim(0, 7.6); ax.axis('off')
def spezza(t, n=15):
    if len(t) <= n or '_' not in t: return t
    parti = t.split('_'); best = None
    for k in range(1, len(parti)):
        a = '_'.join(parti[:k]) + '_'; b = '_'.join(parti[k:])
        sc = abs(len(a) - len(b))
        if best is None or sc < best[0]: best = (sc, a + '\n' + b)
    return best[1]
def box(x, y, w, h, testo, st, size=8.2, bold=False, err=False):
    testo = spezza(testo)
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02,rounding_size=0.06', fc=COL[st], ec='#444', lw=0.8))
    ax.text(x + w/2, y + h/2, testo + (' ᴱ' if err else ''), ha='center', va='center', fontsize=size, fontweight='bold' if bold else 'normal', family='DejaVu Sans')
ax.plot([xs0 + 0.67 for xs0 in (0.15,)] + [0.15 + 8 * 1.44 + 0.67], [7.12, 7.12], color='#999', lw=0.8)
ax.text(0.15, 7.35, 'treeqsm()', fontsize=11, fontweight='bold'); ax.text(1.55, 7.35, '— pipeline of the R implementation (arrows: calls extracted from the source code)', fontsize=9)
xs = [0.15 + i * 1.44 for i in range(9)]; w = 1.34
for i, (f, subs) in enumerate(fasi):
    x = xs[i]; box(x, 6.45, w, 0.6, f, STATO.get(f, 'N'), size=7.2, bold=True, err=f in ERR)
    ax.plot([x + w / 2, x + w / 2], [7.12, 7.05], color='#999', lw=0.8)
    for j, s in enumerate(subs):
        yy = 5.75 - j * 0.62; box(x + 0.14, yy, w - 0.14, 0.5, s, STATO.get(s, 'N'), size=6.6, err=s in ERR)
        ax.plot([x + 0.07, x + 0.14], [yy + 0.25, yy + 0.25], color='#777', lw=0.7)
    if subs: ax.plot([x + 0.07, x + 0.07], [6.45, 5.75 - (len(subs) - 1) * 0.62 + 0.25], color='#777', lw=0.7)
# secondo livello: dettaglio di modify_topology e dei cilindri
ax.text(0.15, 3.05, 'Second level', fontsize=9.5, fontweight='bold')
dettaglio = [('modify_topology', 0.15), ('cylinder_fitting', 3.15), ('least_squares_cylinder', 6.15), ('adjustments', 9.15)]
for f, x0 in dettaglio:
    box(x0, 2.45, 2.7, 0.45, f, STATO.get(f, 'N'), size=8.0, bold=True, err=f in ERR)
    for j, s in enumerate(secondo[f]):
        box(x0 + 0.25, 1.8 - j * 0.58, 2.4, 0.45, s, STATO.get(s, 'N'), size=7.4, err=s in ERR)
        ax.plot([x0 + 0.12, x0 + 0.25], [2.025 - j * 0.58, 2.025 - j * 0.58], color='#777', lw=0.7)
    ax.plot([x0 + 0.12, x0 + 0.12], [2.45, 2.025 - (len(secondo[f]) - 1) * 0.58], color='#777', lw=0.7)
ax.text(0.15, 0.62, 'dbh_cylinder (faithful DBHcyl) is not called by treeqsm: it is called only by tree_data_core, outside the pipeline.   ᴱ = catches errors and continues (0 errors in a run on three trees)', fontsize=7.8, style='italic', color='#444')
# legenda
for k, (st, et) in enumerate(ETI_EN.items()):
    box(0.15 + k * 3.2, 0.05, 0.3, 0.28, '', st); ax.text(0.52 + k * 3.2, 0.19, et, fontsize=7.8, va='center')
plt.savefig('/home/claude/albero_chiamate.png', dpi=210, bbox_inches='tight'); plt.savefig('/home/claude/albero_chiamate.svg', bbox_inches='tight')
print("figura salvata")
# --- albero testuale annotato ---
righe = []; mostrate = set()
def stampa(n, liv, via):
    st = STATO.get(n, 'N'); sigla = {'V': '[verificata]', 'C': '[corretta]', 'D': '[deviazione]', 'N': '[non letta]'}[st] + (' [errori inghiottiti]' if n in ERR else '')
    ripetuta = n in mostrate
    if ripetuta and not E.get(n): return                       # foglia gia' mostrata: non si ripete
    righe.append(f"{'    ' * liv}{n}  ({eff.get(n, ('?'))})  {sigla}" + ('  (gia\' mostrata sopra)' if ripetuta else ''))
    if ripetuta: return
    mostrate.add(n)
    for c in sorted(E.get(n, {})):
        if c in ('set.seed',): continue
        stampa(c, liv + 1, via + [n])
stampa('treeqsm', 0, [])
open('/mnt/user-data/outputs/albero_chiamate.txt', 'w', encoding='utf-8').write(
    "ALBERO DELLE CHIAMATE DI treeqsm (R), generato dal codice sorgente. File = definizione effettiva (l'ultima caricata).\n"
    "Stato: verificata = letta riga per riga e fedele a MATLAB; corretta = difetto trovato e corretto; deviazione = deviazione documentata; non letta = non letta riga per riga.\n\n" + "\n".join(righe) + "\n")
print(len(righe), "righe nell'albero testuale")
# --- inventario per il repository ---
with open('/mnt/user-data/outputs/funzioni_inventario.csv', 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f); w.writerow(['funzione', 'raggiunta_da_treeqsm', 'n_definizioni', 'file_definizione_effettiva', 'tutte_le_definizioni', 'stato_audit'])
    for n in sorted(defs): w.writerow([n, 'si' if n in reach else 'no', len(defs[n]), eff[n], ' | '.join(defs[n]), STATO.get(n, '') ])
print("inventario scritto:", len(defs), "funzioni")
