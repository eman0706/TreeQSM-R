# DBH DIRETTO ALLE TRE QUOTE DI CAMPO (1.30, 2.00, 4.00 m) SULLE DUE NUVOLE, PYTHON
# Il file riferimento_alberi.csv contiene i diametri di campo a 1.30 m (media dei due diametri: ref_cm), a 2.00 m (d200) e a 4.00 m (d400).
# Stesso algoritmo e stessi parametri del DBH diretto a 1.3 m (fetta 0.2 m, Hough + fit geometrico robusto entro +-4 cm), solo la quota cambia.
# Greco et al. (2023, Tabella 9) riportano RMSE e MAPE MEDIATI sulle quote misurate: nuvola 15: 0.88 cm e 2.63%; nuvola 11: 0.92 cm e 2.73%.
# PREVISIONI SCRITTE PRIMA DELL'ESECUZIONE
#  P_Z1: la mediana |errore| del diretto a 2.00 m e a 4.00 m e' <= 3.0% su entrambe le nuvole.
#  P_Z2: l'errore a 4.00 m e' maggiore di quello a 1.30 m (mediana |errore| a 4.00 m >= quella a 1.30 m) su entrambe le nuvole
#        (fusto piu' inclinato, nuvola piu' rada, meno alberi misurati).
#  P_Z3: la MAPE del diretto mediata sulle tre quote e' <= 2.63% (nuvola 15) e <= 2.73% (nuvola 11), cioe' non peggiore della procedura di Greco et al.
import sys, csv, numpy as np
sys.path.insert(0, '/home/claude'); import diretto_py as dp
def centra(P, h, half=0.15):
    # ricentra x, y sul fusto alla quota h: i tronchi sono inclinati e l'estrazione e' centrata sull'asse a 1.3 m, non a 2 o 4 m
    S = P[np.abs(P[:, 2] - h) < half]
    if len(S) < 100: return None
    c = np.median(S[:, :2], axis=0)
    for _ in range(4):
        m = np.hypot(S[:, 0] - c[0], S[:, 1] - c[1]) < 0.45
        if m.sum() < 60: return None
        A = np.column_stack([2 * S[m, 0], 2 * S[m, 1], np.ones(m.sum())]); b = (S[m, :2] ** 2).sum(1); s = np.linalg.lstsq(A, b, rcond=None)[0]; c = s[:2]
    Q = P.copy(); Q[:, 0] -= c[0]; Q[:, 1] -= c[1]; return Q
def val(x):
    try: v = float(x); return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError): return np.nan
risultati = {}
for nome, cart, gref in (('nuvola 15', '/home/claude/pacchetto_matlab/alberi_corretti_v2', (0.88, 2.63)), ('nuvola 11', '/home/claude/pacchetto_matlab/alberi_corretti_las11', (0.92, 2.73))):
    rif = [r for r in csv.DictReader(open(cart + '/riferimento_alberi.csv')) if r['fonte_ref'] == 'campo']
    righe = []
    for r in rif:
        n = int(r['N']); P = np.loadtxt(f'{cart}/albero_{n:03d}.csv', delimiter=',', skiprows=1)[:, :3]
        v = {'N': n}
        for h, k in ((1.3, 'ref_cm'), (2.0, 'd200'), (4.0, 'd400')):
            ref = val(r[k]); v[f'ref{h}'] = ref
            Q = centra(P, h)
            try: d = dp.dbh_diretto(Q, h=h)[0] if Q is not None else np.nan
            except Exception: d = np.nan
            if not (5.0 <= d <= 80.0): d = np.nan   # stima non valida: non entra nelle statistiche, ma e' contata
            v[f'dir{h}'] = d
        righe.append(v)
    risultati[nome] = righe
    print(f"\n=== {nome} ===")
    m_rmse, m_mape = [], []
    for h in (1.3, 2.0, 4.0):
        e = [(v[f'dir{h}'], v[f'ref{h}']) for v in righe if np.isfinite(v[f'ref{h}']) and np.isfinite(v[f'dir{h}'])]
        nref = sum(1 for v in righe if np.isfinite(v[f'ref{h}'])); nval = len(e)
        d = np.array([a for a, b in e]); rf = np.array([b for a, b in e]); er = 100 * (d - rf) / rf
        rmse = np.sqrt(np.mean((d - rf) ** 2)); mape = np.mean(np.abs(er)); m_rmse.append(rmse); m_mape.append(mape)
        print(f"  {h:.2f} m: n = {len(e):>2} su {nref} con riferimento | mediana |err| {np.median(np.abs(er)):.2f}% | media |err| {mape:.2f}% | segno mediano {np.median(er):+.1f}% | RMSE {rmse:.2f} cm | oltre 15%: {int((np.abs(er) > 15).sum())}")
    print(f"  media sulle tre quote (come Greco et al., Tabella 9): RMSE {np.mean(m_rmse):.2f} cm | MAPE {np.mean(m_mape):.2f}%   (Greco et al.: {gref[0]} cm | {gref[1]}%)")
import json; json.dump(risultati, open('/home/claude/quote/dbh_diretto_quote.json', 'w'))
