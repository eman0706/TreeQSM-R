# CONFRONTO DEI DUE SOLUTORI DEL CILINDRO AI MINIMI QUADRATI SU TRATTI DI FUSTO REALI
#   GN : come least_squares_cylinder.m (MATLAB): passo di Gauss-Newton incondizionato, par0 = [0 0 0 0 r], max 50 iterazioni,
#        convergenza se |SS0 - SS1| < 1e-4, rel = false se rcond(-A) < 10000*eps.
#   LM : come la versione effettiva del porting R: smorzamento lambda = 1e-4 * max(diag(A)), passo accettato solo se SS1 <= SS0,
#        par0 = [1e-6 1e-6 1e-6 1e-6 r], convergenza se |SS_prev - SS1| < 1e-5, rel = false se rcond(A_lm) < 1e-12.
# Dati: tratti di fusto (finestre di 0.30 e 0.45 m a quote da 0.6 a 4.0 m) dei 21 alberi della nuvola 15 e dei 20 della nuvola 11,
# con la selezione dei punti di surface_coverage_filtering (porting in Python da surface_coverage_filtering.m) prima del fit,
# come in cylinder_fitting. Asse iniziale dalla linea che unisce i centri (fit di cerchio) del 40% inferiore e superiore.
# PREVISIONI SCRITTE PRIMA DELL'ESECUZIONE
#  P_S1: differenza assoluta mediana dei raggi tra GN e LM < 0.5 mm, 95 percentile < 2 mm.
#  P_S2: nei tratti che contengono 1.3 m (finestra 1.15-1.45 m) la differenza dei raggi e' < 1 mm in almeno 19 alberi su 21 (nuvola 15).
#  P_S3: la quota di fit convergenti e' >= 95% per entrambi i solutori.
#  P_S4: la mediana dell'errore del diametro del solo tratto 1.15-1.45 m contro il campo non differisce tra GN e LM di piu' di 0.3 punti.
import numpy as np, csv, glob, os, sys
EPS = np.finfo(float).eps
def dist_line(P, axis, start):
    V = P - start; h = V @ axis; V = V - np.outer(h, axis); return np.linalg.norm(V, axis=1), V, h
def orthonormal(axis):
    a = axis / np.linalg.norm(axis); t = np.array([1., 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1., 0])
    u = np.cross(a, t); u /= np.linalg.norm(u); w = np.cross(a, u); return u, w
def rot_to_z(axis):
    D = np.cross(axis, [0, 0, 1.]); n = np.linalg.norm(D)
    if n <= 0: return np.eye(3)
    k = D / n; th = np.arccos(np.clip(axis[2], -1, 1)); K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) * np.cos(th) + np.sin(th) * K + (1 - np.cos(th)) * np.outer(k, k)
def form_rot(alpha, beta):
    ca, sa, cb, sb = np.cos(alpha), np.sin(alpha), np.cos(beta), np.sin(beta)
    R = np.array([[cb, 0, sb], [sa*sb, ca, -sa*cb], [-ca*sb, sa, ca*cb]])
    D1 = np.array([[0, 0, 0], [ca*sb, -sa, -ca*cb], [sa*sb, ca, sa*cb]]); D2 = np.array([[-sb, 0, cb], [sa*cb, 0, sa*sb], [-ca*cb, 0, -ca*sb]])
    return R, D1, D2
def fg(par, P):
    x0, y0, al, be, r = par; R, D1, D2 = form_rot(al, be); Q = P - np.array([x0, y0, 0.]); Pt = Q @ R.T
    xt, yt = Pt[:, 0], Pt[:, 1]; rt = np.sqrt(xt**2 + yt**2); d = rt - r; rs = np.maximum(rt, 1e-10); N = np.column_stack([xt/rs, yt/rs])
    J = np.zeros((len(P), 5)); A1 = R @ np.array([-1., 0, 0]); A2 = R @ np.array([0, -1., 0])
    J[:, 0] = N[:, 0]*A1[0] + N[:, 1]*A1[1]; J[:, 1] = N[:, 0]*A2[0] + N[:, 1]*A2[1]
    J[:, 2] = (N * (Q @ D1.T)[:, :2]).sum(1); J[:, 3] = (N * (Q @ D2.T)[:, :2]).sum(1); J[:, 4] = -1
    return d, J
def rcond(A):
    try: return 1.0 / (np.linalg.norm(A, 1) * np.linalg.norm(np.linalg.inv(A), 1))
    except np.linalg.LinAlgError: return 0.0
def fit(P, axis, start, r0, metodo):
    Pt = (P - start) @ rot_to_z(axis).T
    if metodo == 'GN':
        par = np.array([0, 0, 0, 0, r0], float); conv = False; rel = True; it = 0
        while it < 50 and not conv and rel:
            d0, J = fg(par, Pt); SS0 = np.linalg.norm(d0); A = J.T @ J; b = J.T @ d0
            try: p = -np.linalg.solve(A, b)
            except np.linalg.LinAlgError: p = -np.linalg.lstsq(A, b, rcond=None)[0]
            par = par + p
            if rcond(-A) < 10000 * EPS: rel = False
            d1, _ = fg(par, Pt); SS1 = np.linalg.norm(d1)
            if abs(SS0 - SS1) < 1e-4: conv = True
            it += 1
        return par[4], conv, rel, it
    par = np.array([1e-6]*4 + [r0], float); conv = False; rel = True; it = 0; SSp = np.inf
    while it < 50 and not conv:
        d0, J = fg(par, Pt); A = J.T @ J; b = J.T @ d0; SS0 = np.linalg.norm(d0)
        lam = 1e-4 * np.max(np.diag(A)); Al = A + lam * np.eye(5)
        try: p = -np.linalg.solve(Al, b)
        except np.linalg.LinAlgError: rel = False; break
        pn = par + p; d1, _ = fg(pn, Pt); SS1 = np.linalg.norm(d1)
        if SS1 <= SS0: par = pn
        if rcond(Al) < 1e-12: rel = False
        if abs(SSp - SS1) < 1e-5: conv = True
        SSp = SS1; it += 1
    return par[4], conv, rel, it
def scf(P, axis, start, length, lh=0.02, ns=20):
    d, V, h = dist_line(P, axis, start); h = h - h.min(); U, W = orthonormal(axis); Vp = V @ np.column_stack([U, W]); ang = np.arctan2(Vp[:, 1], Vp[:, 0]) + np.pi
    def celle(nl, ns):
        L = np.ceil(h / length * nl).astype(int); L[L == 0] = 1; L[L > nl] = nl
        S = np.ceil(ang / 2 / np.pi * ns).astype(int); S[S == 0] = 1; return L + (S - 1) * nl
    nl = int(np.ceil(length / lh)); lex = celle(nl, ns); Dis = {}
    for c in np.unique(lex):
        D = d[lex == c].min(); Dis[c] = min(1.05 * D, D + 0.02)
    R = np.median(list(Dis.values())); a = max(0.02, 0.2 * R); ns = min(36, max(int(np.ceil(2 * np.pi * R / a)), 8)); nl = max(int(np.ceil(length / a)), 3)
    lex = celle(nl, ns); r = max(0.01, 0.05 * R); Pass = np.zeros(len(P), bool); Dis = {}; k = 0
    for c in np.unique(lex):
        m = lex == c; D = d[m]; Dm = D.min(); Pass[np.where(m)[0][D <= Dm + r]] = True; Dis[c] = min(1.05 * Dm, Dm + 0.02); k += 1
    return Pass, np.median(list(Dis.values())), k / nl / ns
def kasa(xy):
    A = np.column_stack([2*xy[:, 0], 2*xy[:, 1], np.ones(len(xy))]); b = (xy**2).sum(1); s = np.linalg.lstsq(A, b, rcond=None)[0]
    return s[0], s[1], np.sqrt(max(s[2] + s[0]**2 + s[1]**2, 1e-12))
def main():
    righe = []
    for nuvola, cart in (('15', '/home/claude/pacchetto_matlab/alberi_corretti_v2'), ('11', '/home/claude/pacchetto_matlab/alberi_corretti_las11')):
        ref = [r for r in csv.DictReader(open(os.path.join(cart, 'riferimento_alberi.csv'))) if r['fonte_ref'] == 'campo']
        for r in ref:
            n = int(r['N']); P = np.loadtxt(os.path.join(cart, f'albero_{n:03d}.csv'), delimiter=',', skiprows=1)[:, :3]; refd = float(r['ref_cm'])
            sezioni = {}
            for L in (0.30, 0.45):
                for z0 in (0.6, 0.9, 1.15, 1.5, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0):
                    m = (P[:, 2] >= z0) & (P[:, 2] < z0 + L); Q = P[m]
                    if len(Q) < 300: continue
                    # centro e raggio dal fit di cerchio sui punti entro 0.6 m dall'asse approssimativo (centro della fetta), poi selezione dell'anello
                    c0 = np.median(Q[:, :2], axis=0)
                    for _ in range(3):
                        dd = np.hypot(Q[:, 0]-c0[0], Q[:, 1]-c0[1]); sel = dd < 0.45; cx, cy, rr = kasa(Q[sel, :2]); c0 = np.array([cx, cy])
                    dd = np.hypot(Q[:, 0]-c0[0], Q[:, 1]-c0[1]); sel = dd < rr + 0.08; Q = Q[sel]
                    if len(Q) < 200: continue
                    lo = Q[Q[:, 2] < z0 + 0.4*L]; up = Q[Q[:, 2] >= z0 + 0.6*L]
                    if len(lo) < 60 or len(up) < 60: continue
                    cl = np.array([*kasa(lo[:, :2])[:2], lo[:, 2].mean()]); cu = np.array([*kasa(up[:, :2])[:2], up[:, 2].mean()])
                    ax = (cu - cl); ax = ax / np.linalg.norm(ax)
                    if np.degrees(np.arccos(ax[2])) > 35: continue
                    start = np.array([cl[0], cl[1], cl[2]]); h = (Q - start) @ ax; start = start + h.min() * ax; length = h.max() - h.min()
                    Pass, R, sc = scf(Q, ax, start, length)
                    if Pass.sum() < 20 or sc < 0.2: continue
                    out = {'nuvola': nuvola, 'N': n, 'L': L, 'z0': z0, 'n': int(Pass.sum()), 'SC': sc, 'ref': refd}
                    for met in ('GN', 'LM'):
                        rad, conv, rel, it = fit(Q[Pass], ax, start, R, met); out[met] = rad; out[met + '_conv'] = bool(conv and rel); out[met + '_it'] = it
                    righe.append(out)
    return righe
if __name__ == '__main__':
    rs = main(); print("tratti analizzati:", len(rs), "| per nuvola:", {k: sum(1 for r in rs if r['nuvola'] == k) for k in ('15', '11')})
    d = np.array([abs(r['GN'] - r['LM']) * 1000 for r in rs if r['GN_conv'] and r['LM_conv'] and r['GN'] > 0 and r['LM'] > 0])
    print(f"tratti convergenti con entrambi: {len(d)} | differenza assoluta dei raggi (mm): mediana {np.median(d):.3f} | 95 percentile {np.percentile(d, 95):.3f} | massimo {d.max():.2f}")
    rel = np.array([100 * (r['GN'] - r['LM']) / r['LM'] for r in rs if r['GN_conv'] and r['LM_conv'] and r['GN'] > 0 and r['LM'] > 0]); print(f"differenza relativa con segno (GN meno LM, %): mediana {np.median(rel):+.3f} | intervallo {rel.min():+.2f} .. {rel.max():+.2f}")
    for met in ('GN', 'LM'): print(f"convergenza {met}: {100*np.mean([r[met+'_conv'] for r in rs]):.1f}% | iterazioni mediane {np.median([r[met+'_it'] for r in rs]):.0f}")
    q = [r for r in rs if r['nuvola'] == '15' and r['z0'] == 1.15 and r['L'] == 0.30]
    dq = [abs(r['GN'] - r['LM']) * 1000 for r in q if r['GN_conv'] and r['LM_conv']]; print(f"nuvola 15, finestra 1.15-1.45 m (L = 0.30): {len(q)} alberi, {len(dq)} convergenti con entrambi; differenza < 1 mm in {sum(x < 1 for x in dq)}")
    for met in ('GN', 'LM'):
        e = [100 * (2 * 100 * r[met] - r['ref']) / r['ref'] for r in q if r[met + '_conv']]; print(f"  {met}: errore del diametro (solo questo tratto) contro il campo: mediana |err| {np.median(np.abs(e)):.2f}% | segno mediano {np.median(e):+.2f}%  (n = {len(e)})")
    import json; json.dump(rs, open('/home/claude/solutori/risultati_solutori.json', 'w'))
