import numpy as np, pickle, csv, time, openpyxl
from scipy.spatial import cKDTree
from scipy.ndimage import generic_filter
from scipy.optimize import least_squares
from sklearn.cluster import DBSCAN
DT = np.dtype([('X','<i4'),('Y','<i4'),('Z','<i4'),('I','<u2'),('ret','u1'),('cls','u1'),('ang','i1'),('ud','u1'),('src','<u2'),('gps','<f8'),('R','<u2'),('G','<u2'),('B','<u2')])
def carica(nome):
    raw = np.fromfile('/mnt/user-data/uploads/'+nome, dtype=DT, offset=1081)
    return raw['X']*0.001+642000.0, raw['Y']*0.001+4526000.0, raw['Z']*0.001
def dtm_rel(X, Y, Z, cs=2.0):
    ix = ((X-X.min())/cs).astype(int); iy = ((Y-Y.min())/cs).astype(int); nx, ny = ix.max()+1, iy.max()+1
    key = ix*ny+iy; o = np.argsort(key); ks = key[o]; zs = Z[o]; b = np.flatnonzero(np.diff(ks))+1; st = np.r_[0,b]; en = np.r_[b,len(ks)]
    d = np.full(nx*ny, np.nan)
    for s, e in zip(st, en): d[ks[s]] = np.quantile(zs[s:e], 0.02)
    d = d.reshape(nx, ny)
    for _ in range(5):
        m = np.isnan(d)
        if not m.any(): break
        d[m] = generic_filter(np.where(m, np.nan, d), np.nanmean, size=3, mode='nearest')[m]
    return Z - d[ix, iy]
def rileva_tronchi(X, Y, zrel):
    sel = (zrel >= 1.1) & (zrel < 1.5); P2 = np.stack([X[sel], Y[sel]], 1)
    lab = DBSCAN(eps=0.15, min_samples=15).fit(P2).labels_; tr = []
    for k in range(lab.max()+1):
        p = P2[lab == k]
        if len(p) < 100: continue
        c = np.median(p, 0); dm = 2*np.median(np.linalg.norm(p-c, axis=1))*100
        if 8 <= dm <= 70: tr.append((c[0], c[1], len(p), dm))
    return np.array(tr)
def traccia(xy, zr, c0, z_from, direzione, step=0.25, rtrack=0.35, maxstep=0.2, nmin=30, zlim=14.0):
    cen = {}; c = c0.copy(); z0 = z_from
    while 0.0 <= z0 <= zlim:
        lo, hi = (z0, z0+step) if direzione > 0 else (z0-step, z0)
        m = (zr >= lo) & (zr < hi) & (np.hypot(xy[:,0]-c[0], xy[:,1]-c[1]) < rtrack)
        if m.sum() < nmin: break
        cn = np.median(xy[m], 0); dv = cn-c; nd = np.hypot(*dv)
        if nd > maxstep: cn = c + dv/nd*maxstep
        c = cn; cen[(lo+hi)/2] = c.copy(); z0 += direzione*step
    return cen
def estrai(X, Y, zrel, c0, R):
    dx = X-c0[0]; dy = Y-c0[1]
    box = (np.abs(dx) < R+3.0) & (np.abs(dy) < R+3.0) & (zrel > -0.05) & (zrel < 14.0)
    idx = np.flatnonzero(box); xy = np.stack([dx[idx], dy[idx]], 1); zr = zrel[idx]
    cen = {**traccia(xy, zr, np.zeros(2), 1.3, -1), **traccia(xy, zr, np.zeros(2), 1.3, +1)}
    zk = np.array(sorted(cen)); ck = np.array([cen[q] for q in zk])
    if len(zk) < 3: zk = np.array([0.0, 14.0]); ck = np.zeros((2,2))
    cx = np.interp(zr, zk, ck[:,0]); cy = np.interp(zr, zk, ck[:,1])
    keep = np.hypot(xy[:,0]-cx, xy[:,1]-cy) < R
    P = np.stack([xy[keep,0], xy[keep,1], zr[keep]], 1)
    kq = np.floor(P/0.005).astype(np.int64); _, u = np.unique(kq, axis=0, return_index=True)
    return P[np.sort(u)], float(zk.max())
def hough(xy, rmin=0.05, rmax=0.40, db=0.005, off=0.10, step=0.01):
    best = (-1, None); bins = np.arange(rmin, rmax+db, db)
    for dx in np.arange(-off, off+1e-9, step):
        for dy in np.arange(-off, off+1e-9, step):
            h = np.histogram(np.hypot(xy[:,0]-dx, xy[:,1]-dy), bins=bins)[0]; hs = h[:-2]+h[1:-1]+h[2:]; k = hs.argmax()
            if hs[k] > best[0]: best = (hs[k], (dx, dy, bins[k+1]+db/2))
    return best[1]
def dbh_diretto(P, h=1.3, spess=0.2, rmax=0.6, tol=0.04):
    m = (P[:,2] >= h-spess/2) & (P[:,2] < h+spess/2) & (np.hypot(P[:,0], P[:,1]) < rmax); xy = P[m,:2]
    if len(xy) < 50: return np.nan, np.nan
    hx, hy, hr = hough(xy); sel = np.abs(np.hypot(xy[:,0]-hx, xy[:,1]-hy)-hr) < tol
    if sel.sum() >= 20:
        f = lambda p: np.hypot(xy[sel,0]-p[0], xy[sel,1]-p[1]) - p[2]
        hr2 = least_squares(f, [hx, hy, hr], loss='soft_l1', f_scale=0.01).x[2]
    else: hr2 = hr
    return 2*hr2*100, 2*hr*100
