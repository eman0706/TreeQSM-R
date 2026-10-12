import re, os, collections, json
RP = '/home/claude/Rport'
escl = set("""tools.R extract_trees.R estrai_alberi_qsm.R cover_sets_MYFIX_backup.R cover_sets_original.R treeqsm_patched.R debug_segments.R param_extreme_test.R run_new.R stress_test.R test_run.R test_run00.R test_run001.R test_run002.R test_run03.R test_exatto.R test_load.R test_semplice.R retry_weak_trees.R bootstrap_sessione.R ripristino_sessione.R""".split())
override = ['correzioni_sessione_oggi.R', 'sandbox_v5_v6_funzioni.R']          # caricati dopo il bootstrap, nell'ordine
def isola(t):
    t = re.sub(r'#.*', '', t); t = re.sub(r'"(?:[^"\\\n]|\\.)*"', '""', t); return re.sub(r"'(?:[^'\\\n]|\\.)*'", "''", t)
def funzioni(t):
    out = {}
    for m in re.finditer(r'^([A-Za-z_.][A-Za-z0-9_.]*)\s*(?:<-|=)\s*function\s*\(', t, flags=re.M):
        j = t.find('{', m.end())
        if j < 0: continue
        d = 0; k = j
        while k < len(t):
            if t[k] == '{': d += 1
            elif t[k] == '}':
                d -= 1
                if d == 0: break
            k += 1
        out[m.group(1)] = t[m.start():k+1]
    return out
defs = collections.defaultdict(list)       # nome -> [(file, corpo)]
tutti = sorted(f for f in os.listdir(RP) if f.endswith('.R') and f not in escl and not f.startswith('test_'))
ordine = [f for f in tutti if f not in override] + [f for f in override if f in tutti]
for f in ordine:
    t = isola(open(os.path.join(RP, f), encoding='utf-8', errors='ignore').read().replace('\r', ''))
    for n, b in funzioni(t).items(): defs[n].append((f, b))
eff = {n: v[-1] for n, v in defs.items()}                    # definizione effettiva: l'ultima caricata
nomi = set(eff)
edges = {n: collections.Counter(c for c in re.findall(r'\b([A-Za-z_.][A-Za-z0-9_.]*)\s*\(', b) if c in nomi and c != n) for n, (f, b) in eff.items()}
radice = 'treeqsm'
visto = {radice}; coda = [radice]
while coda:
    x = coda.pop()
    for y in edges.get(x, {}):
        if y not in visto: visto.add(y); coda.append(y)
json.dump({'definizioni': {n: [f for f, _ in v] for n, v in defs.items()}, 'effettiva': {n: f for n, (f, _) in eff.items()},
           'archi': {n: dict(c) for n, c in edges.items()}, 'raggiunte': sorted(visto)}, open('/home/claude/albero_chiamate.json', 'w'))
dup = {n: [f for f, _ in v] for n, v in defs.items() if len(v) > 1}
print(f"funzioni definite: {len(defs)} in {len(ordine)} file | raggiunte da treeqsm: {len(visto)} | non raggiunte: {len(defs) - len(visto)} | definite piu' volte: {len(dup)}")
print("\nfunzioni raggiunte da treeqsm, per file di definizione effettiva:")
per_file = collections.defaultdict(list)
for n in sorted(visto): per_file[eff[n][0]].append(n)
for f, ns in sorted(per_file.items()): print(f"  {f:<34} {len(ns):>3}: {', '.join(ns[:9])}{' ...' if len(ns) > 9 else ''}")
print("\nchiamate dirette di treeqsm:", ', '.join(sorted(edges[radice])))
print("\ndefinite piu' volte e raggiunte (prime 14):")
for n, fs in sorted(dup.items()):
    if n in visto: print(f"  {n:<28} {' | '.join(fs)}")
