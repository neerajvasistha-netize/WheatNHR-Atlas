"""Exact-match in-silico PCR of SSR/STS primer pairs on IWGSC RefSeq v1.0 (Ensembl Plants chromosomes).
Usage: python ispcr.py primers.json markers.txt out.json [chrom ...]
Reports every product (F...rc(R) or R...rc(F)) between 40 and 2000 bp."""
import sys, json, gzip, ahocorasick, os
P = json.load(open(sys.argv[1]))
want = [l.strip() for l in open(sys.argv[2]) if l.strip()]
out_path = sys.argv[3]
chroms = sys.argv[4:] or '1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D'.split()
comp = bytes.maketrans(b'ACGTN', b'TGCAN')
rc = lambda s: s.translate(comp)[::-1]
A = ahocorasick.Automaton()
pats = {}
for m in want:
    if m not in P: continue
    F, R = P[m][0].upper().encode(), P[m][1].upper().encode()
    for tag, s in (('F', F), ('Frc', rc(F)), ('R', R), ('Rrc', rc(R))):
        pats.setdefault(s.decode(), []).append((m, tag, len(s)))
for s, v in pats.items(): A.add_word(s, v)
A.make_automaton()
res = json.load(open(out_path)) if os.path.exists(out_path) else {}
for c in chroms:
    fa = '/home/claude/genome/%s.fa.gz' % c
    seq = b''.join(l.strip() for l in gzip.open(fa) if not l.startswith(b'>')).upper().decode()
    hits = {}
    for end, vals in A.iter(seq):
        for m, tag, L in vals:
            hits.setdefault(m, {}).setdefault(tag, []).append(end - L + 1)
    for m, h in hits.items():
        prods = []
        # plus strand product: F at i, rc(R) ending at j
        for i in h.get('F', []):
            for j in h.get('Rrc', []):
                if 40 <= j - i + len(P[m][1]) <= 2000: prods.append((i, j + len(P[m][1]), '+'))
        for i in h.get('R', []):
            for j in h.get('Frc', []):
                if 40 <= j - i + len(P[m][0]) <= 2000: prods.append((i, j + len(P[m][0]), '-'))
        for a, b, st in prods:
            res.setdefault(m, []).append([c, round(a / 1e6, 4), b - a, st])
    print(c, 'done', flush=True)
    json.dump(res, open(out_path, 'w'))
