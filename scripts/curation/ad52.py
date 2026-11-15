"""ad52 (Nov 2026): Zhu et al. 2014 Mol Breed ESM1 Supplementary Table 2 (linked markers of each QTL). Extra DArT/SSR markers placed in silico."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
T = ' Linked markers from Supplementary Table 2 (ESM1), Nov 2026 (ad52.py).'
U = {'S19-118': (['1002801', '990937', '1056249'], 'B', 'Supplementary Table 2 lists 1002801, 1079551, 1123330 and 1056249 as linked markers; 1079551 and 1123330 could not be placed.'),
     'S19-117': (['Xwmc128', '1089962', '1051996'], 'C', 'Supplementary Table 2 lists Xwmc128, Xgwm374, Xwmc31, Xcfd65, 1115814, 983080, 1089962 and 1051996 as linked markers; Xgwm374, Xwmc31, Xcfd65 and 983080 could not be placed, and 1115814 (28.041 Mb) is far from the others and was not used.'),
     'S19-121': (['1104599', '1075529', '1064160'], 'B', 'Supplementary Table 2 lists 1090349, 1016317, 1075529, 1064160, 1136240, 1057878 and 993225 as linked markers; 1075529 and 1064160 also place on distal 6BL, supporting the 6B assignment; the others could not be placed.')}
def fix(r):
    ms, g, t = U[r['id']]; r = dict(r); ok = [(m,) + resolve(m, r['chr']) for m in ms]
    lo, hi = min(x[1] for x in ok), max(x[1] for x in ok); r.pop('peak', None); r.update(start=lo, end=hi, grade=g)
    base = r['note'].split(' Only ')[0]
    r['note'] = base + ' ' + t + ' Marker positions (RefSeq v1.0): ' + '; '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in ok) + '.' + (f' The markers span {hi-lo:.0f} Mb: a chromosome region, not a locus.' if hi - lo > 100 else '') + T
    r['evidence'] = r.get('evidence', '') + '; Supplementary Table 2'; r['asm'] = lift.lift_record(r); return r
def build(rs):
    out, up = [], []
    for r in rs:
        if r['id'] in U: r = fix(r); up.append(r)
        out.append(r)
    return out, up
M = json.load(open('data/merged_v13f.json')); M['reports'], upd = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = build(Sd['reports'])
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
ops = []
for r in upd:
    d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}; d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
    ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'ifVersion': 1, 'data': d})
json.dump(ops, open('/tmp/ops52.json', 'w'), ensure_ascii=False); print(json.dumps(ops, ensure_ascii=False))
