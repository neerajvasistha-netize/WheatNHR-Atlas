"""ad48 (Nov 2026): Slominska-Durdasiak et al. 2020 TAG ESM2 (candidate-gene primers). After removing the T3/T7 tails, exact-match in-silico PCR
(sb19/ispcr.py, RefSeq v1.0) gives single products inside TraesCS5A02G516000 (DCL1 5AL) and TraesCS4D02G267900 (Ara6 4DL)."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19'); import lift
T = ' Gene located by exact-match in-silico PCR of the Online Resource 2 primers (T3/T7 tails removed), Nov 2026 (ad48.py).'
G = {'S172-001': ('TraesCS5A02G516000', 679.655516, 679.667091, 'DCL1(5AL) primers give one 447 bp product at 679.664 Mb'),
     'S172-002': ('TraesCS4D02G267900', 439.12179, 439.124605, 'Ara6(4DL) primers give one 562 bp product at 439.123 Mb')}
def fix(r):
    g, s, e, t = G[r['id']]; r = dict(r); r.pop('peak', None)
    r.update(start=s, end=e, grade='C', gene=g, note=r['note'].replace(' No gene ID given; not placed.', f' {t}, inside {g}.').split(' Read from the paper')[0] + T,
             evidence=r['evidence'] + '; Online Resource 2')
    r['asm'] = lift.lift_record(r); return r
def build(rs):
    out, up = [], []
    for r in rs:
        if r['id'] in G: r = fix(r); up.append(r)
        out.append(r)
    return out, up
M = json.load(open('data/merged_v13f.json')); M['reports'], upd = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = build(Sd['reports'])
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
ops = []
for r in upd:
    d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}; d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
    ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
json.dump(ops, open('/tmp/ops48.json', 'w'), ensure_ascii=False); print(json.dumps(ops, ensure_ascii=False))
