"""ad47 (Nov 2026): Singh et al. 2018 IJMS 19:4054 Supplementary Table S1 (QTL with days to heading and plant height as covariates).
Adds the five covariate-only spot blotch QTL (BC 3AL, 3B, 5B, 7A; WC 2A) not in the main-text table. DArTseq tags placed with the atlas tag BLAST."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
T = ' Read from the paper and Supplementary Table S1, Nov 2026 (ad47.py).'
NEW = [('S198-001', 'BC', '3A', '28.8-32.0', ['987922', '100199979'], 'BARTAI', '2013 2.3, 2014 2.0, mean 2.3', '2.0-2.3', 'linkage group 3AL'),
       ('S198-002', 'BC', '3B', '65.6-71.3', ['1089285', '1091357'], 'BARTAI', '2015 4.1, mean 2.5', '2.5-4.1', ''),
       ('S198-003', 'BC', '5B', '142.1-142.9', ['3028824', '1108969'], 'BARTAI', '2015 3.2', '3.2', ''),
       ('S198-004', 'BC', '7A', '8.3-12.8', ['1008884', '1079791'], 'CIANO T79', '2013 3.1, mean 1.5', '1.5-3.1', 'linkage group 7A2'),
       ('S198-005', 'WC', '2A', '110.5-113.5', ['5579629', '1094923'], 'WUYA', '2013 3.0, mean 1.9', '1.9-3.0', '')]
POP = {'BC': 'BARTAI x CIANO T79', 'WC': 'WUYA x CIANO T79'}
def mk(i, pop, ch, cm, ms, src, pvs, pve, lg):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', curation='hand', name=f'QSb-{pop}-{ch} {cm} cM (Singh 2018, covariate)', disease='Spot blotch', chr=ch,
             source='V_Singh2018', marker='–'.join(ms), pve=pve, evidence=f'Table S1: {pop} {ch} {cm} {ms[0]} {ms[1]}')
    ok = [(m,) + resolve(m, ch) for m in ms]; ok = [x for x in ok if x[1]]
    base = (f'{POP[pop]} F2:7 RIL (232 lines), field spot blotch at Agua Fria, Mexico, 2013-2015; DArTseq map. Detected only with days to heading and plant height as covariates'
            f'{" (" + lg + ")" if lg else ""}; resistance allele from {src}; PVE (%) {pvs}.')
    if len(ok) == 2:
        lo, hi = sorted([ok[0][1], ok[1][1]]); r.update(start=lo, end=hi, grade='B' if hi - lo <= 100 else 'C')
        t = ' Marker positions (RefSeq v1.0): ' + '; '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in ok) + '.'
    elif ok:
        m, p, h = ok[0]; r.update(peak=p, grade='C'); t = f' Only {m} could be placed on {ch} ({p:.3f} Mb, {h}); position is that marker, not a reported peak.'
    else:
        r['grade'] = 'D'; t = f' No flanking marker could be placed on {ch}; left unplaced.'
    r['note'] = base + t + T; r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}
    return r
def build(reports):
    ids = {n[0] for n in NEW}; out = [r for r in reports if r['id'] not in ids]; new = [mk(*n) for n in NEW]; return out + new, new
M = json.load(open('data/merged_v13f.json')); M['reports'], upd = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = build(Sd['reports'])
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
ops = []
for r in upd:
    d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
    d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
    ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d}); print(r['id'], r['grade'], r.get('start'), r.get('end'), r.get('peak'), d['asm'])
json.dump(ops, open('/tmp/ops47.json', 'w'), ensure_ascii=False)
