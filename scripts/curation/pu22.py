"""Position fix (Nov 2026): Pundir et al. 2022 Sci Rep (cereal cyst nematode, DH population).
Supplementary Table S8 gives the RefSeq v1.0 position (Ensembl Plants release 52 BLAST) of the peak marker of every QTL.
- QCcn.ha-3B, -5D, -2B.1, -1D were unplaced: now placed at the printed position (grade A, new IDs).
- The other QTL already had the same position from the marker tables; grade raised B -> A (same IDs).
- QCcn.ha-4A: the paper reports no hit for Kukri_c29142_1046; stays unplaced.
Writes data/merged_v13f.json, web/index.html (seed) and /home/claude/work/pu22_db.json."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19')
import lift
S8 = 'Table S8 (RefSeq v1.0, Ensembl Plants release 52): '
NEW = {'V3186': (789.622967, 'QCcn.ha-3B Tdurum_contig15050_103 3B:789622967-789623051'),
       'V3189': (431.893815, 'QCcn.ha-5D wms0174 5D:431893815-431893836'),
       'V3190': (743.745549, 'QCcn.ha-2B.1 wsnp_RFL_Contig3917_4326857 2B:743745549-743745640'),
       'V3195': (392.554538, 'QCcn.ha-1D wsnp_RFL_Contig2036_1264133 1D:392554538-392554619')}
REGRADE = {'V3178': '2D:650325174', 'V3179': '2A:766400399', 'V3180': '5B:657711648', 'V3181': '7B:174571389', 'V3182': '6B:126145319',
           'V3183': '7D:393578460', 'V3184': '5A:565550080', 'V3185': '4B:169935788', 'V3187': '1A:547963976', 'V3188': '7A:109843375',
           'V3192': '3A:739147873', 'V3193': '6A:81849906', 'V3194': '2B:772053192', 'V3196': '1B:584156264'}
SID = lambda i: 'S82-%03d' % i

def build(reports):
    byid = {r['id']: r for r in reports}
    new, keep, changed = [], [], []
    for n, (i, (p, ev)) in enumerate(NEW.items(), 1):
        r = dict(byid[i]); r.update(id=SID(n), peak=p, start=None, end=None, grade='A', curation='hand', ref='RefSeq v1.0', evidence=S8 + ev)
        r['note'] = ('QTL for H. avenae resistance, ITMI doubled-haploid population (CIM); PVE given as R2 x 100. Peak-marker position printed in Supplementary Table S8 (RefSeq v1.0).'
                     + (' The 90K marker tables place this probe on 2D (609.49 Mb); the paper\'s BLAST hit on 2B is used, as reported.' if i == 'V3190' else '')
                     + ' Corrected from the paper, Nov 2026 (pu22.py).')
        r['asm'] = lift.lift_record(r); new.append(r)
    for r in reports:
        if r['id'] == 'V3191' and 'no hit' not in r['note']:
            r = dict(r); r['note'] += ' The paper reports no BLAST hit for Kukri_c29142_1046 (Table S8).'; changed.append((r['id'], {'note': r['note']}))
        if r['id'] in NEW or r['id'].startswith('S82-'): continue
        if r['id'] in REGRADE and r.get('grade') != 'A':
            r = dict(r); r['grade'] = 'A'; r['evidence'] = S8 + REGRADE[r['id']] + ' (peak marker ' + (r.get('marker') or '') + ')'
            changed.append((r['id'], {'grade': 'A', 'evidence': r['evidence']}))
        keep.append(r)
    return keep + new, new, list(NEW), changed

if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], new, rem, ch = build(M['reports'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</'))
    Sd['reports'], new2, rem2, ch2 = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs, 'removed': rem2, 'changed': ch2}, open('/home/claude/work/pu22_db.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['chr'], r['peak'], r['asm'].get('RefSeq v2.1'), r['note'][-160:])
    print('regraded', len(ch2))
