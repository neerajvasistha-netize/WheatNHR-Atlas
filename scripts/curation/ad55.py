"""ad55 (Nov 2026, user rule): when a preprint was later published and both are in the atlas, remove the preprint.
Removes every record of nine preprint sources (their journal versions stay) and the preprint source keys, incl. the empty V_Ferjaoui2021."""
import json, re, datetime
PRE = {'Aouini2025bR': 'A_STB_039', 'V_Gardner2025': 'Z_1123', 'A_spot_blotch_026': 'Tomar2021FG', 'V_Ren2023': 'Z_1350', 'V_Concepcion2024b': 'V_Concepcion2024',
       'GlobalWW2023': 'A_SNB_001', 'Y_209': 'V_Liu2021b', 'Y_126': 'Z_1463', 'V_Chaudhary2024': 'Chaudhary2025BMC', 'V_Ferjaoui2021': 'Ferjaoui2022BMC'}
VER = {'R052': 2, **{f'S70-0{n}': 1 for n in range(13, 25)}}
def build(D):
    gone = [r for r in D['reports'] if r.get('source') in PRE]
    D['reports'] = [r for r in D['reports'] if r.get('source') not in PRE]
    for k in PRE: D['sources'].pop(k, None)
    return gone
M = json.load(open('data/merged_v13f.json')); gone = build(M); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); g2 = build(Sd)
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
now = datetime.datetime.utcnow().isoformat() + 'Z'
ops = []
for r in gone:
    o = {'op': 'set', 'path': 'reports/' + r['id'], 'data': {'deleted': True, 'deleted_by': 'curation ad55 (preprint superseded by ' + PRE[r['source']] + ')', 'deleted_at': now, 'name': r.get('name', '')}}
    if r['id'] in VER: o['ifVersion'] = VER[r['id']]
    ops.append(o)
json.dump(ops, open('/tmp/ops55.json', 'w'), ensure_ascii=False)
json.dump([[r['id'], r['source'], PRE[r['source']], r['chr'], r.get('name', '')] for r in gone], open('docs/ad55_preprint_removed.json', 'w'), ensure_ascii=False, indent=0)
import collections; print(len(gone), len(g2), collections.Counter(r['source'] for r in gone))
