"""Apply r98_fix.py corrections to data/merged_v13f.json and the seed in web/index.html, and write the
matching live-database merge documents to work/r98_db.json. Run after sb26_apply.py.
Each corrected record gets a new id (R98-nnn) and the old id is deleted: the page keeps precomputed
assembly conversions per seed id, so editing a seed record in place would leave its RefSeq v2.1 position
marked "conversion pending"."""
import json, re, sys
sys.path.insert(0, 'scripts/curation'); sys.path.insert(0, 'scripts/curation/sb19')
import r98_fix, lift

def fix(r, c):
    old = dict(r)
    for k in ('start', 'end', 'peak'): r.pop(k, None)
    r.update({k: v for k, v in c['new'].items()})
    extra = (' The earlier value was the %s position, converted.' % r.get('ref')) if r.get('ref') != 'RefSeq v1.0' else ''
    r['ref'] = 'RefSeq v1.0'
    r['note'] = (str(r.get('note') or '') + ' ' + c['note'] + extra).strip()
    r['curation'] = 'hand'
    for k in ('marker', 'evidence', 'note', 'name'):
        if isinstance(r.get(k), str): r[k] = r[k].replace('\u2014', '\u2013')
    r['asm'] = lift.lift_record(r)
    gone = [k for k in old if k not in r]
    return r, gone

if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    C = r98_fix.run(M['reports'])
    NID = {c['id']: 'R98-%03d' % i for i, c in enumerate(C, 1)}
    CC = {c['id']: c for c in C}
    def redo(reports):
        out, new = [], []
        for r in reports:
            if r['id'] in CC:
                r2, _ = fix(dict(r), CC[r['id']])
                r2['note'] += ' Replaces atlas record %s.' % r['id']
                r2['id'] = NID[r['id']]
                new.append(r2)
            else: out.append(r)
        return out + new, new
    M['reports'], new = redo(M['reports'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    S = json.loads(m.group(2).replace('<\\/', '</'))
    S['reports'], new2 = redo(S['reports'])
    H = H[:m.start(2)] + json.dumps(S, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):]
    open('web/index.html', 'w').write(H)
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs, 'removed': list(NID)}, open('/home/claude/work/r98_db.json', 'w'), ensure_ascii=False)
    json.dump(NID, open('scripts/curation/r98_id_map.json', 'w'), indent=0)
    print('master', len(new), 'seed', len(new2))
