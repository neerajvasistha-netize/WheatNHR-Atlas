"""ad17 (Nov 2026): in-silico placement of unplaced records that name a marker. Of the 762 records without a v2.1 position,
261 name markers; the resolver (90K/55K/660K tables, DArT and DArTseq tag BLAST, exact-match in-silico PCR of SSR primers,
marker names that encode a position) finds a position on the record's own chromosome for 10. Seven of those had already been
reviewed and left unplaced on purpose (arm conflict, or flanks far apart: Z0154, Z0158, Z0419, Z0421, Y0050, V3281, S19-081).
The three Navathe et al. (2023) DArTseq markers are placed here (grade C, single marker)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
OLD = ['R0658', 'R0659', 'R0663']
def build(reports):
    byid = {r['id']: r for r in reports}; keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S94-')]; new = []
    for n, i in enumerate(OLD, 1):
        r = dict(byid[i]); m = r['marker']; p, h = resolve(m, r['chr'])
        r.update(id=f'S94-{n:03d}', peak=round(p, 4), start=None, end=None, grade='C', curation='hand', ref='RefSeq v1.0')
        r['note'] = (f'Marker-trait association for Septoria nodorum blotch, Navathe et al. (2023); the marker is named chromosome:DArTseq-ID. '
                     f'No position is printed; placed in silico at the RefSeq v1.0 position of DArTseq tag {m.split(":")[1]} ({h}). '
                     'Single marker. Placed in silico, Nov 2026 (ad17.py).')
        r['asm'] = lift.lift_record(r); new.append(r)
    return keep + new, new, OLD
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], new, rem = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], new2, rem2 = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    ops = [{'op': 'set', 'path': 'reports/' + i, 'data': {'deleted': True}} for i in rem2]
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    print(json.dumps(ops, ensure_ascii=False))
