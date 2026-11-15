"""ad41 (Nov 2026): Ferjaoui et al. 2022 BMC Genomics 23:372, Additional Table 4 (Svevo RefSeq v1.0 positions of QTL flanking tags).
Svevo positions used as approximate RefSeq v2.1 positions (flag as_reported), as for other Svevo/Zavitan-placed records."""
import json, re
SV = ' Additional Table 4 gives Svevo (RefSeq v1.0) positions from BLAST of the tag sequences; used here as approximate CS v2.1 positions (grade C). Updated Nov 2026 (ad41.py).'
U = {
 'S44-022': dict(note_add=' Additional Table 4: Tag_1694925 at 489.686 Mb on Svevo 1A; Tag_1127081 had no Svevo hit. The CS v1 placement of Tag_1127081 is kept. Updated Nov 2026 (ad41.py).'),
 'S44-023': dict(ref='RefSeq v2.1', start=139.975620, end=145.652386, grade='C', note='Durum Agili39 x Khiar RIL; DArTseq flanking tags Tag_2277193 (52.8 cM, 139.976 Mb) and Tag_100009953 (60.3 cM, 145.652 Mb).' + SV),
 'S44-024': dict(ref='RefSeq v2.1', start=559.062634, end=788.999946, grade='C', marker='Tag_1056626, Tag_111757', note='Durum Agili39 x Khiar RIL; DArTseq flanking tags Tag_1056626 (73.9 cM, 788.999 Mb) and Tag_111757 (77.0 cM, 559.063 Mb). The two Svevo hits are 230 Mb apart although the tags are 3 cM apart, so one hit is probably paralogous; the record spans both (a chromosome region, not a locus).' + SV),
 'S44-025': dict(ref='RefSeq v2.1', peak=783.725692, grade='C', note='Durum Agili39 x Khiar RIL; DArTseq flanking tags Tag_100031118 (111.9 cM, no Svevo hit) and Tag_3027184 (123.2 cM, 783.726 Mb); position is Tag_3027184.' + SV),
}
def apply(R):
    out = []
    for r in R:
        if r['id'] in U:
            u = dict(U[r['id']])
            if 'note_add' in u: r['note'] = r['note'] + u.pop('note_add')
            else:
                for k in ('start', 'end', 'peak'): r.pop(k, None)
                r.update(u)
                r['asm'] = {'RefSeq v2.1': {'chr': r['chr'], 'flag': 'as_reported', **({'start': r['start'], 'end': r['end']} if 'start' in r else {'peak': r['peak']})}}
            out.append(r)
    return out
M = json.load(open('data/merged_v13f.json')); ch = apply(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
Sd = json.loads(m.group(2).replace('<\\/', '</')); apply(Sd['reports'])
open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
ops = []
for r in ch:
    d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
    d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
    for k in ('start', 'end', 'peak'):
        if k not in d: d[k] = {'__delete__': True}
    ops.append({'op': 'update', 'path': 'reports/' + r['id'], 'ifVersion': 1, 'data': d})
json.dump(ops, open('/tmp/ops41.json', 'w'), ensure_ascii=False)
print(json.dumps(ops, ensure_ascii=False))
