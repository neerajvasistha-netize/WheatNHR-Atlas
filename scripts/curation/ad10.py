"""ad10 (Nov 2026): Li et al. 2010 TAG 121:941 (crown rot and FHB, DH Batavia/Ernie, SSR+DArT map).
Replaces two auto records (V1407 3B unplaced; V1408 '2D', which is the 2-49 locus cited in the Discussion, not a result of this
paper) with the four loci of Table 4. Positions from the paper's map markers (Figs 1-2 order on 3BL) at RefSeq v1.0."""
import json, re, sys
sys.path[:0] = ['scripts/curation', 'scripts/curation/sb19']
import lift, stb44
OLD = ['V1407', 'V1408']
P = lambda m, c: stb44.P(m, c)[0]
E = ' Corrected from the paper, Nov 2026 (ad10.py).'
def recs():
    a, b = P('wPt-9368', '3B'), P('wPt-1834', '3B')
    return [
     dict(name='QCr.3BL (Batavia/Ernie)', disease='Fusarium crown rot', chr='3B', start=min(a, b), end=max(a, b), peak=None, grade='B',
          marker='wPt-1151 (closest); wPt-9368-wPt-1834 (placed flanks)', pve='16.7-34.6', stat='LOD 3.3-7.6',
          note=('Major QTL for crown rot resistance on 3BL in the Batavia/Ernie doubled-haploid population, resistance allele from Ernie; detected '
                'in all six seedling trials, three with F. graminearum CS3005 (LOD 3.3-6.1, R2 16.7-28.2%) and three with F. pseudograminearum '
                'CS3096 (LOD 4.3-7.6, R2 18.7-34.6%), interval 136.4-141.3 cM, closest marker wPt-1151. Not changed when plant height was '
                'fitted as a covariate (CR05-PH LOD 6.0, R2 27.5%; CR96-PH LOD 8.8, R2 35.8%). No other CR QTL detected; none of the FHB QTL '
                'coincides with it. wPt-1151 has no physical position; the interval is the RefSeq v1.0 span of the nearest placed map '
                f'markers on either side, wPt-9368 ({a:.3f} Mb) and wPt-1834 ({b:.3f} Mb) (Figs 1-2 map order).' + E),
          evidence='Table 4: CR05-1 to CR96-3, 3BL, 136.4-141.3 cM, wPt-1151, LOD 3.3-7.6, R2 16.7-34.6%, Ernie'),
     dict(name='QFhb.3A (Batavia/Ernie)', disease='Fusarium head blight', chr='3A', start=None, end=None, peak=P('wPt-1353', '3A'), grade='C',
          marker='wPt-1353', pve='15.9', stat='LOD 3.2',
          note=('QTL for FHB resistance (point inoculation) in the Batavia/Ernie DH population, trial FHB-1 only, 4.2-5.7 cM, resistance allele '
                'from Ernie. Single marker; position is the RefSeq v1.0 position of the closest marker wPt-1353.' + E),
          evidence='Table 4: FHB-1, 3A, 4.2-5.7 cM, wPt-1353, LOD 3.2, R2 15.9%'),
     dict(name='QFhb.6A (Batavia/Ernie)', disease='Fusarium head blight', chr='6A', start=None, end=None, peak=P('wPt-6268', '6A'), grade='C',
          marker='wPt-6268', pve='17.4-20.2', stat='LOD 3.5-3.8',
          note=('QTL for FHB resistance in the Batavia/Ernie DH population, both trials (FHB-1 LOD 3.5, R2 17.4%; FHB-2 LOD 3.8, R2 20.2%), '
                '0.0-1.0 cM at the end of the linkage group, resistance allele from Ernie. Single marker; position is the RefSeq v1.0 position of '
                'the closest marker wPt-6268.' + E),
          evidence='Table 4: FHB-1 and FHB-2, 6A, 0.0-1.0 cM, wPt-6268'),
     dict(name='QFhb.5B (Batavia/Ernie)', disease='Fusarium head blight', chr='5B', start=None, end=None, peak=None, grade='D',
          marker='wPt-3661', pve='22.6-23.5', stat='LOD 4.2-4.4',
          note=('Largest QTL for FHB resistance in the Batavia/Ernie DH population, both trials (FHB-1 LOD 4.4, R2 23.5%; FHB-2 LOD 4.2, R2 '
                '22.6%), 66.8-74.4 cM, resistance allele from Ernie; near the 5B plant-height QTL of PH-2. The closest marker wPt-3661 has no '
                'physical position, so the locus is left unplaced.' + E),
          evidence='Table 4: FHB-1 and FHB-2, 5B, 66.8-74.4 cM, wPt-3661'),
    ]
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S88-')]
    new = []
    for n, f in enumerate(recs(), 1):
        r = dict(id=f'S88-{n:03d}', cls='QTL', ref='RefSeq v1.0', source='V_Li2010', curation='hand', **f)
        r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}
        new.append(r)
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
    json.dump(ops, open('/home/claude/work/ad10.json', 'w'), ensure_ascii=False)
    print(json.dumps(ops, ensure_ascii=False))
