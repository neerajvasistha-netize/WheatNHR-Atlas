"""ad21 (Nov 2026): Steiner et al. 2004 TAG 109:215 (FHB, Frontana x Remus DH). Replaces 7 unplaced auto records (V1383-V1389)
with the FHB QTL of Tables 4 and 6 (severity, incidence, spread; regions shared by both tables merged). No positions printed;
SSR markers placed in silico; AFLP-only QTL stay unplaced (grade D). Plant height, anthesis and spike-compactness QTL not added;
the spread QTL on an unassigned linkage group is left out (no chromosome)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
OLD = ['V%04d' % i for i in range(1383, 1390)]
Q = [  # name, chr, markers, resistance source, detections (LOD, VE%), max VE, extra
 ('1B', '1B', ['Xs12m25_14', 'Xs24m17_2'], 'Remus', 'severity 3-year LOD 2.1 (5.5%), 2002 LOD 2.5 (6.5%); incidence 2002 LOD 2.8 (7.2%)', 7.2, ''),
 ('2A', '2A', ['Xs13m26_4', 'Xs12m15_7'], 'Remus', 'severity 3-year LOD 2.9 (7.9%), 1999 3.1 (8.4%), 2001 3.0 (8.2%); spread 2-year LOD 2.5 (6.5%), 2002 2.4 (6.1%)', 8.4, ''),
 ('2B.1', '2B', ['Xs13m25_8', 'Xs24m15_6'], 'Frontana', 'severity 3-year LOD 2.3 (6.1%), 1999 3.2 (8.7%); incidence 2002 LOD 2.3 (5.9%); spread 2002 LOD 3.0 (7.9%)', 8.7, ''),
 ('2B.2', '2B', ['Xs25m15_2', 'Xs24m19_11'], 'Frontana', 'spread 2-year LOD 3.0 (8.1%), 2001 2.8 (7.4%), 2002 3.2 (8.7%)', 8.7, ' Consistent spread (type II) QTL in both years.'),
 ('3A', '3A', ['Xdupw227', 'Xgwm720'], 'Frontana', 'severity 3-year LOD 6.6 (16.2%), 1999 4.8 (12.0%), 2001 4.9 (12.1%), 2002 5.2 (12.9%); incidence 2002 LOD 4.6 (11.5%)', 16.2,
  ' The major QTL, described as near the centromere, most likely in the Xgwm720-Xdupw227 interval; not associated with spread; stable in all years. Lines with Frontana alleles at both 3A and 5A had about 40% less severity.'),
 ('3B', '3B', ['Xs23m15_3'], 'Remus', 'severity 2001 LOD 2.8 (7.2%)', 7.2, ''),
 ('4B', '4B', ['Xs13m25_9'], 'Frontana', 'incidence 2002 LOD 2.8 (7.4%)', 7.4, ' Overlaps a plant-height QTL.'),
 ('5A', '5A', ['Xgwm129', 'Xbarc197'], 'Frontana', 'severity 3-year LOD 3.5 (8.8%), 1999 3.6 (9%), 2001 2.1 (5.4%), 2002 3.1 (7.8%)', 9.0,
  ' Stable over years; composite interval mapping suggests a single QTL linked to Xgwm129; co-located with a plant-height QTL at Xbarc197.'),
 ('5AL', '5A', ['Xs23m20_8'], 'Frontana', 'severity 1999 LOD 2.5 (6.5%)', 6.5, ''),
 ('5A (spread)', '5A', ['Xs12m25_9'], 'Remus', 'spread 2001 LOD 2.7 (7.3%)', 7.3, ''),
 ('6A', '6A', ['Xgwm169'], 'Remus', 'spread 2-year LOD 2.9 (7.6%), 2001 2.9 (7.5%)', 7.6, ''),
 ('6B', '6B', ['Xs23m14_4'], 'Frontana', 'severity 3-year LOD 2.5 (6.7%), 1999 2.5 (6.6%), 2002 2.7 (5.8%); incidence 2002 LOD 2.0 (5.5%)', 6.7, ''),
]
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S98-')]; new = []
    for n, (nm, ch, ms, par, det, mx, ex) in enumerate(Q, 1):
        pl = [(m,) + resolve(m, ch) for m in ms]; pl = [x for x in pl if x[1]]
        r = dict(id=f'S98-{n:03d}', name=f'QFhb {nm} (Frontana x Remus)', cls='QTL', disease='Fusarium head blight', chr=ch, ref='RefSeq v1.0',
                 start=None, end=None, peak=None, marker='-'.join(ms), pve=f'up to {mx}', stat='', source='V_Steiner2004', curation='hand')
        if len(pl) == 2:
            r['start'], r['end'] = min(x[1] for x in pl), max(x[1] for x in pl); r['grade'] = 'B' if r['end'] - r['start'] <= 100 else 'C'
            place = ' Placed in silico: ' + ', '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in pl) + '.'
            if r['grade'] == 'C': place += f' Interval width {r["end"] - r["start"]:.0f} Mb: a chromosome region, not a locus.'
        elif pl:
            m, p, h = pl[0]; r['peak'] = p; r['grade'] = 'C'
            place = f' Only {m} could be placed ({p:.3f} Mb, {h}); position is that marker.'
            if nm == '3A': place += ' This puts the locus on 3AL, not at the centromere as the paper describes; treat the position with care.'
        else:
            r['grade'] = 'D'; place = ' The map markers are AFLPs with no physical position, so the locus is unplaced.'
        r['note'] = (f'QTL for FHB resistance in the Frontana x Remus DH population (Tulln, Austria; spray and single-spikelet inoculation), resistance allele from {par}. '
                     f'Detections: {det}.' + ex + place + ' No physical position printed. Corrected from the paper, Nov 2026 (ad21.py).')
        r['evidence'] = ('Table 6' if 'spread' in det and 'severity' not in det else 'Table 4') + f': {"-".join(ms)} {ch}; {det}'
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
    print(json.dumps(ops, ensure_ascii=False, separators=(',', ':')))
