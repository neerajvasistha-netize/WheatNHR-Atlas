"""ad16 (Nov 2026): Thambugala et al. 2020 TAG 133:2775 (native FHB; Cutler/AC Barrie RIL and AC Barrie/Reeder DH).
Replaces the 21 unplaced auto records (V0332-V0352) with the 21 FHB QTL of Tables 5 and 6. The paper gives cM only (no CI); each
QTL is placed from its population's linkage map (Tables S5, S11): 90K map markers within the range of detected peak positions
near the strongest peak (+/- 1 cM) are given their RefSeq v1.0 positions (marker resolver), outliers > 50 Mb from the median dropped. Grade B, C if > 100 Mb."""
import json, re, sys, statistics
sys.path[:0] = ['scripts/curation', 'scripts/curation/sb19']
import lift, stb44
OLD = ['V%04d' % i for i in range(332, 353)]
T = json.load(open('/tmp/t3631/t56.json')); MAPS = json.load(open('/tmp/t3631/maps.json'))
POP = {'CB': 'Cutler/AC Barrie RIL population (Table 5)', 'BR': 'AC Barrie/Reeder DH population (Table 6)'}
TR = {'INC': 'incidence', 'SEV': 'severity', 'VRI': 'visual rating index', 'DON': 'DON'}
EX = {'QFhb.mcb-3B': ' Near Fhb1 (3BS); Fhb1 candidate genes TaHRC and TaPFT were characterised in the parents.', 'QFhb.mcb-4D': ' Co-located with the plant-height QTL QHt.mcb-4D at Rht-D1.',
      'QFhb.mcb-2D.1': ' Near Ppd-D1, co-located with anthesis-date and height QTL.', 'QFhb.mbr-5A': ' Near Vrn-A1 (587.4 Mb) per the paper.', 'QFhb.mcb-5A.1': ' Fhb5 region of 5AS per the paper.'}
_cache = {}
def P(m, c):
    if (m, c) not in _cache: _cache[(m, c)] = stb44.P(m, c)[0]
    return _cache[(m, c)]
def recs():
    out = []
    for key, rows in T.items():
        pop, nm = key.split('|')
        if 'Fhb' not in nm: continue
        lg = rows[0][1]; ch = lg[:2]
        det = []
        for loc, _, nums, _ in rows:
            if len(nums) == 8: det += [(nums[0], nums[1], nums[2], loc, 'IM'), (nums[4], nums[5], nums[6], loc, 'ICIM')]
            elif len(nums) == 4:
                det.append((nums[0], nums[1], nums[2], loc, '?'))
        best = max(det, key=lambda d: d[1]); cms = [d[0] for d in det if abs(d[0] - best[0]) <= 5]
        for w in (1, 2, 3, 5):
            lo, hi = min(cms) - w, max(cms) + w
            mk = [(m, c) for m, l, c in MAPS[pop] if l == lg and lo <= c <= hi]
            pl = [(m, c, P(m, ch)) for m, c in mk]; pl = [x for x in pl if x[2]]
            if pl: break
        if not pl:
            out.append((pop, nm, ch, lg, None, None, None, best, rows, lo, hi, 0, 0)); continue
        med = statistics.median(x[2] for x in pl); core = [x for x in pl if abs(x[2] - med) <= 50]
        pk = min(core, key=lambda x: abs(x[1] - best[0]))
        out.append((pop, nm, ch, lg, min(x[2] for x in core), max(x[2] for x in core), pk, best, rows, lo, hi, len(core), len(pl) - len(core)))
    return out
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S93-')]
    new = []
    for n, (pop, nm, ch, lg, a, b, pk, best, rows, lo, hi, k, drop) in enumerate(recs(), 1):
        traits = sorted({TR.get(r[0].split('_')[0], r[0]) for r in rows}); envs = len(rows)
        meth = 'ICIM' if best[4] == 'ICIM' else 'IM'
        r = dict(id=f'S93-{n:03d}', name=nm, cls='QTL', disease='Fusarium head blight', chr=ch, ref='RefSeq v1.0', start=a, end=b, peak=None,
                 marker=(f'{pk[0]} (nearest placed map marker to the strongest peak)' if pk else ''), pve=str(best[2]), stat=f'LOD {best[1]} ({best[3]})',
                 source='V_Thambugala2020', curation='hand')
        if a is None:
            r['grade'] = 'D'; place = ' No map marker in the detected range could be placed physically, so the locus is unplaced.'
        else:
            r['grade'] = 'B' if b - a <= 100 else 'C'
            place = (f' Physical interval: RefSeq v1.0 span of the {k} placed map marker{"s" if k > 1 else ""} within {lo:.1f}-{hi:.1f} cM on '
                     f'{lg} (range of the reported peaks within 5 cM of the strongest one, widened until a placed marker was found)' + (f', {drop} outlying marker(s) dropped' if drop else '') + '.'
                     + (f' Interval width {b - a:.0f} Mb: a chromosome region, not a locus.' if r['grade'] == 'C' else ''))
        r['note'] = (f'Native FHB resistance QTL in the {POP[pop]}, field FHB nurseries in Canada; traits: {", ".join(traits)}; '
                     f'{envs} trait-environment detections. Strongest: {best[3]}, LOD {best[1]}, PVE {best[2]}% at {best[0]} cM on {lg}.'
                     + EX.get(nm, '') + place + ' The paper gives map positions only. Corrected from the paper, Nov 2026 (ad16.py).')
        r['evidence'] = f'Table {"5" if pop == "CB" else "6"}: {nm}, {lg}, {envs} trait-environment detections'
        r['asm'] = lift.lift_record(r) if a is not None else {'_note': 'no physical position in source'}
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
    json.dump(ops, open('/home/claude/work/ad16.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['chr'], r['start'], r['end'], r['grade'], r['stat'], r['marker'][:25], (r['asm'].get('RefSeq v2.1') or {}).get('flag'))
