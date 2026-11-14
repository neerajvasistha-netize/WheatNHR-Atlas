"""ad22 (Nov 2026): Jia et al. 2005 Euphytica 146:183 (FHB, Wangshuibai x Alondra's DH, Jianyang, natural infection 2001-2003).
Replaces 7 unplaced auto records (V1457-V1463) with the CIM QTL of Table 2 (detections of the same interval merged) and the
epistasis-only 6A QTL of Table 3. 1D and 7D were single-marker hits not significant in CIM and are not kept. No positions printed;
SSR markers placed in silico."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
OLD = ['V%04d' % i for i in range(1457, 1464)]
Q = [
 ('2D', '2D', ['Xgwm261', 'Xgwm484'], 'Wangshuibai', '2001 LOD 2.15 (R2 8.2%), 2002 2.90 (8.0%), 2001-2003 mean 2.46 (10.6%)', 10.6, ''),
 ('3A', '3A', ['Xgwm369', 'Xbarc045'], 'Wangshuibai', '2001 LOD 2.59 (R2 12.3%)', 12.3, ' Also in significant additive x additive epistasis with 4B and 6A (Table 3).'),
 ('3BS', '3B', ['Xgwm533', 'Xgwm493'], 'Wangshuibai', '2001 LOD 2.80 (R2 14.4%), 2002 2.33 (12.7%), 2003 2.55 (13.0%), 2001-2003 mean 2.54 (11.0%)', 14.4,
  ' The only QTL stable in all three years; Xgwm533 had the largest marker-class effect and is proposed for marker-assisted selection.'),
 ('3D', '3D', ['Xgwm645', 'Xgwm383'], 'Wangshuibai', '2003 LOD 2.66 (R2 11.6%)', 11.6, ''),
 ('4B', '4B', ['Xgwm368', 'Xgwm149'], 'Wangshuibai', '2001 LOD 2.70 (R2 10.5%), 2001-2003 mean 2.47 (9.9%)', 10.5, ''),
 ('5A', '5A', ['Xgwm129', 'Xgwm156'], 'Wangshuibai', '2003 LOD 3.54 (R2 14.5%)', 14.5, ''),
 ('5B.1', '5B', ['Xgwm335', 'Xgwm371'], 'Wangshuibai', '2001 LOD 3.15 (R2 15.7%), 2003 2.18 (12.2%), 2001-2003 mean 2.05 (10.8%)', 15.7,
  ' Largest additive effect in 2001; Xgwm335 proposed for marker-assisted selection.'),
 ('5B.2', '5B', ['Xgwm443', 'Xbarc32'], 'Wangshuibai', '2002 LOD 2.70 (R2 13.1%), 2001-2003 mean 2.86 (13.3%)', 13.3, ''),
 ('5D', '5D', ['Xgwm190', 'Xgwm358'], "Alondra's", '2002 LOD 2.60 (R2 11.6%)', 11.6, " The Wangshuibai allele increased severity, so resistance comes from Alondra's."),
 ('6B', '6B', ['Xgwm133', 'Xgwm191'], 'Wangshuibai', '2001 LOD 2.31 (R2 8.8%)', 8.8, ''),
 ('7A', '7A', ['Xgwm276', 'Xgwm282'], 'Wangshuibai', '2001 LOD 2.84 (R2 12.0%), 2003 2.40 (7.7%), 2001-2003 mean 2.75 (12.6%)', 12.6, ''),
 ('7B.1', '7B', ['Xgwm297', 'Xgwm644'], 'Wangshuibai', '2002 LOD 2.67 (R2 7.9%)', 7.9, ''),
 ('7B.2', '7B', ['Xgwm146', 'Xgwm611'], "Alondra's", '2003 LOD 2.35 (R2 11.4%)', 11.4, " The Wangshuibai allele increased severity, so resistance comes from Alondra's."),
 ('6A (epistatic)', '6A', ['Xgwm169'], 'Wangshuibai', 'epistasis with 3A only: LOD 3.56, AA -2.78 scabbed spikelets, R2 12.0%', 12.0,
  ' No additive effect; detected only as additive x additive epistasis with the 3A QTL (Table 3), 4 cM from Xgwm169.'),
]
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S99-')]; new = []
    for n, (nm, ch, ms, par, det, mx, ex) in enumerate(Q, 1):
        pl = [(m,) + resolve(m, ch) for m in ms]; pl = [x for x in pl if x[1]]
        r = dict(id=f'S99-{n:03d}', name=f'QFhb {nm} (Wangshuibai x Alondra\'s)', cls='QTL', disease='Fusarium head blight', chr=ch, ref='RefSeq v1.0',
                 start=None, end=None, peak=None, marker='-'.join(ms), pve=f'up to {mx}', stat='', source='V_Jia2005', curation='hand')
        if len(pl) == 2:
            r['start'], r['end'] = min(x[1] for x in pl), max(x[1] for x in pl); r['grade'] = 'B' if r['end'] - r['start'] <= 100 else 'C'
            place = ' Placed in silico: ' + ', '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in pl) + '.'
            if r['grade'] == 'C': place += f' Interval width {r["end"] - r["start"]:.0f} Mb: a chromosome region, not a locus.'
        elif pl:
            m, p, h = pl[0]; r['peak'] = p; r['grade'] = 'C'
            place = f' Only {m} could be placed ({p:.3f} Mb, {h}); position is that marker.'
            
        else:
            r['grade'] = 'D'; place = ' Neither flanking SSR could be placed in silico, so the locus is unplaced.'
        r['note'] = (f'QTL for FHB resistance (scabbed spikelets per spike) in the Wangshuibai x Alondra\'s DH population (Jianyang, Fujian, natural infection 2001-2003; QTLMapper CIM), resistance allele from {par}. '
                     f'Detections: {det}.' + ex + place + ' No physical position printed. Corrected from the paper, Nov 2026 (ad22.py).')
        if 'epistasis' in det: r['note'] = r['note'].replace(', resistance allele from Wangshuibai', '')
        r['evidence'] = ('Table 3' if 'epistasis' in det else 'Table 2') + f': {"-".join(ms)} {ch}; {det}'
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
