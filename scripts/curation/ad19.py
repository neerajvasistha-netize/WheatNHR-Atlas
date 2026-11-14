"""ad19 (Nov 2026): He et al. 2016 PLoS ONE 11:e0158052 (FHB, Soru#1 x Naxos RIL). Replaces 7 unplaced auto records (V1552-V1558)
with the 9 disease QTL of Tables 4-6 (spray/spawn, with DH/PH covariates, point). Positions are not printed; flanking 90K markers
are placed with the marker resolver (paper abbreviations expanded: Ex_ = Excalibur_/wsnp_Ex_, w_ = wsnp_, D_c = D_contig, GENEU = GENE-);
Rht-D1 and Vrn-A1 at their gene positions."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
OLD = ['V%04d' % i for i in range(1552, 1559)]
GENE = {'Rht-D1': (18.781, 'Rht-D1 gene, TraesCS4D01G040400'), 'Vrn-A1': (587.416, 'Vrn-A1, position given by Thambugala et al. 2020')}
def alias(m):
    c = [m]
    if m.startswith('Ex_c'): c += ['Excalibur_' + m[3:], 'wsnp_' + m]
    if m.startswith('w_Ku_'): c += ['wsnp_' + m[2:]]
    if m.startswith('wExc'): c += ['wsnp_Ex_c' + m[4:]]
    if m.startswith('wC11_r_'): c += ['wsnp_CAP11_rep_' + m[7:]]
    if m.startswith('D_c'): c += ['D_contig' + m[3:]]
    if m.startswith('GENEU'): c += ['GENE-' + m[5:]]
    return c
def P(m, ch):
    if m in GENE: return GENE[m] + (m,)
    for a in alias(m):
        p, h = resolve(a, ch)
        if p: return p, h, a
    return None, '', m
KASP = [{'m': 'Kukri_c36639_186 (IWB44589, KASP)', 'F': 'TCTGCTAGAGGAGACATCACG (FAM) / TCTGCTAGAGGAGACATCACA (VIC)', 'R': 'TGTGCTCTGACTGACTTGATT', 'type': 'KASP', 'src': 'He et al. 2016 PLoS ONE (S3 Table)'},
        {'m': 'Excalibur_c7282_512 (IWB28643, KASP)', 'F': 'CGGTGACAAGGCTCTCTTCG (FAM) / CGGTGACAAGGCTCTCTTCA (VIC)', 'R': 'TGTTGCTGGAGATGTTCTCGT', 'type': 'KASP', 'src': 'He et al. 2016 PLoS ONE (S3 Table)'}]
# name arm, chr, flank sets [(a,b,label)], detections, max R2, traits, source parent, extra
Q = [
 ('2DLc', '2D', [('Kukri_c36639_186', 'Ex_c7282_512', 'Table 4, 55.3-58.2 cM'), ('w_Ku_c8712_14751858', 'GENEU0808_728', 'Table 5 (covariates), 53.2-57.8 cM')],
  'FHB index spray 2011 18.1, 2012 15.8, mean 15.9; FDK 7.4; DON 2011 20.1, 2012 13.6, mean 16.5; spawn severity 2011 15, 2012 22.4, mean 18.2; still 9.7-17.1 with DH and PH as covariates', 22.4,
  'FHB index, FDK, DON and spawn severity', 'Soru#1', ' The major QTL, stable across experiments and after correction for heading date and plant height; two flanking KASP markers developed (S3 Table).'),
 ('2DL', '2D', [('BS00021881_51', 'GENEU3588_374', 'Table 4, 162.6-163.4 cM'), ('BS00021881_51', 'Kukri_c31995_1948', 'Table 5, 162.6-164.4 cM')],
  'FHB index 2011 6.3, 2012 8.6, mean 4.4; DON 2011 9.4, 2012 7.5, mean 6.1; spawn 2012 5.3, mean 3.2', 9.4, 'FHB index, DON and spawn severity', 'Naxos', ' Associated with anther extrusion.'),
 ('4AL', '4A', [('Ex_c11968_204', 'RAC875_c35979_263', 'Tables 4 and 6, 101.2-102.6 cM')],
  'FHB index 2011 7.1; FDK 5.4; spawn 2011 6.8, 2012 4, mean 5.9; point 2011 6.6', 7.1, 'FHB index, FDK, spawn and point severity', 'Soru#1', ' Associated with days to heading.'),
 ('4DS', '4D', [('Rht-D1', 'D_c56766_278', 'Table 4, 0.0-15.9 cM')],
  'FDK 19; with covariates FHB index 4.6, FDK 3.2, DON 5.4', 19.0, 'FDK', 'Naxos', ' At Rht-D1; associated with plant height.'),
 ('5AL', '5A', [('Vrn-A1', 'Ex_c31769_793', 'Table 4, 167.2-171.2 cM'), ('Vrn-A1', 'Ex_c7729_144', 'Table 6 (point), 167.2-167.6 cM'), ('Ku_c12469_983', 'IAAV4799', 'Table 5 (covariates), 173.4-175.2 cM')],
  'FHB index 2011 36.8, 2012 11.2, mean 22.8; FDK 5.7; DON 2011 17.5, 2012 14.1, mean 17.6; spawn 2011 8, 2012 7.5, mean 7.8; point 2013 12.3, mean 10.9; with covariates 4.3-6.8', 36.8,
  'FHB index, FDK, DON, spawn and point severity', 'Soru#1', ' At Vrn-A1; associated with days to heading and plant height; the effect falls sharply with DH and PH as covariates.'),
 ('2AS', '2A', [('Ex_c18324_390', 'BS00022331_51', 'Table 5 (covariates), 7.2-8.6 cM')],
  'with covariates: FHB index mean 4.5; spawn 2011 5.9, 2012 4.3, mean 5.2', 5.9, 'FHB index and spawn severity (with DH and PH covariates)', 'Naxos', ''),
 ('5DL', '5D', [('Ra_c27043_437', 'IAAV2323', 'Table 5 (covariates), 5.8-10.6 cM')],
  'with covariates: DON 2011 4.7, 2012 2.7; spawn 2011 6.4, mean 3.8', 6.4, 'DON and spawn severity (with DH and PH covariates)', 'Soru#1', ''),
 ('2DS', '2D', [('D_F1BEJMU02GB94Z_188', 'wExc1408_2704613', 'Table 6 (point), 16.0-22.0 cM')], 'point 2011 13.5', 13.5, 'FHB severity (point inoculation)', 'Naxos', ''),
 ('3AS', '3A', [('Kukri_c96747_274', 'wC11_r_c4157_1965583', 'Table 6 (point), 12.5-18.0 cM')], 'point 2013 9.9', 9.9, 'FHB severity (point inoculation)', 'Soru#1', ''),
]
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S96-')]; new = []
    for n, (arm, ch, sets, det, mx, tr, par, ex) in enumerate(Q, 1):
        a, b, lab = sets[0]; pa = P(a, ch); pb = P(b, ch)
        r = dict(id=f'S96-{n:03d}', name=f'QFhb {arm} (Soru#1 x Naxos)', cls='QTL', disease='Fusarium head blight', chr=ch, ref='RefSeq v1.0',
                 start=None, end=None, peak=None, marker='; '.join(f'{x}-{y}' for x, y, _ in sets), pve=f'up to {mx}', stat='', source='V_He2016', curation='hand')
        pts = [(m, x) for m, x in ((a, pa), (b, pb)) if x[0]]
        if arm == '5AL':  # main and point interval share Vrn-A1
            c = P('Ex_c7729_144', ch); pts.append(('Ex_c7729_144', c))
        if len(pts) >= 2:
            v = [x[0] for _, x in pts]; r['start'], r['end'] = min(v), max(v); r['grade'] = 'B' if r['end'] - r['start'] <= 100 else 'C'
            place = ' Placed in silico: ' + ', '.join(f'{m} {x[0]:.3f} Mb ({x[1]})' for m, x in pts) + '.'
        elif pts:
            m, x = pts[0]; r['peak'] = x[0]; r['grade'] = 'C'; other = b if m == a else a
            place = f' Only {m} could be placed on {ch} ({x[0]:.3f} Mb, {x[1]}); {other} maps elsewhere or has no position.'
        else:
            r['grade'] = 'D'; place = ' Neither flanking marker maps to this chromosome in Chinese Spring, so the locus is unplaced.'
        for x, y, l in sets[1:]:
            px, py = P(x, ch), P(y, ch)
            if arm != '5AL' or l.startswith('Table 5'):
                place += f' {l}: ' + ', '.join(f'{m} {q[0]:.3f} Mb' if q[0] else f'{m} not placed' for m, q in ((x, px), (y, py))) + '.'
        r['note'] = (f'QTL for {tr} in the Soru#1 x Naxos RIL population (131 lines; Mexico, Norway and China), resistance from {par}. '
                     f'R2 (%): {det}.' + ex + ' Flanks: ' + '; '.join(f'{x}-{y} ({l})' for x, y, l in sets) + '.' + place
                     + (f' Interval width {r["end"] - r["start"]:.0f} Mb: a chromosome region, not a locus.' if r.get('grade') == 'C' and r['start'] else '')
                     + ' No physical position printed. Corrected from the paper, Nov 2026 (ad19.py).')
        r['evidence'] = '; '.join(l for _, _, l in sets) + f'; R2 {det}'
        if arm == '2DLc': r['primers'] = KASP
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
