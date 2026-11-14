"""ad18 (Nov 2026): Lu et al. 2013 TAG 126:317 (FHB Type I, SHA3/CBRD x Naxos RIL; spawn, spray and point inoculation).
Replaces six auto records (V1421-V1426) with the 23 QTL of Tables 4 and 5. Positions are not printed; flanking markers are placed
in silico (exact-match PCR of SSR primers, Zhao et al. 2019 SSR table, DArT tag BLAST); Rht-B1 at its gene position.
Two flanks on the chromosome: grade B (C if > 100 Mb); one flank: grade C at that marker; none: grade D."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
OLD = ['V1421', 'V1422', 'V1423', 'V1424', 'V1425', 'V1426']
GENE = {'Rht-B1': (30.8614, 'Rht-B1 gene, TraesCS4B01G043100')}
# arm, flank1, flank2, R2 by trial (text), max R2, traits, source parent, associations
Q = [
 ('1AL.1', 'wPt-8797', 'wPt-7030', 'DON spray 2009 1.3, 2010 2.6, mean 2.0', 2.6, 'DON', 'Naxos', ''),
 ('1BS', 'gwm550', 'wmc619', 'severity spawn 2008 5.7; spray 2009 3.5, 2010 4.1, mean 4.7', 5.7, 'FHB severity', 'Naxos', ''),
 ('2AS', 'gwm636', 'barc124', 'FDK 2010 9.3; DON 2009 7.3, 2010 1.8, mean 4.1', 9.3, 'FDK and DON', 'Naxos', ''),
 ('2BL', 'wmc441', 'gwm1267b', 'DON spray 2009 5.2, mean 2.8', 5.2, 'DON', 'Naxos', ''),
 ('2DL', 'gwm265', 'mag3616', 'severity spawn 2008 3.1; spray 2009 2.0, 2010 2.9, mean 2.0', 3.1, 'FHB severity', 'Naxos', ''),
 ('3AS', 'wmc489b', 'wmc695b', 'DON spray 2010 4.4, mean 2.2', 4.4, 'DON', 'Naxos', ''),
 ('4BS', 'Rht-B1', 'gwm368', 'severity spawn 2008 11.2, 2011 1.6, mean 3.4; spray 2009 4.3, 2010 10.8, mean 6.5', 11.2, 'FHB severity', 'Naxos',
  ' At the Rht-B1 locus; associated with anther extrusion and plant height (Rht-B1b lines more susceptible after spawn and spray inoculation).'),
 ('5AS', 'wmc489d', 'wPt-8226', 'FDK 2010 11.5; DON 2009 7.1, 2010 14, mean 10.8', 14.0, 'FDK and DON', 'Naxos', ' The most important FDK/DON QTL, close to Xgwm156.'),
 ('5BL', 'barc275', 'barc232', 'severity spawn 2008 7.0, mean 6.1; DON 2009 2.9, mean 1.8', 7.0, 'FHB severity and DON', 'Naxos', ' Associated with anther extrusion.'),
 ('5DL', 'gwm174', 'wPt-1400', 'severity spray 2010 3.3, mean 3.5', 3.5, 'FHB severity', 'Naxos', ''),
 ('7AL.1', 'wmc603', 'barc292', 'DON spray 2009 16.2, 2010 7.8, mean 11.5', 16.2, 'DON', 'Naxos', ' Major DON QTL near the centromere close to Xwmc603, stable over two years.'),
 ('1AL.2', 'wPt-8016', 'wPt-2847', 'severity spawn 2011 8.3, mean 9.3', 9.3, 'FHB severity', 'SHA3/CBRD', ''),
 ('2DLc', 'wmc18', 'wmc41', 'severity spawn 2008 7.5, 2011 22.3, mean 24.3; spray 2009 2.0, 2010 12.4, mean 5.0; DON 2009 7.4; point 2010 4.2', 24.3, 'FHB severity and DON', 'SHA3/CBRD',
  ' The major and consistent QTL, near the centromere close to Xgwm539; effective with all three inoculation methods; associated with anther extrusion and plant height.'),
 ('3DL', 'cfd9', 'barc323', 'severity spawn 2008 3.8, mean 1.7; spray 2009 2.6, 2010 0.9, mean 2.0', 3.8, 'FHB severity', 'SHA3/CBRD', ' Associated with anther extrusion.'),
 ('4AL', 'gwm160', 'wPt-5172', 'severity spray 2010 10.5', 10.5, 'FHB severity', 'SHA3/CBRD', ' Associated with plant height.'),
 ('5AL', 'gwm617', 'gwm291', 'severity spawn 2011 5.7, mean 6.8; spray 2009 2.4, 2010 6.5, mean 3.6', 6.8, 'FHB severity', 'SHA3/CBRD', ' Associated with plant height.'),
 ('6AS', 'wPt-0832', 'wPt-6904', 'severity spawn 2008 7.1; spray 2009 5.3, mean 4.2', 7.1, 'FHB severity', 'SHA3/CBRD', ' Associated with plant height.'),
 ('6ASc', 'barc37', 'wmc748a', 'DON spray 2009 2.6, mean 2.1', 2.6, 'DON', 'SHA3/CBRD', ' Associated with plant height.'),
 ('7AL.2', 'barc121', 'wPt-8399', 'severity spawn 2011 1.8, mean 2.8', 2.8, 'FHB severity', 'SHA3/CBRD', ' Associated with anther extrusion.'),
 ('1DS', 'wmc432', 'barc152', 'point severity 2009 3.7', 3.7, 'FHB severity (point)', 'Naxos', ''),
 ('2AL', 'gwm328', 'gdm93', 'point severity 2010 3.9, mean 4.4', 4.4, 'FHB severity (point)', 'Naxos', ''),
 ('2BL (point)', 'wmc441', 'mag548a', 'point severity 2009 2.8, 2010 8.9, mean 9.4', 9.4, 'FHB severity (point)', 'Naxos', ' Same region as the spray-inoculation DON QTL on 2BL.'),
 ('2DS', 'gwm296', 'wPt-11625', 'point severity 2009 8.9, 2010 4.8, mean 10.3; DON 2009 4.1', 10.3, 'FHB severity and DON (point)', 'Naxos', ''),
 ('4DL', 'barc98', 'cfd71', 'point severity 2009 4.6, 2010 5.3, mean 5.9; DON 2009 3.8', 5.9, 'FHB severity and DON (point)', 'SHA3/CBRD', ''),
]
def R(m, ch):
    if m in GENE: return GENE[m]
    p, h = resolve(m, ch)
    if not p and re.search(r'\d[a-z]$', m):
        p, h = resolve(m[:-1], ch)
        if p: h += '; locus-copy suffix dropped'
    return p, h
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S95-')]; new = []
    for n, (arm, a, b, det, mx, tr, par, ex) in enumerate(Q, 1):
        ch = arm[:2]; pa, ha = R(a, ch); pb, hb = R(b, ch)
        r = dict(id=f'S95-{n:03d}', name=f'QFhb {arm} (SHA3/CBRD x Naxos)', cls='QTL', disease='Fusarium head blight', chr=ch, ref='RefSeq v1.0',
                 start=None, end=None, peak=None, marker=f'{a}-{b}', pve=f'up to {mx}', stat='', source='V_Lu2012', curation='hand')
        pl = [(m, p, h) for m, p, h in ((a, pa, ha), (b, pb, hb)) if p]
        if len(pl) == 2:
            r['start'], r['end'] = min(pa, pb), max(pa, pb); r['grade'] = 'B' if r['end'] - r['start'] <= 100 else 'C'
            place = f' Flanks placed in silico: {a} {pa:.3f} Mb ({ha}), {b} {pb:.3f} Mb ({hb}).'
            if r['grade'] == 'C': place += f' Interval width {r["end"] - r["start"]:.0f} Mb: a chromosome region, not a locus.'
        elif pl:
            m, p, h = pl[0]; r['peak'] = p; r['grade'] = 'C'
            place = f' Only {m} could be placed ({p:.3f} Mb, {h}); position is that flanking marker.'
        else:
            r['grade'] = 'D'; place = ' Neither flanking marker could be placed in silico, so the locus is unplaced.'
        r['note'] = (f'QTL for {tr} in the SHA3/CBRD x Naxos RIL population (Norway; spawn, spray and point inoculation), resistance from {par}. '
                     f'R2 (%): {det}.' + ex + place + ' No physical position printed. Corrected from the paper, Nov 2026 (ad18.py).')
        r['evidence'] = f'Table {"5" if "point" in tr else "4"}: {arm} {a}-{b}; {det}'
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
    json.dump(ops, open('/home/claude/work/ad18.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['start'], r['end'], r['peak'], r['grade'], (r['asm'].get('RefSeq v2.1') or {}).get('flag'))
