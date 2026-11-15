"""ad35 (Nov 2026): Puskas et al. 2006 CRC (FHB, Ning8331) and Kalih et al. 2014 TAG (triticale FHB, Ddw1)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Puskas2006CRC': ('Puskas K, Vida G, Komaromi J, Buerstmayr H, Lemmens M, Bedo Z, Veisz O (2006) Study of Fusarium head blight resistance in wheat using microsatellite markers. Cereal Res Commun 34(1):629', 'https://doi.org/10.1556/CRC.34.2006.1.157', 'Puskas et al. (2006)'),
 'Kalih2014TAG': ('Kalih R, Maurer HP, Hackauf B, Miedaner T (2014) Effect of a rye dwarfing gene on plant height, heading stage, and Fusarium head blight in triticale (xTriticosecale Wittmack). Theor Appl Genet 127:1527-1536', 'https://doi.org/10.1007/s00122-014-2316-9', 'Kalih et al. (2014)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad35.py).'
def mk(i, **k):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand'); r.update(k); r['note'] += T; return r
CEN = {'1A': 213, '1B': 241, '1D': 166, '2A': 214, '2B': 343, '2D': 264, '3A': 316, '3B': 346, '3D': 237, '5A': 253, '5B': 200, '5D': 186, '6A': 283, '6B': 323, '6D': 211, '7A': 360, '7B': 300, '7D': 336, '4B': 303, '4D': 184}
def place(r, ms, ch, arm=None):
    pl = [(m,) + resolve(m, ch) for m in ms]; pl = [x for x in pl if x[1]]
    drop = []
    if arm and ch in CEN:
        drop = [x for x in pl if (x[1] > CEN[ch]) != (arm == 'L')]; pl = [x for x in pl if x not in drop]
    if len(pl) >= 2:
        lo, hi = min(x[1] for x in pl), max(x[1] for x in pl); r['start'], r['end'] = lo, hi; r['grade'] = 'B' if hi - lo <= 100 else 'C'
        t = ' Placed in silico: ' + ', '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in pl) + '.' + (f' The markers span {hi - lo:.0f} Mb: a chromosome region, not a locus.' if r['grade'] == 'C' else '')
    elif pl: m, p, h = pl[0]; r['peak'], r['grade'] = p, 'C'; t = f' Only {m} could be placed ({p:.3f} Mb, {h}); position is that marker.'
    else: r['grade'] = 'D'; t = ' None of the flanking markers could be placed' + (' on the stated arm.' if drop else '.')
    if drop: t += ' Placements off the stated chromosome arm were discarded (approximate centromere ' + str(CEN[ch]) + ' Mb): ' + ', '.join(f'{m} {p:.1f} Mb' for m, p, h in drop) + '.'
    r['note'] = r['note'][:-len(T)] + t + T
def build(reports):
    keep = [r for r in reports if not re.match(r'S15[56]-', r['id'])]; N = []
    PU = 'Field FHB resistance (F. culmorum spray, mist irrigation, Martonvasar) in 228 F5-derived RILs of Ning8331/Martonvasari 17, single-marker regression with 63 SSRs (two-year means); all effects from Ning8331. '
    r = mk('S155-001', name='3BS FHB QTL (Ning8331)', disease='Fusarium head blight', chr='3B', source='Puskas2006CRC', marker='Xgwm493, Xgwm389', pve='up to 17.7', stat='P 2.7e-7 (Gwm493, kernel infection)',
           note=PU + 'Loci on 3BS were most closely related to FHB resistance: Gwm493 with spike infection (R2 6.5%) and kernel infection (14.7%); Gwm389 with loss in kernel weight per head (14.7%), test weight (15.1%) and thousand kernel weight (17.7%).', evidence='Table 1; Results')
    place(r, ['Xgwm493', 'Xgwm389'], '3B'); N.append(r)
    r = mk('S155-002', name='2BS FHB QTL (Ning8331)', disease='Fusarium head blight', chr='2B', source='Puskas2006CRC', marker='Xgwm429, Xgwm257', pve='up to 14.0', stat='P 1.9e-4 (Gwm429, kernel infection)',
           note=PU + 'A substantial effect on 2BS: Gwm429 with kernel infection (R2 8.3%) and loss in test weight (9.7%); Gwm257 with loss in kernel weight per head (14.0%) and thousand kernel weight (8.9%).', evidence='Table 1; Results')
    place(r, ['Xgwm429', 'Xgwm257'], '2B'); N.append(r)
    r = mk('S155-003', name='6A FHB field resistance (Ning8331, Gwm169)', disease='Fusarium head blight', chr='6A', source='Puskas2006CRC', marker='Xgwm169',
           note=PU + 'Confirmed an earlier report of a field resistance effect linked to Gwm169 on 6A, far weaker than the 3B and 2B effects. Minor effects on 5A, 5B, 6D, 4B, 7A, 2D, 1D, 3A and 6B were reported without markers and are not entered.', evidence='Results')
    place(r, ['Xgwm169'], '6A'); N.append(r)
    KA = 'FHB severity (F. culmorum spray) in 199 triticale DH lines of HeTi117-06 x Pigmej at two locations in two years, composite interval mapping on a DArT/SSR map. '
    N.append(mk('S156-001', name='C_5A FHB QTL (triticale HeTi117-06 x Pigmej)', disease='Fusarium head blight', chr='5A', source='Kalih2014TAG', grade='D', marker='wPt-9204, wPt-6825', pve='10.17', stat='LOD 3.34',
        note=KA + 'QTL on 5A at 7 cM (CI 2-11), additive effect 1.141, 10.2% of genotypic variance, found in 21% of cross-validation runs. The flanking DArT markers could not be placed.', evidence='Table 3'))
    N.append(mk('S156-002', name='C_5R FHB QTL at Ddw1 (triticale, rye 5RL)', disease='Fusarium head blight', chr='5R', source='Kalih2014TAG', grade='D', ref='alien (rye 5RL)', marker='Xiac130, Xiac128', pve='48.40', stat='LOD 14.29',
        note=KA + 'Major QTL on rye chromosome 5RL at 112 cM (CI 111-114), 48.4% of genotypic variance, 95% cross-validation frequency; coincides with the rye dwarfing gene Ddw1, whose dwarfing allele (from Pigmej) increased FHB severity, reduced height and delayed heading. Rye chromosome: no wheat RefSeq position.', evidence='Table 3; Table 4; Fig. 2'))
    N.append(mk('S156-003', name='C_6R FHB QTL (triticale, rye 6R)', disease='Fusarium head blight', chr='6R', source='Kalih2014TAG', grade='D', ref='alien (rye 6R)', marker='rPt-401114, tPt-513005', pve='9.58', stat='LOD 4.06',
        note=KA + 'QTL on rye chromosome 6R at 89 cM (CI 88-90), additive effect -1.113, 9.6% of genotypic variance, with a 5R x 6R epistatic interaction (6.1%). Rye chromosome: no wheat RefSeq position.', evidence='Table 3'))
    for r in N: r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}
    return keep + N, N
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], new = build(M['reports']); M['sources'].update(SRC); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], new2 = build(Sd['reports'])
    if 'sources' in Sd: Sd['sources'].update(SRC)
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    ops = [{'op': 'set', 'path': 'sources/' + k, 'data': v} for k, v in SRC.items()]
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/tmp/ops35.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], (r['asm'].get('RefSeq v2.1') or {}).get('flag'))
    print(len(ops), len(json.dumps(ops)))
