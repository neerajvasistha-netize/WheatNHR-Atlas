"""ad40 (Nov 2026): Buerstmayr and Buerstmayr 2015 TAG 128:1519, Capo x Arina, from Online Resource 3 (Fig. S1 linkage maps). Replaces auto records V1429-V1434 and adds 2A and 4A QTL. Anther retention and plant height (5AS) not entered."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {}
T = ' Read from Online Resource 3 (Fig. S1), Nov 2026 (ad40.py).'
def mk(i, **k):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand'); r.update(k); r['note'] += T; return r
CEN = {'1A': 213, '1B': 241, '1D': 166, '2A': 214, '2B': 343, '2D': 264, '3A': 316, '3B': 346, '3D': 237, '5A': 253, '5B': 200, '5D': 186, '6A': 283, '6B': 323, '6D': 211, '7A': 360, '7B': 300, '7D': 336, '4B': 303, '4D': 184}
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
IDS = ['V1429', 'V1430', 'V1431', 'V1432', 'V1433', 'V1434', 'S195-001', 'S195-002']
def build(reports):
    keep = [r for r in reports if r['id'] not in IDS]; N = []
    BU = 'FHB resistance (AUDPC, spray inoculation, 2011-2013) and anther retention in the Capo x Arina winter wheat population; MQM mapping on a DArT/SSR map. Flanking markers read from the LOD profiles in Online Resource 3 (Fig. S1); LOD and effect sizes from the main paper (not uploaded) are not entered. '
    Q = [('V1433', 'QFhb.ifa-1BS (Capo x Arina)', '1B', ['wPt-3103', 'wPt-6117', 'wPt-666564'], 'S', 'FHB QTL on 1BS near 20 cM (wPt-3103 to wPt-666564), LOD above 3 in 2011; coincides with the AUDPC QTL of Draeger et al. (2007).'),
         ('S195-001', 'QFhb.ifa-2A (Capo x Arina)', '2A', ['wPt-7721', 'wPt-8490', 'wPt-0568', 'wPt-8068'], None, 'FHB QTL at the distal end of 2A (0-10 cM, wPt-7721 to wPt-8068), LOD about 5.5 in 2012 and above 3 for the mean.'),
         ('V1432', 'QFhb.ifa-3B (Capo x Arina)', '3B', ['wPt-10323', 'wPt-2372', 'wPt-7688', 'wPt-4608'], None, 'FHB QTL on 3B peaking at wPt-2372 (about 93 cM), between wPt-10323 and wPt-7688, LOD about 4 in 2012; coincides with the AUDPC QTL of Paillard et al. (2004). wPt-4608 lies about 15 cM proximal of the peak.'),
         ('S195-002', 'QFhb.ifa-4A (Capo x Arina; with anther retention)', '4A', ['wPt-1082', 'Xbarc70', 'wPt-7327', 'wPt-2903', 'wPt-9901'], None, 'FHB QTL on 4A at about 63-72 cM (wPt-1082 to wPt-9901, near Xbarc70), LOD about 3.5 for FHB and about 7 for anther retention in the same interval; coincides with the AUDPC QTL of Paillard et al. (2004).'),
         ('V1431', 'QFhb.ifa-6BL (Capo x Arina; with anther retention)', '6B', ['rPt-1040', 'wPt-4648', 'tPt-8161', 'Xbarc24', 'tPt-3689'], None, 'FHB QTL on distal 6BL (about 80-91 cM, rPt-1040 to tPt-3689, near Xbarc24), LOD about 3.5 for the FHB mean and up to 6 for anther retention; coincides with the AUDPC, DON and relative spikelet weight QTL of Draeger et al. (2007).'),
         ('V1430', 'QFhb.ifa-3D (Capo x Arina)', '3D', ['wPt-741038', 'wPt-4544'], None, 'FHB QTL on the small 3D linkage group (wPt-741038 to wPt-4544), LOD about 4 in 2011 only.'),
         ('V1429', 'QFhb.ifa-5AL at B1 (Capo x Arina)', '5A', ['wPt-4203', 'Xgwm291', 'B1'], None, 'FHB QTL on distal 5AL at about 5-12 cM, at the awn inhibitor B1 and Xgwm291, LOD about 4.5 in 2012; coincides with the AUDPC QTL of Paillard et al. (2004). B1 itself is in the atlas at 698.53 Mb (v1).'),
         ('V1434', 'QFhb.ifa-7D (Capo x Arina)', '7D', ['wPt-743857', 'wPt-7508', 'wPt-664469'], None, 'FHB QTL at the distal end of the 7D linkage group (0-3 cM), LOD about 4.5 in 2011 only.')]
    for i, nm, ch, ms, arm, nt in Q:
        r = mk(i, name=nm, disease='Fusarium head blight', chr=ch, source='V_Buerstmayr2015', marker=', '.join(ms), note=BU + nt, evidence='Online Resource 3, Fig. S1')
        place(r, ms, ch, arm); N.append(r)
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
    json.dump(ops, open('/tmp/ops40.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], (r['asm'].get('RefSeq v2.1') or {}).get('flag'))
    print(len(ops), len(json.dumps(ops)))
