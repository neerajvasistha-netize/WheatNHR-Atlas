"""ad45 (Nov 2026): Martin et al. 2015 TAG 128:377 crown rot QTL refined with Supplementary Figure 1 (ESM1 linkage maps of the four populations;
marker names and QTL bars read from the embedded EMF maps). Markers inside each QTL bar, or the nearest flanking markers, placed in silico.
Adds the two Janz-derived QTL (2BS in 2-49/Janz, 2DS in CPI133814/Janz) that the main paper reports."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
T = ' Read from the paper and Supplementary Figure 1, Nov 2026 (ad45.py).'
RHTB1 = 30.861268
def pl(ms, ch):
    out = []
    for m in ms:
        p, h = resolve(m, ch)
        if p: out.append((m, p, h))
    return out
def txt(ok): return ', '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in ok)
U = {  # id: (markers used, which, map text)
 'S150-001': (['Xwmc24', 'Xbarc148'], 'range', 'In Fig. S1 the 1A QTL bars lie at about 35-41 cM (2-49/Janz; barc148-gwm164-wmc120) and 46-52 cM (Sunco/2-49); wmc24 is the distal flank (24 cM). gwm164 (280.6 Mb) and wmc120 (344.7 Mb) place on 1AL and were not used.'),
 'S150-002': (['Xgwm413', 'Xgwm11'], 'rangeC', 'In Fig. S1 the 1B QTL bar lies at about 15-17 cM (2-49/Janz), at cfd65/gwm11, between gwm413 (0 cM) and gwm131 (23 cM); cfd59, cfd65 and gwm131 could not be placed, so the QTL is placed between gwm413 and gwm11.'),
 'S150-003': (['Xwmc429', 'Xwmc216'], 'range', 'In Fig. S1 the 1D seedling QTL bar lies at cfd19-wmc216 (about 46-49 cM in 2-49/Janz, about 67-74 cM in Sunco/2-49), with wmc429 just distal.'),
 'S150-004': (['wPt-0615'], 'peak', 'In Fig. S1 the 2B seedling and field QTL bars (Sunco/2-49) lie at about 13-15 cM, at barc55, wPt-2907 and wPt-5287; none of these could be placed, so the QTL is placed at the nearest placeable marker, wPt-0615 (about 9 cM).'),
 'S150-006': (['Xgwm285', 'Xgwm376'], 'rangeC', 'In Fig. S1 the 3B seedling and field QTL bars (Sunco/2-49) lie at about 68-77 cM, from gwm131/gwm376 to gwm285, near the centromere: gwm376 places at 389.5 Mb and gwm285 at 273.1 Mb.'),
 'S150-007': (['Rht-B1'], 'rht', 'In Fig. S1 the 4B field QTL bar (2-49/Janz) lies at about 18-22 cM (wmc467, gwm513, wPt-7569PCR, gwm165) on a map anchored by Rht-B1b at 0 cM; gwm513, gwm165 and wmc349 place on 4BL in silico, contradicting the map, so the QTL is placed at Rht-B1 (TraesCS4B02G043100, 30.861 Mb), the map anchor about 20 cM distal.'),
 'S150-008': (['Xgwm340', 'Xgwm181', 'wPt-0365'], 'range', 'In Fig. S1 the 3BL seedling QTL bar (IRN497/Janz) lies at about 186-198 cM, ending at wPt-0365PCR and just distal of gwm340/gwm181 at the end of 3BL.'),
 'S150-013': (['wPt-5947', 'Xgwm299'], 'range', 'In Fig. S1 the 3BL seedling QTL bar (CPI133814/Janz) lies at about 132-145 cM, from wPt-0021PCR to gwm299 (wPt-3822, wPt-5947, wPt-8559, wPt-4194 inside).'),
 'S150-014': (['Xcfd76', 'Xbarc273', 'Xbarc175'], 'range', 'In Fig. S1 the 6D field QTL bar (CPI133814/Janz) lies at about 98-108 cM, over wmc416, barc196, barc175, barc273 and cfd76.'),
}
NEW = [('S150-015', 'QCr 2BS (Janz; 2-49/Janz; Martin 2015)', '2B', ['Xwmc272'], 'Janz', 'peak',
        'Seedling QTL from Janz in 2-49/Janz, gdm86-cfa2278, 8.8% (main paper).|In Fig. S1 the bar lies at about 89-95 cM (gdm86, gwm271, cfa2278), none placeable, so the QTL is placed at the distal flank wmc272 (85 cM).', '8.8'),
       ('S150-016', 'QCr 2DS (Janz; CPI133814/Janz; Martin 2015)', '2D', ['Xgwm484', 'Xgwm102'], 'Janz', 'range',
        'Field QTL from Janz in CPI133814/Janz, 2012, 12.1% (main paper).|In Fig. S1 it lies at about 33-40 cM between gwm484 and gwm102, the same region as the Sunco 2DS QTL.', '12.1')]
MR = 'Crown rot (F. pseudograminearum) QTL in four DH populations from partially resistant sources 2-49, Sunco, IRN497 and CPI133814 (crossed with Janz or each other), seedling and field trials in Queensland. '
def apply(r, ms, mode, maptxt):
    base = r['note'].split(' Only ')[0].split(' Placed in silico')[0].split(' None of the flanking')[0]
    for k in ('start', 'end', 'peak'): r.pop(k, None)
    if mode == 'rht':
        r['peak'], r['grade'] = RHTB1, 'C'; t = ''
    else:
        ok = pl(ms, r['chr'])
        if mode == 'peak' or len(ok) == 1:
            m, p, h = ok[0]; r['peak'], r['grade'] = p, 'C'; t = f' Placed at {m} ({p:.3f} Mb, {h}).'
        else:
            lo, hi = min(x[1] for x in ok), max(x[1] for x in ok); r['start'], r['end'] = lo, hi
            r['grade'] = 'B' if (hi - lo <= 100 and mode == 'range') else 'C'; t = ' Placed in silico: ' + txt(ok) + '.' + (f' The markers span {hi - lo:.0f} Mb: a chromosome region, not a locus.' if hi - lo > 100 else '')
    r['note'] = base + ' ' + maptxt + t + T
    r['marker'] = r['marker'] + ' (Fig. S1: ' + ', '.join(ms) + ')' if 'Fig. S1' not in r['marker'] else r['marker']
    r['evidence'] = r['evidence'].split('; Fig. S1')[0] + '; Fig. S1 (ESM1)'
    r['asm'] = lift.lift_record(r)
def build(reports):
    out, upd = [], []
    for r in reports:
        if r['id'] in ('S150-015', 'S150-016'): continue
        if r['id'] in U:
            r = dict(r); apply(r, *U[r['id']]); upd.append(r)
        out.append(r)
    for i, nm, ch, ms, par, mode, nt, pv in NEW:
        r = dict(id=i, cls='QTL', ref='RefSeq v1.0', curation='hand', name=nm, disease='Fusarium crown rot', chr=ch, source='Martin2015TAG', marker=', '.join(ms), pve=pv,
                 note=MR + f'Resistance from {par}. ' + nt.split('|')[0], evidence='Results; Fig. S1')
        apply(r, ms, mode, nt.split('|')[1]); r['marker'] = ', '.join(ms); out.append(r); upd.append(r)
    return out, upd
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], upd = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], _ = build(Sd['reports'])
    open('web/index.html', 'w').write(H[:m.start(2)] + json.dumps(Sd, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/') + H[m.end(2):])
    ops = []
    for r in upd:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', '_note') if k in r['asm']}
        o = {'op': 'set', 'path': 'reports/' + r['id'], 'data': d}
        if r['id'] in U: o['ifVersion'] = 1
        ops.append(o); print(r['id'], r['grade'], r.get('start'), r.get('end'), r.get('peak'), d['asm'].get('RefSeq v2.1', {}).get('flag'))
    json.dump(ops, open('/tmp/ops45.json', 'w'), ensure_ascii=False)
