"""ad33 (Nov 2026): eight uploads. Added: Li et al. 2026 TAG 139:303 (SbWE35 and Sb3), Ma et al. 2010 TAG 120:1119 (Qcrs.cpi-3B,
Qcrs.cpi-4B), Ma et al. 2012 Euphytica 183:147 (Qcrs.cpi-3B NILs).
Not added: Lorrain et al. preprint and Nature Plants 2026 (pathogen genes), Mackintosh et al. 2006 Plant Cell Rep 25:313 (transgenic
overexpression), Ma et al. 2015 MGG 290:1955 (MDA content QTL under drought; links to Qcrs.cpi-3B but is not a disease locus),
Loffler et al. 2009 Mol Breed 23:473 (meta-QTL by other workers, excluded by atlas rule)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Li2026TAG': ('Li B, Xie J, Guo G, Zhang P, Zhang H, Zhu K, Huang B, Wu Q, Chen Y, Yuan C, Dong L, Liu Z, Lu P (2026) Two independent dominant resistance genes confer resistance to spot blotch in wheat. Theor Appl Genet 139:303', 'https://doi.org/10.1007/s00122-026-05412-5', 'Li et al. (2026)'),
 'Ma2010TAG': ('Ma J, Li HB, Zhang CY, Yang XM, Liu YX, Yan GJ, Liu CJ (2010) Identification and validation of a major QTL conferring crown rot resistance in hexaploid wheat. Theor Appl Genet 120:1119-1128', 'https://doi.org/10.1007/s00122-009-1239-3', 'Ma et al. (2010)'),
 'Ma2012Euph': ('Ma J, Yan GJ, Liu CJ (2012) Development of near-isogenic lines for a major QTL on 3BL conferring Fusarium crown rot resistance in hexaploid wheat. Euphytica 183:147-152', 'https://doi.org/10.1007/s10681-011-0414-1', 'Ma et al. (2012)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad33.py).'
def mk(i, **k):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand'); r.update(k); r['note'] += T; return r
def build(reports):
    keep = [r for r in reports if not re.match(r'S14[6-8]-', r['id'])]; N = []
    N.append(mk('S146-001', name='SbWE35', cls='DG', disease='Spot blotch', chr='4B', source='Li2026TAG', start=581.159216, end=582.113681, grade='A', marker='YK12918-YK13104 (co-segregating YK4BM4, YK12920, YK12922)',
        note='Dominant spot blotch resistance gene SbWE35 from WE35 (BSR-Seq on WE35 x 621-18-3; 4,500 F2 plants for fine mapping) delimited to 0.07 cM between YK12918 and YK13104, a 1.08 Mb segment of Chinese Spring 4BL with 12 high-confidence genes (TraesCS4B01G295600 to TraesCS4B01G296700); five are expressed: a cation/calcium exchanger, a DUF1677 protein, a NEDD4-binding protein 2, an alpha/beta hydrolase and a glutathione-related gene. The interval overlaps Sb4, so SbWE35 and Sb4 may be the same locus (allelism not yet tested). Placed across the RefSeq v1.0 positions of the interval genes (581.16-582.11 Mb).',
        evidence='Results; Fig. 4; Table 3'))
    N.append(mk('S146-002', name='Sb3 (WE35)', cls='DG', disease='Spot blotch', chr='3B', source='Li2026TAG', start=6.0887, end=6.4501, grade='C', marker='Sb3 markers (Lu et al. 2016)',
        note='In WE35 x S2199 F2 and F2:3 families resistance segregated as two independent dominant genes, the known Sb3 (3BS) and SbWE35 (4BL); WE35 carries Sb3. Placed at the Sb3 interval in the atlas (Lu et al. 2016).',
        evidence='Results: two-gene segregation, Table 1'))
    MA = 'Crown rot resistance QTL in Lang x CSCR6 RILs (seedling assays with F. pseudograminearum CS3096 and F. graminearum CS3005; DArT and SSR map). '
    N.append(mk('S147-001', name='Qcrs.cpi-3B (CSCR6; original mapping)', disease='Fusarium crown rot', chr='3B', source='Ma2010TAG', start=820.501, end=822.202, grade='C', marker='wPt11143/wPt10505-wPt2277; gwm0181 (2.2 cM)', pve='29.7-48.8', stat='LOD up to 12.6',
        note=MA + 'Major QTL on 3BL from CSCR6 (T. spelta) detected in all six trials, explaining up to 48.8% of the variance (IM), with a strong interaction with Qcrs.cpi-4B; its effect was validated with gwm0181 in four other populations. The DArT flanks cannot be placed; placed at the fine-mapped Qcrs.cpi-3B interval in the atlas (Zheng et al. 2015).',
        evidence='Table 3: Qcrs.cpi-3B, wPt10505-wPt2277, LOD 10.8-12.6, R2 up to 48.8%'))
    pa, ha = resolve('Xbarc199', '4B'); pb, hb = resolve('Xgwm495', '4B')
    N.append(mk('S147-002', name='Qcrs.cpi-4B (Lang)', disease='Fusarium crown rot', chr='4B', source='Ma2010TAG', start=min(pa, pb), end=max(pa, pb), grade='C', marker='barc0199-gwm0495 (wPt5334-wPt4918)', pve='13.5-22.8', stat='LOD 3.7-4.7',
        note=MA + f'Minor QTL from Lang on 4B (0-11.8 cM), detected with CS3005 in two trials, up to 22.8% (IM). Placed between barc0199 ({pa:.1f} Mb, {ha}) and gwm0495 ({pb:.1f} Mb, {hb}): a chromosome region, not a locus.',
        evidence='Table 3: Qcrs.cpi-4B, barc0199-gwm0495, R2 18.1% (combined)'))
    N.append(mk('S148-001', name='Qcrs.cpi-3B NILs (3BL)', disease='Fusarium crown rot', chr='3B', source='Ma2012Euph', start=820.501, end=822.202, grade='C', marker='gwm0181',
        note='Near-isogenic lines for the 3BL crown rot QTL (Qcrs.cpi-3B) developed by heterogeneous inbred family analysis from three populations using the linked SSR gwm0181: resistance differed significantly between isolines in nine of 13 putative NIL pairs, confirming the locus. Placed at the fine-mapped Qcrs.cpi-3B interval in the atlas.',
        evidence='Results: nine of 13 NIL pairs differed'))
    for r in N: r['asm'] = lift.lift_record(r)
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
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1',) if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/tmp/ops33.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    print(open('/tmp/ops33.json').read())
