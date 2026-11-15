"""ad34 (Nov 2026): eight uploads. Added: McCartney et al. 2003 TAG 107:1181 (Stb7), Martin et al. 2015 TAG 128:377 (crown rot QTL in
four sources), McCartney et al. 2007 Mol Breed 20:209 (FHB QTL alleles in Canadian backgrounds), Miedaner et al. 2006 TAG 112:562 (QTL
stacking), Meyer et al. 2011 Euphytica 177:267 (Pch1 markers), Matus-Cadiz et al. 2006 Can J Plant Sci 86:711 (Fhb1 multiplex PCR).
Not added: Mahboubi et al. preprint (already in atlas from the published J Appl Genet 2022 version), McCartney et al. 2004 TAG 109:261
(haplotype survey of known QTL markers; no new loci)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'McCartney2003TAG': ('McCartney CA, Brule-Babel AL, Lamari L, Somers DJ (2003) Chromosomal location of a race-specific resistance gene to Mycosphaerella graminicola in the spring wheat ST6. Theor Appl Genet 107:1181-1186', 'https://doi.org/10.1007/s00122-003-1359-0', 'McCartney et al. (2003)'),
 'Martin2015TAG': ('Martin A, Bovill WD, Percy CD, Herde D, Fletcher S, Kelly A, Neate SM, Sutherland MW (2015) Markers for seedling and adult plant crown rot resistance in four partially resistant bread wheat sources. Theor Appl Genet 128:377-385', 'https://doi.org/10.1007/s00122-014-2437-1', 'Martin et al. (2015)'),
 'McCartney2007MB': ('McCartney CA, Somers DJ, Fedak G, DePauw RM, Thomas J, Fox SL, Humphreys DG, Lukow O, Savard ME, McCallum BD, Gilbert J, Cao W (2007) The evaluation of FHB resistance QTLs introgressed into elite Canadian spring wheat germplasm. Mol Breed 20:209-221', 'https://doi.org/10.1007/s11032-007-9084-z', 'McCartney et al. (2007)'),
 'Miedaner2006TAG': ('Miedaner T, Wilde F, Steiner B, Buerstmayr H, Korzun V, Ebmeyer E (2006) Stacking quantitative trait loci (QTL) for Fusarium head blight resistance from non-adapted sources in an European elite spring wheat background and assessing their effects on deoxynivalenol (DON) content and disease severity. Theor Appl Genet 112:562-569', 'https://doi.org/10.1007/s00122-005-0163-4', 'Miedaner et al. (2006)'),
 'Meyer2011Euph': ('Meyer N, Lind V, Heindorf M, Korzun V, Friedt W, Ordon F (2011) Diagnostic value of molecular markers linked to the eyespot resistance gene Pch1 in wheat. Euphytica 177:267-275', 'https://doi.org/10.1007/s10681-010-0266-0', 'Meyer et al. (2011)'),
 'MatusCadiz2006': ('Matus-Cadiz MA, Pozniak CJ, Hughes GR, Hucl P (2006) Multiplexed PCR to screen for a major QTL carrying fusarium head blight resistance in Sumai-3 wheat. Can J Plant Sci 86:711-716', 'https://doi.org/10.4141/P05-153', 'Matus-Cadiz et al. (2006)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad34.py).'
FHB1, FHB5 = 8.528, (105.4338, 214.2)
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
MQ = [('1AS', '1A', ['Xbarc148', 'Xgwm164'], '2-49', '2-49/Janz seedling 2001 LOD 6.0 (12.8%), field 2004 3.2 (16.4%), 2006 5.1 (11.3%); Sunco/2-49 field 2009 3.4 (8.5%)', '16.4'),
      ('1BS', '1B', ['Xgwm11', 'Xcfd65'], '2-49', 'Sunco/2-49 seedling 2007 LOD 3.3 (7.2%); 2-49/Janz field 2005 5.0 (12.6%)', '12.6'),
      ('1DL', '1D', ['Xcfd19', 'Xwmc216', 'Xbarc162'], '2-49', '2-49/Janz seedling 2001 LOD 8.0 (17.4%), 2006 4.8 (13.6%); Sunco/2-49 seedling 2009 3.3 (12.6%); seedling only', '17.4'),
      ('2BS', '2B', ['Xgwm630', 'Xcfa2278'], 'Sunco', 'Sunco/2-49 seedling 2007 LOD 3.2 (7.3%), field 2009 4.7 (12.2%); on a T. timopheevii introgression; a QTL in the same region came from Janz in 2-49/Janz (gdm86-cfa2278, 8.8%)', '12.2'),
      ('2DS', '2D', ['Xgwm484', 'Xgwm102'], 'Sunco', 'Sunco/2-49 field 2009 LOD 3.4 (8.6%); field only; a QTL in the same region came from Janz in CPI133814/Janz (field 2012, 12.1%)', '12.1'),
      ('3BS', '3B', ['Xgwm131'], '2-49', 'Sunco/2-49 field 2009 LOD 4.2 (11.1%), flanked by gwm131 and wPt-9310PCR', '11.1'),
      ('4BS', '4B', ['Xwmc467', 'Xgwm165', 'Xwmc657'], '2-49', '2-49/Janz field 2004 LOD 3.1 (13.8%), 2005 4.4 (11.0%), 2006 8.2 (18.3%); Sunco/2-49 seedling 2008 7.1 (20.4%); flanks wPt-7569PCR with wmc467, gwm165 or wmc657', '20.4'),
      ('3BL (IRN497)', '3B', ['Xwmc236'], 'IRN497', 'IRN497/Janz seedling 2003 LOD 5.8 (18.8%), flanked by wmc236 and wPt-0365PCR; seedling only', '18.8'),
      ('6DL (IRN497)', '6D', ['Xcfd188', 'Xcfd47'], 'IRN497', 'IRN497/Janz seedling 2003 LOD 2.6 (7.8%), field 2006 2.3 (7.1%); suggestive', '7.8'),
      ('2AL (IRN497)', '2A', ['Xgwm95', 'Xcfa2043'], 'IRN497', 'IRN497/Janz field 2006 LOD 5.6 (19.7%); adult plant only', '19.7'),
      ('4BS (IRN497)', '4B', ['Xbarc193', 'Xwmc349'], 'IRN497', 'IRN497/Janz field 2007 LOD 2.8 (8.6%); flanks 25 cM apart', '8.6'),
      ('3AL (CPI133814)', '3A', ['Xcfa2134', 'Xcfa2262'], 'CPI133814', 'CPI133814/Janz seedling 2008 LOD 3.2 (11.9%)', '11.9'),
      ('3BL (CPI133814)', '3B', ['Xgwm299'], 'CPI133814', 'CPI133814/Janz seedling 2008 LOD 3.4 (12.1%), flanked by wPt-0021PCR and gwm299; same region as the IRN497 3BL QTL', '12.1'),
      ('6DL (CPI133814)', '6D', ['Xbarc196', 'Xbarc273'], 'CPI133814', 'CPI133814/Janz field 2011 LOD 4.6 (18.6%)', '18.6')]
def build(reports):
    keep = [r for r in reports if not re.match(r'S1(49|5[0-4])-', r['id'])]; N = []
    p, h = resolve('Xwmc313', '4A')
    N.append(mk('S149-001', name='Stb7 (ST6)', cls='DG', disease='Septoria tritici blotch', chr='4A', source='McCartney2003TAG', peak=p, grade='C', marker='Xwmc313 (0.3-0.5 cM), Xwmc219, Xgwm160',
        note='Race-specific resistance of spring wheat ST6 to M. graminicola isolate MG2, a single gene mapped on 4AL in ST6/Erik and ST6/Katepwa F2:3 families and designated Stb7; Xwmc313, Xwmc219 and Xgwm160 mapped within 3.5 cM, none flanking; Xwmc313 was closest (0.3 and 0.5 cM) and its ST6 allele was absent from other tested lines. '
             f'Placed at Xwmc313 ({p:.3f} Mb, {h}).', evidence='Fig. 2: 4AL linkage map, Xwmc313 0.3-0.5 cM from Stb7'))
    MR = 'Crown rot (F. pseudograminearum) QTL in four DH populations from partially resistant sources 2-49, Sunco, IRN497 and CPI133814 (crossed with Janz or each other), seedling and field trials in Queensland. '
    for n, (nm, ch, ms, par, det, mx) in enumerate(MQ, 1):
        r = mk(f'S150-{n:03d}', name=f'QCr {nm.split(" (")[0]} ({par}; Martin 2015)', disease='Fusarium crown rot', chr=ch, source='Martin2015TAG', marker=', '.join(ms), pve=f'up to {mx}',
               note=MR + f'Resistance from {par}. Detections: {det}.', evidence=f'Table 5: {nm}')
        a = re.search(r'[1-7][ABD]([SL])', nm); place(r, ms, ch, a.group(1) if a else None); N.append(r)
    MC = 'FHB resistance QTL alleles of Nyuubai, Sumai-3 and Wuhan-1 introgressed into three elite Canadian spring wheat backgrounds (BC2 populations, FHB nurseries and DON, FDK, agronomic and quality traits). '
    r = mk('S151-001', name='Wuhan-1 4B FHB allele (Canadian backgrounds)', disease='Fusarium head blight', chr='4B', source='McCartney2007MB', marker='Xwmc238, Xgwm149',
           note=MC + 'The Wuhan-1 4B allele was the most effective FHB resistance allele but increased plant height by 9.3 cm.'); place(r, ['Xwmc238', 'Xgwm149'], '4B'); r['evidence'] = 'Abstract; markers Xwmc238, Xgwm149'; N.append(r)
    r = mk('S151-002', name='Wuhan-1 2D FHB allele (Canadian backgrounds)', disease='Fusarium head blight', chr='2D', source='McCartney2007MB', marker='Xwmc245, Xgwm608',
           note=MC + 'The Wuhan-1 2D allele was an effective FHB resistance allele.'); place(r, ['Xwmc245', 'Xgwm608'], '2D'); r['evidence'] = 'Abstract; markers Xwmc245, Xgwm608'; N.append(r)
    r = mk('S151-003', name='Nyuubai and Sumai-3 3BSc FHB alleles (Canadian backgrounds)', disease='Fusarium head blight', chr='3B', source='McCartney2007MB', marker='Xgwm566, Xwmc231 (also Xwmc625, Xwmc693, Xwmc307, Xwmc418)',
           note=MC + 'The Nyuubai and Sumai-3 alleles at the centromeric 3BS QTL (3BSc) were effective FHB resistance alleles. The other 3BSc markers place far apart on 3B, so the pericentromeric pair Xgwm566/Xwmc231 is used.'); place(r, ['Xgwm566', 'Xwmc231'], '3B'); r['evidence'] = 'Abstract; 3BSc markers'; N.append(r)
    r = mk('S151-004', name='Sumai-3 and Nyuubai 5AS FHB alleles (Canadian backgrounds)', disease='Fusarium head blight', chr='5A', source='McCartney2007MB', marker='Xgwm154, Xwmc705, Xgwm304',
           note=MC + 'The Nyuubai and Sumai-3 5AS alleles were effective; the Sumai-3 5AS allele was significantly associated with lower grain protein content (same trend for Nyuubai, not significant).'); place(r, ['Xgwm154', 'Xwmc705', 'Xgwm304'], '5A'); r['evidence'] = 'Abstract; 5AS markers'; N.append(r)
    N.append(mk('S151-005', name='Nyuubai and Sumai-3 3BS (Fhb1) alleles (Canadian backgrounds)', disease='Fusarium head blight', chr='3B', source='McCartney2007MB', peak=FHB1, grade='C', marker='Xgwm533, Xgwm493',
        note=MC + 'The Nyuubai and Sumai-3 3BS (Fhb1) alleles were the least effective FHB resistance alleles in the nursery tests. Placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Abstract'))
    MI = 'Donor QTL from CM82036 (Sumai 3/Thornbird; 3B, 5A) and Frontana (3A) introgressed and stacked by marker selection in German elite spring wheat; eight QTL classes tested in F5 at four locations with F. culmorum spray inoculation. All three donor alleles significantly reduced DON content, and all but the Frontana 3A allele (low, non-significant effect on FHB severity) reduced FHB severity; 3B and 5A together reduced DON by 78% and FHB rating by 55%. '
    N.append(mk('S152-001', name='Fhb1 (CM82036; stacked in German spring wheat)', disease='Fusarium head blight', chr='3B', source='Miedaner2006TAG', peak=FHB1, grade='C', marker='WMS533',
        note=MI + 'Placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Abstract; Table 2'))
    N.append(mk('S152-002', name='Fhb5 / Qfhs.ifa-5A (CM82036; stacked in German spring wheat)', disease='Fusarium head blight', chr='5A', source='Miedaner2006TAG', start=FHB5[0], end=FHB5[1], grade='C', marker='WMS156, WMS293',
        note=MI + 'WMS156 places on distal 5AL, away from the pericentromeric Fhb5; placed at the Fhb5 interval in the atlas (105.4-214.2 Mb).', evidence='Abstract; Table 2'))
    N.append(mk('S152-003', name='Frontana 3A FHB QTL (stacked in German spring wheat)', disease='Fusarium head blight', chr='3A', source='Miedaner2006TAG', grade='D', marker='WMS720',
        note=MI + 'The Frontana 3A QTL maps near WMS720, which could not be placed.', evidence='Abstract'))
    pl = [(m,) + resolve(m, '7D') for m in ('Xcfd175', 'Xwmc14')]
    N.append(mk('S153-001', name='Pch1 (German DH populations and cultivars)', cls='DG', disease='Eyespot', chr='7D', source='Meyer2011Euph', start=min(x[1] for x in pl), end=max(x[1] for x in pl), grade='B',
        marker='Xorw1, Xorw5, Xorw6, Xcfd175, Xbarc76, Xwmc14, Xcfa2040',
        note='In 127 DH lines from three populations, Xorw1, Xorw5, Xorw6, Xcfd175, Xbarc76, Xwmc14 and Xcfa2040 were linked to Pch1 (Ep-D1) on 7DL; in 104 German winter wheat cultivars Xorw1, Xorw6 and Xcfd175 were the best markers for MAS, whereas XustSSR2001-7DL was population-specific. '
             'Placed between ' + ' and '.join(f'{m} ({p:.3f} Mb, {h})' for m, p, h in pl) + '.', evidence='Results: linkage in 127 DH lines; validation on 104 cultivars'))
    N.append(mk('S154-001', name='Fhb1 / Qfhs.ndsu-3BS (multiplex PCR screen)', disease='Fusarium head blight', chr='3B', source='MatusCadiz2006', peak=FHB1, grade='C', marker='Xgwm493, Xgwm533',
        note='A high-throughput marker-assisted selection method for Qfhs.ndsu-3BS (Fhb1) using 96-well DNA extraction, multiplexed gwm493 and gwm533 and agarose gels, verified on six resistant (Sumai-3, ND2710, McVey, BacUp, HY644, Alsen) and 52 susceptible parents and used to screen 5,567 breeding lines from Sumai-3 derivatives. Placed at the cloned Fhb1 (TaHRC) position in the atlas.',
        evidence='Abstract'))
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
    json.dump(ops, open('/tmp/ops34.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], (r['asm'].get('RefSeq v2.1') or {}).get('flag'))
    print(len(ops), len(json.dumps(ops)))
