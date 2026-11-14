"""ad29 (Nov 2026): ten uploads. Added: Jiang et al. 2015 Heredity 114:318 (Rht-D1 and Ppd-D1 associations with FHB), Islam et al. 2024
Phytopathology 114:1878 (Rmg11), Zhang et al. 2025 preprint (leaf rust APR QTL, Xinong1163-4 x Thatcher), Jin et al. 2020 TAG 133:2063
(crown rot 5DL region), Kage et al. 2018 Eur J Plant Pathol 151:125 (QTL-2DL candidate genes).
Not added: Joshi & Chand 2002 Euphytica 124:283 and Joshi et al. 2007 Euphytica 153:59 (trait inheritance and association, no mapping),
J Bangladesh Agric Univ 2019 17:161 (in silico domain survey, no locus), the TAG 2017 130:471 study (prediction accuracy of MAS/GS,
no loci reported), Hu et al. 2019 BMC Genet 20:76 (copyright correction; the original article is already in the atlas as V_Hu2019b)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Jiang2015Hered': ('Jiang Y, Zhao Y, Rodemann B, Plieske J, Kollers S, Korzun V, Ebmeyer E, Argillier O, Hinze M, Ling J, Roder MS, Ganal MW, Mette MF, Reif JC (2015) Potential and limits to unravel the genetic architecture and predict the variation of Fusarium head blight resistance in European winter wheat (Triticum aestivum L.). Heredity 114:318-326', 'https://doi.org/10.1038/hdy.2014.104', 'Jiang et al. (2015)'),
 'Islam2024Rmg11': ('Islam MT, Nago C, Yoshioka M, Vy TTP, Tosa Y, Asuke S (2024) Identification of Rmg11 in tetraploid wheat as a new blast resistance gene with tolerance to high temperature. Phytopathology 114:1878-1883', 'https://doi.org/10.1094/PHYTO-02-24-0074-R', 'Islam et al. (2024)'),
 'Zhang2025pre': ('Zhang JQ, Kang ZH, Li X, Li M, Xue L, Li X (2025) QTL mapping of adult plant resistance to wheat leaf rust in the Xinong1163-4 x Thatcher RIL population. Preprints.org preprint', 'https://doi.org/10.20944/preprints202505.0686.v1', 'Zhang et al. (2025)'),
 'Jin2020TAG': ('Jin J, Duan S, Qi Y, Yan S, Li W, Li B, Xie C, Zhen W et al. (2020) Identification of a novel genomic region associated with resistance to Fusarium crown rot in wheat. Theor Appl Genet 133:2063-2073', 'https://doi.org/10.1007/s00122-020-03577-1', 'Jin et al. (2020)'),
 'Kage2018EJPP': ('Kage U, Hukkeri S, Kushalappa AC (2018) Liquid chromatography and high resolution mass spectrometry-based metabolomics to identify quantitative resistance-related metabolites and genes in wheat QTL-2DL against Fusarium head blight. Eur J Plant Pathol 151:125-139', 'https://doi.org/10.1007/s10658-017-1362-y', 'Kage et al. (2018)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad29.py).'
RHTD1, PPDD1 = (18.781062, 18.782933), (33.952048, 33.956269)
def mk(i, **k):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand'); r.update(k); r['note'] += T; return r
def place(r, ms, ch):
    pl = [(m,) + resolve(m, ch) for m in ms]; pl = [x for x in pl if x[1]]
    if len(pl) >= 2:
        lo, hi = min(x[1] for x in pl), max(x[1] for x in pl); r['start'], r['end'] = lo, hi; r['grade'] = 'B' if hi - lo <= 100 else 'C'
        t = ' Placed in silico: ' + ', '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in pl) + '.' + (f' The markers span {hi - lo:.0f} Mb: a chromosome region, not a locus.' if r['grade'] == 'C' else '')
    elif pl: m, p, h = pl[0]; r['peak'], r['grade'] = p, 'C'; t = f' Only {m} could be placed ({p:.3f} Mb, {h}); position is that marker.'
    else: r['grade'] = 'D'; t = ' Neither flanking marker could be placed.'
    r['note'] = r['note'][:-len(T)] + t + T
def build(reports):
    keep = [r for r in reports if not re.match(r'S1(3[6-9]|40)-', r['id'])]; N = []
    JI = 'Association mapping in 372 European winter wheat varieties (FHB scores from multi-environment trials; SSR, 9k and 90k SNP data): the functional marker was significantly (P < 0.005) associated with FHB resistance, and excluding it and linked markers lowered genomic prediction accuracy by 14-30%. '
    N.append(mk('S136-001', name='Rht-D1 (FHB; European winter wheat varieties)', cls='MTA', disease='Fusarium head blight', chr='4D', source='Jiang2015Hered', start=RHTD1[0], end=RHTD1[1], grade='C', marker='Rht-D1 functional marker',
        note=JI + 'Rht-D1b is associated with higher FHB; Rht-B1 was not significant (mutant allele in only 26 varieties). Placed at the Rht-D1 gene (TraesCS4D02G040400).', evidence='Results; Supplementary Table S2'))
    N.append(mk('S136-002', name='Ppd-D1 (FHB; European winter wheat varieties)', cls='MTA', disease='Fusarium head blight', chr='2D', source='Jiang2015Hered', start=PPDD1[0], end=PPDD1[1], grade='C', marker='Ppd-D1 functional marker',
        note=JI + 'The association reflects the link between heading time and FHB. Placed at the Ppd-D1 gene (TraesCS2D02G079600).', evidence='Results; Supplementary Table S2'))
    p, h = resolve('Xgwm635', '7A')
    N.append(mk('S137-001', name='Rmg11 (St19, T. dicoccum)', cls='DG', disease='Wheat blast', chr='7A', source='Islam2024Rmg11', ref='RefSeq v2.1', start=5.3, end=5.7, peak=5.5, grade='C', marker='IMT_5-IMT_6/IMT_7 (Xgwm635)',
        note='New blast resistance gene Rmg11 in tetraploid accession St19 (T. dicoccum KU-114), effective against Br48 and its AVR-Rmg8 disruptant even at 30 C and against Bangladeshi isolates at the seedling stage; intermediate head resistance to Brazilian isolates. Single dominant gene (91 F2:3 lines of St19 x St30) on 7AS between IMT_5 (1.0 cM distal) and IMT_6/IMT_7 (1.1 cM proximal), 2.2 cM distal of Xgwm635. '
             f'Fig. 2 places the interval at about 5.3-5.7 Mb on durum Svevo v1.0; that position is used directly as an approximate RefSeq v2.1 position (Xgwm635 places at {p:.1f} Mb on RefSeq v1.0).', evidence='Fig. 2: genetic map of 7AS and Svevo positions'))
    ZP = 'Preprint (not peer reviewed). Leaf rust adult-plant resistance QTL (maximum disease severity) in 195 Xinong1163-4 x Thatcher RILs, Baoding 2018-2020 (BSA with SSR markers, ICIM); resistance from Xinong1163-4. '
    for n, (nm, ch, ms, det, mx, ex) in enumerate([('1BL.1', '1B', ['Xgwm140', 'Xgpw1077'], '2018 LOD 4.20 (10.7%)', '10.7', ' The 1BL QTL are predicted to correspond to Lr46/Yr29.'),
            ('1BL.2', '1B', ['Xgpw1077', 'Xcfa2147.1'], '2018 LOD 2.82 (3.4%), 2020 3.79 (5.4%)', '5.4', ''), ('1BL.3', '1B', ['Xgwm728', 'Xgwm140'], '2019 LOD 6.85 (21.9%), 2020 2.68 (6.2%)', '21.9', ' Shares flanking markers with Lr46.'),
            ('3A', '3A', ['Xwmc215', 'Xcfa2076'], '2019 LOD 9.57 (25.0%), 2020 3.13 (6.6%)', '25.0', ''), ('4BL', '4B', ['Xwmc657', 'Xgwm251'], '2018 LOD 4.13 (8.1%), 2019 2.95 (5.3%), 2020 6.03 (12.2%)', '12.2', ' Detected in all three seasons; region reported for leaf rust, stripe rust, FHB and powdery mildew.'),
            ('4BL.1', '4B', ['Xgwm193', 'Xwmc692'], '2018 LOD 3.46 (12.2%)', '12.2', '')], 1):
        r = mk(f'S138-{n:03d}', name=f'QLr.hbau-{nm}', disease='Leaf rust', chr=ch, source='Zhang2025pre', marker='-'.join(ms), pve=f'up to {mx}',
               note=ZP + f'Detections: {det}.' + ex, evidence=f'Table 3: QLr.hbau-{nm} {"-".join(ms)}')
        place(r, [m.split('.')[0] for m in ms], ch); N.append(r)
    N.append(mk('S139-001', name='QFcr 5DL region (Chinese germplasm GWAS)', cls='MTA', disease='Fusarium crown rot', chr='5D', source='Jin2020TAG', start=218.636512, end=232.419393, grade='A',
        marker='Affx-110560650, Affx-110363452, Affx-111302661, Affx-109780353, Affx-110167226',
        note='GWAS of seedling crown rot resistance in 358 Chinese wheat germplasms (55K SNP array, three trials): 104 SNPs on 1BS, 1DS, 2AL, 5AL, 5DS, 5DL, 6BS and 7BL were significant (P < 0.001; listed only in Table S6, not uploaded); five SNPs on 5DL detected in all three trials define a novel 13.78 Mb region (218.64-232.42 Mb), confirmed in a biparental population; TIR-NBS-LRR and DON-detoxification genes in the region were induced faster in resistant 04 Zhong 36 than in susceptible Xinmai 26. Positions as printed (RefSeq v1.0).',
        evidence='Table 3: five 5DL SNPs 218,636,512-232,419,393 bp'))
    r = mk('S140-001', name='QTL-2DL candidate genes GPAT, CER5, PAL (Wuhan 1 NILs)', disease='Fusarium head blight', chr='2D', source='Kage2018EJPP', marker='Xwmc245, Xgpw8003, Xgwm539, Xgwm608',
           note='Metabolomics of near-isogenic lines for the Wuhan 1 QTL-2DL (spikelet inoculation): phenylpropanoid, fatty acid and glycerophospholipid metabolites were induced in NIL-R, and genes for glycerol-3-phosphate acyltransferase (GPAT), the ABC transporter CER5 and phenylalanine ammonia-lyase (PAL) located in the QTL-2DL interval (defined by wmc245, gpw8003, gwm539 and gwm608) were proposed as candidates reinforcing cell walls.')
    place(r, ['Xwmc245', 'Xgwm539'], '2D'); r['evidence'] = 'Results: candidate genes in QTL-2DL'; N.append(r)
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
    json.dump(ops, open('/tmp/ops29.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], r['asm'].get('RefSeq v2.1'))
