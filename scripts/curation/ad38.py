"""ad38 (Nov 2026). Not added: Wang et al. preprint (Elymus repens 1StL telosome addition), Zheng et al. 2021 (meta-QTL by others), Cardozo Tellez et al. (germplasm screening), Zhong et al. 2021 (pathogen study), Wilhelm et al. 2013 (Rht haplotypes), Yao et al. 2019 (HL2 head length), Zhang et al. 2022 (already in atlas), Zhao et al. 2018 (UGT functional study)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Tessmann2018Agron': ('Tessmann EW, Van Sanford DA (2018) GWAS for Fusarium head blight related traits in winter wheat (Triticum aestivum L.) in an artificially warmed treatment. Agronomy 8:68', 'https://doi.org/10.3390/agronomy8050068', 'Tessmann and Van Sanford (2018)'),
 'Tessmann2019CS': ('Tessmann EW, Dong Y, Van Sanford DA (2019) GWAS for Fusarium head blight traits in a soft red winter wheat mapping panel. Crop Sci 59:1823-1837', 'https://doi.org/10.2135/cropsci2018.08.0492', 'Tessmann et al. (2019)'),
 'Tadesse2010JAG': ('Tadesse W, Schmolke M, Hsam SLK, Mohler V, Wenzel G, Zeller FJ (2010) Chromosomal location and molecular mapping of a tan spot resistance gene in the winter wheat cultivar Red Chief. J Appl Genet 51:235-242', 'https://doi.org/10.1007/BF03208853', 'Tadesse et al. (2010)'),
 'Waldron1999CS': ('Waldron BL, Moreno-Sevilla B, Anderson JA, Stack RW, Frohberg RC (1999) RFLP mapping of QTL for Fusarium head blight resistance in wheat. Crop Sci 39:805-811', 'https://doi.org/10.2135/cropsci1999.0011183X003900030032x', 'Waldron et al. (1999)'),
 'Tidakbi2026G3': ('Tidakbi L, Wang H, Bernardo A, St. Amand P, Fritz A, Zhang G, Bruce M, Rupp Noller JL, Brown-Guedira G, Bai G, Fellers J, Guttieri MJ, Jordan KW (2026) A winter wheat genome-wide association study identifies a large effect Fusarium head blight resistance haplotype using diverse lines from US breeding programs. G3 jkag291', 'https://doi.org/10.1093/g3journal/jkag291', 'Tidakbi et al. (2026)'),
 'Yoshioka2024Phyto': ('Yoshioka M, Kishii M, Singh PK, Inoue Y, Vy TTP, Tosa Y, Asuke S (2024) Rmg10, a novel wheat blast resistance gene derived from Aegilops tauschii. Phytopathology 114(9):2113', 'https://doi.org/10.1094/PHYTO-01-24-0018-R', 'Yoshioka et al. (2024)'),
 'Tagle2015Phyto': ('Tagle AG, Chuma I, Tosa Y (2015) Rmg7, a new gene for resistance to Triticum isolates of Pyricularia oryzae identified in tetraploid wheat. Phytopathology 105:495-499', 'https://doi.org/10.1094/PHYTO-06-14-0182-R', 'Tagle et al. (2015)'),
 'Wang2024PB': ('Wang X, Li G, Jia H, Cheng R, Shi J, Chen R, Wen Y, Ma Z et al. (2024) Breeding evaluation and precise mapping of Fhb8 for Fusarium head blight resistance in wheat (Triticum aestivum). Plant Breed 143:26-33', 'https://doi.org/10.1111/pbr.13113', 'Wang et al. (2024)'),
 'Zhang2011TAG': ('Zhang X, Shen X, Hao Y, Cai J, Ohm HW, Kong L (2011) A genetic map of Lophopyrum ponticum chromosome 7E, harboring resistance genes to Fusarium head blight and leaf rust. Theor Appl Genet 122:263-270', 'https://doi.org/10.1007/s00122-010-1441-3', 'Zhang et al. (2011)'),
 'Zanke2017TAG': ('Zanke CD, Rodemann B, Ling J, Muqaddasi QH, Plieske J, Polley A, Kollers S, Ebmeyer E, Korzun V, Argillier O, Stiewe G, Zschackel T, Ganal MW, Roder MS (2017) Genome-wide association mapping of resistance to eyespot disease (Pseudocercosporella herpotrichoides) in European winter wheat (Triticum aestivum L.) and fine-mapping of Pch1. Theor Appl Genet 130:505-514', 'https://doi.org/10.1007/s00122-016-2830-z', 'Zanke et al. (2017)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad38.py).'
FHB1, RMG7, PCH1 = 8.528, 761.8578, (631.678, 635.512)
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
def build(reports):
    keep = [r for r in reports if not re.match(r'S1(7[6-9]|8[0-5])-', r['id'])]; N = []
    TA = 'GWAS (90K array, 3,919 SNPs, GAPIT with known QTL as covariates) in 238 soft red winter wheat lines under control and soil-warmed (+3 C) treatments, Lexington KY, 2017. Markers are panel-internal SNP IDs; positions as printed (unit and assembly not stated), so not placed. '
    for n, (nm, ch, ms, nt) in enumerate([
        ('1B (243.6-255.1)', '1B', 'M1528, M1563, M1591', 'Control: M1528 (243.59) DON effect +0.87%; M1563 (250.98) and M1591 (255.13) FDK effects +2.11% and +1.99%.'),
        ('1B (382.3)', '1B', 'M1905', 'Control: FDK effect +2.40% (P 8.7e-4).'),
        ('2A (174.0)', '2A', 'M2508', 'Control: FDK effect -2.32% (P 5.2e-4).'),
        ('3B (285.3)', '3B', 'M6241', 'Control: FDK effect -1.82% (P 7.2e-4).'),
        ('1A (256.0)', '1A', 'M480', 'Control: FDK effect +1.72% (P 8.0e-4).'),
        ('5A (311.8)', '5A', 'M8313', 'Control: FDK effect -1.73% (P 9.4e-4).'),
        ('4A (356.8)', '4A', 'M6959', 'Control: DON effect -0.74% (P 1.8e-4); also found by Tessmann et al. (2019).'),
        ('6B (375.2)', '6B', 'M11423', 'Control: DON effect -0.96% (P 3.4e-4); also found by Tessmann et al. (2019).'),
        ('6A (97.5)', '6A', 'M9821', 'Control: DON effect +0.63% (P 9.0e-4).'),
        ('6B (226.6)', '6B', 'M11046', 'Control: DON effect +0.64% (P 9.4e-4).'),
        ('3B (218.7-224.5), warmed', '3B', 'M5743, M5744', 'Warmed: DON effects +0.94% and +1.13% (P 7.4e-4, 1.6e-4).'),
        ('4A (564.6), warmed', '4A', 'M7150', 'Warmed: DON effect +0.90% (P 1.0e-3).')], 1):
        N.append(mk(f'S176-{n:03d}', name=f'FHB MTA {nm} (Tessmann 2018)', disease='Fusarium head blight', chr=ch, source='Tessmann2018Agron', grade='D', marker=ms, note=TA + nt, evidence='Table 6'))
    TB = 'GWAS (90K array) in 250 soft red winter wheat lines of the TCAP panel, Lexington KY, 2015-2016, for FHB rating, severity, incidence, index, FDK and DON. Markers are panel-internal SNP IDs; positions as printed, so not placed. '
    for n, (nm, ch, ms, nt) in enumerate([
        ('7B (360.0-379.9)', '7B', 'M12955, M12957, M12960', 'Pleiotropic cluster for FHB rating, incidence and index (average): incidence effects +3.10% to +5.02% (P down to 8e-5).'),
        ('5B (433.1)', '5B', 'M9385', 'Pleiotropic (2015): severity -5.04% (P 1.0e-4), index -5.30%, FDK -1.09%.'),
        ('3A (284.2)', '3A', 'M5248', 'Pleiotropic (2015): severity -3.93% (P 2.9e-4) and index -4.26%.'),
        ('2B (240.2)', '2B', 'M3268', 'Pleiotropic (2015): incidence -3.91% (P 2.0e-4) and index -3.91%.'),
        ('5B (447.7)', '5B', 'M9432', 'Pleiotropic: DON -1.42% (average, P 2e-5) and -1.47% (2016); FDK -0.64% (average).'),
        ('4A (356.8)', '4A', 'M6959', 'DON -1.17% (average, P 8.2e-4); also found by Tessmann and Van Sanford (2018).'),
        ('6B (375.2)', '6B', 'M11423', 'FDK -0.86% (average, P 8.4e-4); also found by Tessmann and Van Sanford (2018).')], 1):
        N.append(mk(f'S177-{n:03d}', name=f'FHB MTA {nm} (Tessmann 2019)', disease='Fusarium head blight', chr=ch, source='Tessmann2019CS', grade='D', marker=ms, note=TB + nt, evidence='Table 3; Table 5'))
    r = mk('S178-001', name='tsr4 (Red Chief)', cls='DG', disease='Tan spot', chr='3A', source='Tadesse2010JAG', marker='Xgwm2a',
           note='Monosomic analysis (21 Chinese Spring monosomics x Red Chief) showed tan spot resistance in Red Chief is a single recessive gene on 3A; linkage analysis in Red Chief/Chinese Spring F2 placed tsr4 around Xgwm2a on 3AS. The Xgwm2a 130-bp allele was found only in Red Chief and Dashen (the 120-bp allele marks tsr3 on 3D in synthetic lines).', evidence='Abstract; Results')
    place(r, ['Xgwm2'], '3A'); N.append(r)
    WA = 'Type II FHB resistance (greenhouse point inoculation, two experiments) in 112 F5-derived RILs of Sumai 3 x Stoa; RFLP map, single-marker and interval analysis. '
    N.append(mk('S179-001', name='QFhs.ndsu-3B (Sumai 3; first report of Fhb1)', disease='Fusarium head blight', chr='3B', source='Waldron1999CS', peak=FHB1, grade='C', marker='Xcdo981, Xbcd907', pve='15.4',
        note=WA + '3BS QTL from Sumai 3 (LOD over 3.0), designated QFhs.ndsu-3B; best marker Xcdo981 explained 15.4% (Xbcd907 recommended as more reliable). RILs with the 3BS and 2AL QTL had median severity 20.9% vs 36.2% overall. RFLP markers not placed; placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Table 2'))
    for n, (nm, ch, m, org, pv, nt) in enumerate([('QFhs.ndsu-2A', '2A', 'XksuH16', 'Stoa', '14.3', '2AL QTL (LOD over 3.0).'),
                                                 ('4BL FHB region (Stoa)', '4B', 'Xwg909', 'Stoa', '7.2', '4BL region (P 0.007), below the interval-mapping threshold.'),
                                                 ('6BS FHB region (Sumai 3)', '6B', 'Xbcd331, Xcdo524', 'Sumai 3', '6.0', '6BS region (P 0.010 and 0.049), below the interval-mapping threshold.')], 2):
        N.append(mk(f'S179-{n:03d}', name=nm + ' (Sumai 3 x Stoa)', disease='Fusarium head blight', chr=ch, source='Waldron1999CS', grade='D', marker=m, pve=pv,
            note=WA + f'Resistance from {org}. {nt} RFLP marker, not placed.', evidence='Table 2'))
    TI = 'GWAS (exome-capture SNPs, CS v2.1) in 269 US winter wheat lines evaluated for FHB incidence, severity, AUDPC, FDK and DON in the Rocky Ford nursery (Manhattan, KS) over 3 seasons; 22 prioritized loci (Table 3). '
    N.append(mk('S180-001', name='Jamestown 1BL FHB resistance haplotype', disease='Fusarium head blight', chr='1B', source='Tidakbi2026G3', ref='RefSeq v2.1', start=166.444551, end=435.466844, grade='A', marker='1B SNPs at 166.44, 332.26 and 435.47 Mb; KASP assays',
        note=TI + 'Three 1B loci (166.44 Mb, near TraesCS1B03G0350400/0350500; 332.26 Mb, missense in TraesCS1B03G0534700; 435.47 Mb, missense in TraesCS1B03G0687800) affecting AUDPC, FDK and DON form a haplotype carried by Jamestown and 21 other lines, reducing AUDPC and FDK by about 25% and DON by 18% (P 1.9e-7 or less); corresponds to the Jamestown native QTL of Carpenter et al. (2020). A chromosome region, not a locus.', evidence='Table 3; Fig. 4'))
    L = [('1B', 458.665282, 'FDK', 'missense, TraesCS1B03G0724900'), ('1B', 569.16332, 'DON', 'upstream of TraesCS1B03G0919200'), ('2A', 28.823494, 'severity', 'intron/splice region, TraesCS2A03G0111400'),
         ('2A', 665.979623, 'severity, AUDPC', 'intron, TraesCS2A03G0976100'), ('2B', 98.006278, 'DON', 'intron, TraesCS2B03G0298900'), ('2B', 162.824813, 'severity, AUDPC', "5' UTR, TraesCS2B03G0437700"),
         ('2B', 580.030865, 'AUDPC', 'missense, TraesCS2B03G1028200'), ('2B', 782.993929, 'FDK', "3' UTR, TraesCS2B03G1470300"), ('3A', 11.532385, 'FDK', 'downstream of TraesCS3A03G00368'),
         ('3B', 843.038669, 'DON', "3' UTR, TraesCS3B03G1503100"), ('5B', 509.272985, 'severity, AUDPC', "3' UTR, TraesCS5B03G0810500"), ('5B', 701.047429, 'severity, AUDPC', 'intron, TraesCS5B03G1326300'),
         ('5D', 550.632004, 'severity, AUDPC', 'intron, TraesCS5D03G1160400'), ('6A', 18.452189, 'severity', "3' UTR, TraesCS6A03G0074800"),
         ('6A', 553.038439, 'DON, FDK', 'intron, TraesCS6A03G0818700; within the NC-Neuse QTL region, NC-Neuse allele lowered DON by 9% and FDK by 14%'),
         ('7A', 18.049084, 'DON', 'missense, TraesCS7A03G0083000'), ('7A', 21.559496, 'severity, AUDPC', 'missense, TraesCS7A03G0101600'), ('7A', 26.275779, 'severity', 'missense, TraesCS7A03G0118400'),
         ('7A', 612.449697, 'severity, AUDPC', 'synonymous, TraesCS7A03G1010300')]
    for n, (ch, p, tr, g) in enumerate(L, 2):
        N.append(mk(f'S180-{n:03d}', name=f'FHB GWAS locus {ch} {p:.2f} Mb (US winter wheat)', disease='Fusarium head blight', chr=ch, source='Tidakbi2026G3', ref='RefSeq v2.1', peak=p, grade='A', marker=f'SNP {ch}:{round(p*1e6)}',
            note=TI + f'Prioritized SNP at {p:.6f} Mb (CS v2.1) for {tr} ({g}).', evidence='Table 3'))
    N.append(mk('S181-001', name='Rmg10 (Ae. tauschii KU-2097)', cls='DG', disease='Wheat blast', chr='2D', source='Yoshioka2024Phyto', ref='RefSeq v2.1', start=6.76, end=24.73, grade='C', marker='MS-4, MS-12',
        note='Of 199 Ae. tauschii accessions screened with Brazilian MoT isolate Br48, 14 were resistant; synthetic Ldn/KU-2097 was resistant at seedling and spike stages, and a single dominant gene, Rmg10, was mapped in Chinese Spring x Ldn/KU-2097 (86 F2:3 lines) to a 16.1-cM interval on 2DS between MS-4 (7.8 cM distal) and MS-12 (8.3 cM proximal), about 5.1 Mb in CS v2.1. Effective against Brazilian but not Bangladesh isolates; a Gladius/KU-2097 synthetic octoploid had lower blast index in Bolivian field tests. Placed at the mapped 2DS region shown in Fig. 5 (6.76-24.73 Mb, v2.1); the exact 5.1-Mb interval is not printed.', evidence='Results; Fig. 5'))
    N.append(mk('S182-001', name='Rmg7 (T. dicoccum St24; original report)', cls='DG', disease='Wheat blast', chr='2A', source='Tagle2015Phyto', peak=RMG7, grade='C', marker='',
        note='A single gene, Rmg7, identified in tetraploid T. dicoccum St24 (KU120) against Triticum isolate Br48; St17 and St25 carry the same gene (no susceptible F2 from intercrosses), and a single corresponding avirulence gene was detected in Br48 x MZ5-1-6 progeny. No map position in this paper; placed at the Rmg7 position in the atlas (2AL, allele of Pm4; Asuke et al. 2024).', evidence='Abstract'))
    N.append(mk('S183-001', name='Fhb8 / Qfdk.nau-7D (Wangshuibai)', cls='DG', disease='Fusarium head blight', chr='7D', source='Wang2024PB', ref='RefSeq v2.1', start=93.9, end=96.5, grade='A', marker='Xwgrb1500, Xwgrb1559, Xwgrb1587 (co-segregating)',
        note='Qfdk.nau-7D from Wangshuibai (FDK) was transferred into susceptible PH691 by marker-assisted backcrossing; NILs showed it lowers FDK through resistance to infection and spread, and differ from PH691 in head length. Fine-mapped in a BC3F2 population (97 plants) with recombinant testing to a 1.0-cM interval Xwgrb1500-Xwgrb1559 (93.9-96.5 Mb in CS v2.1), co-segregating with Xwgrb1587, and designated Fhb8.', evidence='Abstract; Results'))
    N.append(mk('S184-001', name='FhbLoP (L. ponticum 7EL)', disease='Fusarium head blight', chr='7E', source='Zhang2011TAG', grade='D', ref='alien (L. ponticum 7EL)', marker='Xcfa2240, Xswes19', pve='30.46',
        note='A 7E genetic map (64 markers, 95.76 cM) from 237 F7:8 RILs of two Thatcher-L. ponticum substitution lines, K11463 (7el1(7D)) x K2620 (7el2(7D)); the major FHB resistance locus FhbLoP mapped to the distal 7EL in a 3.71-cM interval Xcfa2240-Xswes19, explaining 30.46% of the variance. Alien chromosome: no wheat RefSeq position.', evidence='Abstract'))
    ZA = 'GWAS for eyespot (field tests in three environments) in 168 European winter and 3 spring wheat varieties with 732 SSRs and 7,761 mapped 90K SNPs; 108 SSR and 235 SNP MTAs (-log10P 3.0 or more). '
    N.append(mk('S185-001', name='Pch1 (European winter wheat; fine-mapping)', cls='DG', disease='Eyespot', chr='7D', source='Zanke2017TAG', start=PCH1[0], end=PCH1[1], grade='C', marker='Xorw1, AX-94525717',
        note=ZA + 'Pch1 (T. ventricosum introgression on 7DL) was the main resistance source, present in 18 varieties (eyespot score 1.20-2.67 vs 3.43-5.73 without). Varieties with recombination in the introgression allowed fine-mapping to the distal 7DL region Traes_7DL_59D114555-Traes_4D32A4BE1, which contains the CC-NBS-LRR gene Traes_7DL_973A33763 (marker AX-94525717), a candidate for Pch1. Placed at the Pch1 interval in the atlas.', evidence='Table 2; Results'))
    r = mk('S185-002', name='Eyespot MTA 5A (XBARC303; European winter wheat)', disease='Eyespot', chr='5A', source='Zanke2017TAG', marker='Xbarc303',
           note=ZA + 'Besides the 7D markers, BARC303 (12.6 cM on 5A) was the only marker beyond the Bonferroni threshold; distant from the WMS639 QTL of Burt et al. (2011).', evidence='Discussion'); place(r, ['Xbarc303'], '5A'); N.append(r)
    N.append(mk('S185-003', name='Eyespot QTL region 2AS (European winter wheat)', disease='Eyespot', chr='2A', source='Zanke2017TAG', grade='D',
        note=ZA + 'SSR and SNP data indicated a QTL region on 2AS supported by numerous markers, most alleles reducing eyespot score (smaller effects than Pch1); proposed as a novel resistance source. Markers listed only in the online resources (not uploaded); not placed.', evidence='Discussion'))
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
    json.dump(ops, open('/tmp/ops38.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], (r['asm'].get('RefSeq v2.1') or {}).get('flag'))
    print(len(ops), len(json.dumps(ops)))
