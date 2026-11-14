"""ad28 (Nov 2026): ten uploads. Added: Grausgruber et al. 1997 (substitution lines, chromosome level), Gyawali et al. 2023 Crops 3:28
(Fhb1 and Fhb7 pyramiding), Guo et al. 2003 TAG 106:1011 (3BS QTL, AFLP/STS), Hayden et al. 2004 TAG 109:1641 (Sr2 STM markers),
Gilsinger et al. 2005 TAG 110:1218 (FHB incidence and flower opening), Guo et al. 2015 TAG 128:2301 (Fhb7 high-density map), Herter et al.
2018 TAG 131:1263 (FHB QTL, Solitar x Bussard), Herter et al. 2019 TAG 132:1121 (GWAS hits used in genomic prediction), He et al. 2014
Euphytica 195:383 (haplotyping of Chinese germplasm). Not added: Hiebert et al. 2010 Mol Breed 26:681 (gene stacking with telocentric
chromosomes; no mapping)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Grausgruber1997CRC': ('Grausgruber H, Lemmens M, Buerstmayr H, Ruckenbauer P (1997) Genetic analysis of Fusarium head blight resistance and toxin tolerance in wheat using inter-varietal chromosome substitution lines. Cereal Res Commun 25:743-744', 'https://doi.org/10.1007/BF03543831', 'Grausgruber et al. (1997)'),
 'Gyawali2023Crops': ('Gyawali B, Scofield SR, Mohammadi M (2023) Marker development and pyramiding of Fhb1 and Fhb7 for enhanced resistance to Fusarium head blight in soft red winter wheat. Crops 3:28', 'https://doi.org/10.3390/crops3040028', 'Gyawali et al. (2023)'),
 'Guo2003TAG': ('Guo PG, Bai GH, Shaner GE (2003) AFLP and STS tagging of a major QTL for Fusarium head blight resistance in wheat. Theor Appl Genet 106:1011-1017', 'https://doi.org/10.1007/s00122-002-1129-4', 'Guo et al. (2003)'),
 'Hayden2004TAG': ('Hayden MJ, Kuchel H, Chalmers KJ (2004) Sequence tagged microsatellites for the Xgwm533 locus provide new diagnostic markers to select for the presence of stem rust resistance gene Sr2 in bread wheat (Triticum aestivum L.). Theor Appl Genet 109:1641-1647', 'https://doi.org/10.1007/s00122-004-1787-5', 'Hayden et al. (2004)'),
 'Gilsinger2005TAG': ('Gilsinger J, Kong L, Shen X, Ohm H (2005) DNA markers associated with low Fusarium head blight incidence and narrow flower opening in wheat. Theor Appl Genet 110:1218-1225', 'https://doi.org/10.1007/s00122-005-1953-4', 'Gilsinger et al. (2005)'),
 'Guo2015TAG': ('Guo J, Zhang X, Hou Y, Cai J, Shen X, Zhou T, Xu H, Ohm HW, Wang H, Li A, Han F, Wang H, Kong L (2015) High-density mapping of the major FHB resistance gene Fhb7 derived from Thinopyrum ponticum and its pyramiding with Fhb1 by marker-assisted selection. Theor Appl Genet 128:2301-2316', 'https://doi.org/10.1007/s00122-015-2586-x', 'Guo et al. (2015)'),
 'Herter2018TAG': ('Herter CP, Ebmeyer E, Kollers S, Korzun V, Leiser WL, Wurschum T, Miedaner T (2018) Rht24 reduces height in the winter wheat population Solitar x Bussard without adverse effects on Fusarium head blight infection. Theor Appl Genet 131:1263-1272', 'https://doi.org/10.1007/s00122-018-3076-8', 'Herter et al. (2018)'),
 'Herter2019TAG': ('Herter CP, Ebmeyer E, Kollers S, Korzun V, Wurschum T, Miedaner T (2019) Accuracy of within- and among-family genomic prediction for Fusarium head blight and Septoria tritici blotch in winter wheat. Theor Appl Genet 132:1121-1135', 'https://doi.org/10.1007/s00122-018-3264-6', 'Herter et al. (2019)'),
 'He2014Euph': ('He X, Singh PK, Schlang N, Duveiller E, Dreisigacker S, Payne T, He Z (2014) Characterization of Chinese wheat germplasm for resistance to Fusarium head blight at CIMMYT, Mexico. Euphytica 195:383-395', 'https://doi.org/10.1007/s10681-013-1002-3', 'He et al. (2014)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad28.py).'
FHB1, RHTD1, SR2 = 8.528, 18.781, 1.17
def mk(i, **k):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand'); r.update(k); r['note'] += T; return r
def two(r, a, b, ch):
    pa, ha = resolve(a, ch); pb, hb = resolve(b, ch)
    if pa and pb: r['start'], r['end'] = min(pa, pb), max(pa, pb); return f' Placed between {a} ({pa:.3f} Mb, {ha}) and {b} ({pb:.3f} Mb, {hb}).'
    p, m, h = (pa, a, ha) if pa else (pb, b, hb)
    if p: r['peak'] = p; return f' Only {m} could be placed ({p:.3f} Mb, {h}).'
    r['grade'] = 'D'; return ' Neither marker could be placed.'
def build(reports):
    keep = [r for r in reports if not re.match(r'S1(2[7-9]|3[0-5])-', r['id'])]; N = []
    G = 'Inter-varietal chromosome substitution lines tested for FHB in mist-irrigated field trials (Austria 1994-1996, F. graminearum and F. culmorum spray inoculation; SAUDPC and relative ear weight) and for toxin tolerance (Petritox test). '
    for n, (ch, don, bg, txt) in enumerate([('6B', 'T. macha', 'Hobbit sib', 'significantly higher resistance (SAUDPC) and higher DON and 3-AcDON tolerance'),
            ('4B', 'T. macha', 'Hobbit sib', 'significantly higher resistance (SAUDPC) and higher toxin tolerance'), ('3A', 'T. macha', 'Hobbit sib', 'significantly higher resistance (SAUDPC)'),
            ('5A', 'T. macha', 'Hobbit sib', 'tolerance: higher relative ear weight with no effect on symptoms, and higher toxin tolerance'),
            ('6D', 'T. spelta', 'Chinese Spring', 'significantly better FHB index'), ('7A', 'Lutescens 62', 'Chinese Spring', 'significantly better FHB index')], 1):
        N.append(mk(f'S127-{n:03d}', name=f'FHB chromosome effect {ch} ({don} in {bg})', disease='Fusarium head blight', chr=ch, source='Grausgruber1997CRC', grade='D',
            note=G + f'Chromosome {ch} of {don} substituted into {bg}: {txt}. Whole-chromosome effect: no position.', evidence=f'Table 1: {bg} [{don}] {ch}'))
    GY = 'Gene-based markers for Fhb1 (TaHRC-R/S) and Fhb7 (the GST ORF) screened in 57 soft red winter wheat breeding lines from crosses introducing Fhb7 into Fhb1 lines (Purdue); 15 lines carried both resistant alleles and two years of inoculated trials showed the Fhb1+Fhb7 combination lowered DON. '
    N.append(mk('S128-001', name='Fhb1 (SRW wheat, Fhb1+Fhb7 pyramids)', disease='Fusarium head blight', chr='3B', source='Gyawali2023Crops', peak=FHB1, grade='C', marker='TaHRC gene-based marker',
        note=GY + 'Placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Results: marker genotypes and DON of 57 lines'))
    N.append(mk('S128-002', name='Fhb7 (SRW wheat, Fhb1+Fhb7 pyramids)', disease='Fusarium head blight', chr='7D', ref='alien (7E/7DL translocation)', source='Gyawali2023Crops', grade='D', marker='GST (Fhb7) gene-based marker',
        note=GY + 'Fhb7 lies on a Thinopyrum segment translocated to 7DL and has no counterpart position in Chinese Spring.', evidence='Results: marker genotypes and DON of 57 lines'))
    N.append(mk('S129-001', name='Fhb1 (Ning 7840 x Clark; AFLP/STS)', disease='Fusarium head blight', chr='3B', source='Guo2003TAG', peak=FHB1, grade='C', marker='SRST.3B1 (STS from pAG/mCAA244), pAGT/mCTG57, pACT/mCTG136', pve='38-51',
        note='Major 3BS FHB QTL in Ning 7840/Clark RILs (type II resistance): five PstI-AFLPs associated, three explaining 38-51% of the variation; pAG/mCAA244 (in repulsion) was converted into the STS marker SRST.3B1, whose resistant band matched the pedigrees of Ning 7840 and Sumai 3 derived cultivars. Placed at the cloned Fhb1 (TaHRC) position in the atlas.',
        evidence='Results: R2 0.508, 0.459 and 0.38 for the three AFLPs'))
    N.append(mk('S130-001', name='Sr2 (STM markers at Xgwm533)', cls='DG', disease='Stem rust', chr='3B', source='Hayden2004TAG', ref='RefSeq v2.1', peak=SR2, grade='C', marker='Xgwm533 STMs (stm559tgag, stm598tcac)',
        note='Sequence-tagged microsatellites developed at the Xgwm533 locus discriminate the 120-bp allele linked to adult-plant stem rust resistance gene Sr2 from a size-homoplasious allele, and were diagnostic for Sr2 across a wide range of germplasm; the locus is also associated with a major FHB resistance gene (Fhb1). Placed at the Sr2 position in the atlas.',
        evidence='Results: STM genotypes of Sr2 carriers and non-carriers'))
    GI = 'Patterson x Goldfield RILs (100 lines; six environments, misted and non-misted field tests): narrow, brief flower opening of Goldfield associated with low FHB incidence (passive resistance). '
    r = mk('S131-001', name='QFhb.pur-2BS (low incidence, narrow flower opening)', disease='Fusarium head blight', chr='2B', source='Gilsinger2005TAG', grade='C', marker='Xbarc200-Xgwm210', pve='29', stat='LOD 5.21',
           note=GI + 'QTL in the Xbarc200-Xgwm210 interval on 2BS, 29% of the variation for incidence over six environments (LOD 5.21), also associated with narrow flower opening in three of five tests; Goldfield allele.')
    r['note'] = r['note'][:-len(T)] + two(r, 'Xbarc200', 'Xgwm210', '2B') + T; r['evidence'] = 'Table 1: Xbarc200, Xgwm210 2B'; N.append(r)
    p, h = resolve('Xwmc149', '2B')
    N.append(mk('S131-002', name='FHB incidence association at Xwmc149', cls='MTA', disease='Fusarium head blight', chr='2B', source='Gilsinger2005TAG', peak=p, grade='C', marker='Xwmc149', pve='12',
        note=GI + f'Xwmc149 (not linked to Xbarc200-Xgwm210) explained 12% of the variation for incidence. Placed at Xwmc149 ({p:.3f} Mb, {h}).', evidence='Table 1: Xwmc149 2B'))
    N.append(mk('S131-003', name='FHB incidence association at Xgwm344', cls='MTA', disease='Fusarium head blight', chr='7B', source='Gilsinger2005TAG', grade='D', marker='Xgwm344', pve='7',
        note=GI + 'Xgwm344 (7B) explained 7% of the variation for incidence and was not associated with flower opening. Xgwm344 could not be placed.', evidence='Table 1: Xgwm344 7B'))
    N.append(mk('S132-001', name='Fhb7 (high-density map; shortened translocations)', disease='Fusarium head blight', chr='7D', ref='alien (7E/7DL translocation)', source='Guo2015TAG', grade='D', marker='XsdauK66-Xcfa2240',
        note='Fhb7 from Th. ponticum mapped in a 1.7 cM interval flanked by XsdauK66 and Xcfa2240 (SSR, DArT, EST and COS markers) on the 7el2L arm of the 7DS.7el2L translocation in KS24-2; ph1b-induced secondary translocation lines SDAU1881 and SDAU1886 carry Fhb7 on shortened alien segments (about 16-17% of the chromosome), and lines with both Fhb1 and Fhb7 were developed by MAS. Fhb7 lies on an alien segment and has no counterpart position in Chinese Spring.',
        evidence='Results: Fhb7 between XsdauK66 and Xcfa2240'))
    HB = 'FHB severity QTL in 157 Solitar x Bussard RILs (four German locations, spray inoculation); the population is fixed for Rht-1 wild-type alleles and segregates for Rht24, which had no effect on FHB. Flanking markers are given only as codes (Supplementary Table 2, not uploaded), so the QTL cannot be placed. '
    for n, (ch, cm, l, rr, don, eff, pg, lod) in enumerate([('5B', '118 (116-124)', 'M2441', 'M2290', 'Solitar', '1.01', '13.86', '3.81'), ('2D', '92 (86-98)', 'M425', 'M615', 'Solitar', '1.31', '14.31', '4.10'),
                                                            ('7D', '140 (138-154)', 'M2491', 'M791', 'Bussard', '0.92', '12.14', '3.25')], 1):
        N.append(mk(f'S133-{n:03d}', name=f'QFhb {ch} (Solitar x Bussard)', disease='Fusarium head blight', chr=ch, source='Herter2018TAG', grade='D', marker=f'{l}-{rr} (coded)', pve=pg, stat=f'LOD {lod}',
            note=HB + f'QTL at {cm} cM, resistance from {don}, additive effect {eff}%, {pg}% of genotypic variance, LOD {lod}.', evidence=f'Table 3: FHB {ch} {cm} cM'))
    HE = 'GWAS within two winter wheat data sets (438 and 585 lines from six and eight bi-parental families) used to fit QTL as fixed effects in genomic prediction. '
    N.append(mk('S134-001', name='Rht-D1 (FHB susceptibility; winter wheat families)', cls='MTA', disease='Fusarium head blight', chr='4D', source='Herter2019TAG', peak=RHTD1, grade='C', marker='Rht-D1 diagnostic marker',
        note=HE + 'The Rht-D1 locus accounted for 62.3% of the genotypic variance for FHB severity in data set 1 (allele substitution effect -10.86%) and about 13% in data set 2; Rht-D1b increases susceptibility. Placed at the Rht-D1 gene (TraesCS4D01G040400; position entered from the atlas curator, Nov 2026).',
        evidence='Results; Supplementary Table 4'))
    N.append(mk('S134-002', name='Rht-D1 (STB susceptibility; winter wheat families)', cls='MTA', disease='Septoria tritici blotch', chr='4D', source='Herter2019TAG', peak=RHTD1, grade='C', marker='Rht-D1 diagnostic marker',
        note=HE + 'In data set 2 the Rht-D1 marker was also significantly associated with STB severity (about 13% of genotypic variance). Placed at the Rht-D1 gene (TraesCS4D01G040400).', evidence='Results'))
    p, h = resolve('wsnp_BE444720A_TA_2_1', '5A')
    N.append(mk('S134-003', name='Fhb5 (Qfhs.ifa-5A; winter wheat families)', cls='MTA', disease='Fusarium head blight', chr='5A', source='Herter2019TAG', peak=p, grade='C', marker='wsnp_BE444720A_TA_2_1', pve='22.8',
        note=HE + f'In data set 2 wsnp_BE444720A_TA_2_1, in LD with Fhb5 crossed in from the exotic parents, explained 22.8% of the genotypic variance for FHB severity; Fhb1, also carried by the parents, was not detected. Placed at the marker ({p:.3f} Mb, {h}).',
        evidence='Results: DS2 marker-trait association on 5A'))
    N.append(mk('S134-004', name='FHB association on 4B (winter wheat families)', cls='MTA', disease='Fusarium head blight', chr='4B', source='Herter2019TAG', grade='D', marker='see Supplementary Figure 4',
        note=HE + 'A minor FHB association on 4B in data set 1 (not linked to Rht-B1, fixed for Rht-B1a); the marker is given only in the supplement, so the locus is unplaced.', evidence='Results; Supplementary Figure 4'))
    HC = '102 elite Chinese bread wheat lines from the CIMMYT gene bank (FHB nurseries at El Batan, Mexico, 2009-2012) haplotyped with markers for 10 validated FHB QTL. '
    r = mk('S135-001', name='Wuhan 1 2D QTL (Chinese germplasm at CIMMYT)', disease='Fusarium head blight', chr='2D', source='He2014Euph', grade='C', marker='Xwmc144, Xwmc245',
           note=HC + '57% of lines carried the 2DL QTL alleles of Wuhan 1 or CJ 9306; Wuhan 1 QTL markers Xwmc144 and Xwmc245.')
    r['note'] = r['note'][:-len(T)] + two(r, 'Xwmc144', 'Xwmc245', '2D') + T; r['evidence'] = 'Haplotyping table: WU_2D'; N.append(r)
    r = mk('S135-002', name='CJ 9306 2DL QTL (Chinese germplasm at CIMMYT)', disease='Fusarium head blight', chr='2D', source='He2014Euph', grade='C', marker='Xgwm157, Xgwm539',
           note=HC + '57% of lines carried the 2DL QTL alleles of Wuhan 1 or CJ 9306; CJ 9306 QTL markers Xgwm157 and Xgwm539.')
    r['note'] = r['note'][:-len(T)] + two(r, 'Xgwm157', 'Xgwm539', '2D') + T; r['evidence'] = 'Haplotyping table: CJ_2D'; N.append(r)
    N.append(mk('S135-003', name='Fhb1 (Chinese germplasm at CIMMYT)', disease='Fusarium head blight', chr='3B', source='He2014Euph', peak=FHB1, grade='C', marker='UMN10',
        note=HC + '26.5% of lines had the Sumai 3 allele at the 3BS (Fhb1) marker UMN10. Placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Haplotyping table: SU_3B'))
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
    json.dump(ops, open('/tmp/ops28.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], r['asm'].get('RefSeq v2.1'))
    print(len(ops), len(json.dumps(ops)))
