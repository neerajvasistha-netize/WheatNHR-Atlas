"""ad37 (Nov 2026). Not added: Song et al. preprint (TaNDC1) and Song et al. 2025 IJMS (TaU11/U12-35K), functional gene studies; Stewart et al. 2016 (phenotyping method); Steed et al. 2002 (conference version of Steed et al. 2005, in atlas); Shi et al. 2008 (haplotype survey). Erratum to Szabo-Hever et al. 2014: author names corrected in V_Agnes2014 citation."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Simon2004Euph': ('Simon MR, Ayala FM, Cordo CA, Roder MS, Borner A (2004) Molecular mapping of quantitative trait loci determining resistance to septoria tritici blotch caused by Mycosphaerella graminicola in wheat. Euphytica 138:41-48', 'https://doi.org/10.1023/B:EUPH.0000047059.57839.98', 'Simon et al. (2004)'),
 'Simon2007Euph': ('Simon MR, Ayala FM, Cordo CA, Roder MS, Borner A (2007) The use of wheat/goatgrass introgression lines for the detection of gene(s) determining resistance to septoria tritici blotch (Mycosphaerella graminicola). Euphytica 154:249-254', 'https://doi.org/10.1007/s10681-006-9293-2', 'Simon et al. (2007)'),
 'Somers2005TAG': ('Somers DJ, Thomas J, DePauw R, Fox S, Humphreys G, Fedak G (2005) Assembling complex genotypes to resist Fusarium in wheat (Triticum aestivum L.). Theor Appl Genet 111:1623-1631', 'https://doi.org/10.1007/s00122-005-0094-0', 'Somers et al. (2005)'),
 'Srinivasachary2008TAG': ('Srinivasachary, Gosman N, Steed A, Simmonds J, Leverington-Waite M, Wang Y, Snape J, Nicholson P (2008) Susceptibility to Fusarium head blight is associated with the Rht-D1b semi-dwarfing allele in wheat. Theor Appl Genet 116:1145-1153', 'https://doi.org/10.1007/s00122-008-0742-2', 'Srinivasachary et al. (2008)'),
 'Sheng2012TAG': ('Sheng H, See DR, Murray TD (2012) Mapping QTL for resistance to eyespot of wheat in Aegilops longissima. Theor Appl Genet 125:355-366', 'https://doi.org/10.1007/s00122-012-1838-2', 'Sheng et al. (2012)'),
 'Sheng2014TAG': ('Sheng H, See DR, Murray TD (2014) Mapping resistance genes for Oculimacula acuformis in Aegilops longissima. Theor Appl Genet 127:2085-2093', 'https://doi.org/10.1007/s00122-014-2361-4', 'Sheng et al. (2014)'),
 'Su2018TAG': ('Su Z, Jin S, Zhang D, Bai G (2018) Development and validation of diagnostic markers for Fhb1 region, a major QTL for Fusarium head blight resistance in wheat. Theor Appl Genet 131:2371-2380', 'https://doi.org/10.1007/s00122-018-3159-6', 'Su et al. (2018)'),
 'Slominska2020TAG': ('Slominska-Durdasiak KM, Kollers S, Korzun V, Nagel M, Schweizer P, Djamei A, Reif JC (2020) Association mapping of wheat Fusarium head blight resistance-related regions using a candidate-gene approach and their verification in a biparental population. Theor Appl Genet 133:341-351', 'https://doi.org/10.1007/s00122-019-03463-5', 'Slominska-Durdasiak et al. (2020)'),
 'Tadesse2008Euph': ('Tadesse W, Hsam SLK, Wenzel G, Zeller FJ (2008) Chromosome location of a gene conferring resistance to Pyrenophora tritici-repentis in Ethiopian wheat cultivars. Euphytica 162:423-430', 'https://doi.org/10.1007/s10681-007-9593-1', 'Tadesse et al. (2008)'),
 'ShenOhm2007MB': ('Shen X, Ohm H (2007) Molecular mapping of Thinopyrum-derived Fusarium head blight resistance in common wheat. Mol Breed 20:131-140', 'https://doi.org/10.1007/s11032-007-9079-9', 'Shen and Ohm (2007)'),
 'Singh2010MB': ('Singh PK, Mergoum M, Adhikari TB, Shah T, Ghavami F, Kianian SF (2010) Genetic and molecular analysis of wheat tan spot resistance effective against Pyrenophora tritici-repentis races 2 and 5. Mol Breed 25:369-379', 'https://doi.org/10.1007/s11032-009-9336-1', 'Singh et al. (2010)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad37.py).'
FHB1, RHTD1, TSN1, TSC2 = 8.528, 18.781062, (546.698, 546.843), (24.09, 25.02)
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
    keep = [r for r in reports if not re.match(r'S1(6[4-9]|7[0-5])-', r['id'])]; N = []
    r = mk('S164-001', name='Qfhs.ndsu-3AS (T. dicoccoides; validated in Argentine durum)', disease='Fusarium head blight', chr='3A', source='V_Soresi2017', marker='Xgwm2',
           note='Langdon(Dic-3A)-10 carrying Qfhs.ndsu-3AS (type II resistance from T. turgidum ssp. dicoccoides) was crossed with Argentine durum cultivars; F4 progeny screened with the tightly linked Xgwm2 and point-inoculated with F. graminearum. Carriers (heterozygous or homozygous) had significantly lower severity at 21 dpi and slower disease progress than the cultivated parents, showing stable dominant expression in durum backgrounds.', evidence='Abstract; Results')
    place(r, ['Xgwm2'], '3A'); N.append(r)
    SI = 'STB resistance in 65 ITMI RILs (W7984 x Opata 85) inoculated with Argentine isolates IPO 92067 and IPO 93014 at seedling and adult stages; RFLP-based ITMI map. '
    for n, (nm, ch, ms, nt) in enumerate([
        ('QStb.ipk-1D', '1D', 'Xbcd1434, Xksud14d', 'Seedling QTL on distal 1DS, both isolates (LOD 3.0 or more with IPO 92067), resistance from W7984.'),
        ('QStb.ipk-2D', '2D', 'Xbcd260a, Xbcd111', 'Seedling QTL on 2DS near the centromere, both isolates (LOD 1.5-3.0), resistance from Opata 85.'),
        ('QStb.ipk-6B', '6B', 'Xksuh14b, Xbcd1383', 'Seedling QTL on 6BS, both isolates (LOD 3.0 or more with IPO 92067), resistance from W7984.'),
        ('QStb.ipk-3D', '3D', 'Xbcd22a, Xbcd1555a', 'Adult-plant QTL on 3DL specific to IPO 92067 (early milk, dough stage and AUDPC), resistance from W7984.'),
        ('QStb.ipk-7B', '7B', 'Xwg686, Xksud2a', 'Adult-plant QTL on 7BL specific to IPO 93014 (flowering, early milk, dough stage and AUDPC), resistance from W7984.')], 1):
        N.append(mk(f'S165-{n:03d}', name=nm, disease='Septoria tritici blotch', chr=ch, source='Simon2004Euph', grade='D', marker=ms, pve='18-25 (seedling, IPO 93014); 10-16 (seedling, IPO 92067); 16-21 (adult)',
            note=SI + nt + ' Markers read from Fig. 2; RFLP loci, not placed.', evidence='Results; Fig. 2'))
    r = mk('S166-001', name='Stb5 (Ae. tauschii 7D introgression lines)', cls='DG', disease='Septoria tritici blotch', chr='7D', source='Simon2007Euph', marker='Xgwm44',
           note='Thirteen Chinese Spring/Ae. tauschii chromosome 7D introgression lines, CS and CS(Syn 7D) tested with Argentine isolates IPO 92067 and IPO 93014 at seedling and adult stages (two locations, 2003-2004); resistance was effective against both isolates at both stages and mapped to the centromeric region of 7DS; from its relationship with Xgwm44 the gene is likely Stb5.', evidence='Abstract; Fig. 2')
    place(r, ['Xgwm44'], '7D'); N.append(r)
    SO = 'Molecular breeding (high-throughput SSR genotyping and half-seed analysis) of four crosses to introduce six FHB resistance QTL from Sumai 3, Wuhan and Nyubai into Canadian spring wheat; QTL tracked with three markers each (Table 2). '
    for n, (nm, ch, ms, arm, nt) in enumerate([
        ('3BS (Fhb1)', '3B', ['Xgwm533', 'Xgwm493', 'Xwmc808'], None, 'type II, distal 3BS, 33-cM interval'),
        ('3BSc', '3B', ['Xwmc625', 'Xgwm566', 'Xwmc418'], 'S', 'type I, near the 3B centromere, 13-cM interval'),
        ('4B', '4B', ['Xwmc710', 'Xwmc238', 'Xgwm149'], None, 'type II, 20-cM interval'),
        ('2DL', '2D', ['Xwmc245', 'Xgwm608', 'Xgwm301'], 'L', 'type I, 40-cM interval'),
        ('5AS', '5A', ['Xgwm293', 'Xwmc705', 'Xgwm154'], 'S', 'types I and II, 18-cM interval'),
        ('6B', '6B', ['Xwmc104', 'Xwmc397', 'Xgwm219'], 'S', 'types I and II, 45-cM interval (6BS)')], 1):
        r = mk(f'S167-{n:03d}', name=f'FHB QTL {nm} (Canadian molecular breeding)', disease='Fusarium head blight', chr=ch, source='Somers2005TAG', marker=', '.join(ms), note=SO + f'QTL {nm}: {nt}.', evidence='Table 2')
        place(r, ms, ch, arm); N.append(r)
    SR = 'FHB resistance (AUDPC, spray inoculation) in the Spark x Rialto DH population in three field trials (A2005, N2005, J2006); MQM mapping with DArT and SSR markers (LOD threshold 2.8). '
    for n, (nm, ch, ms, org, det, pv) in enumerate([
        ('Qfhs.jic-1b', '1B', ['wPt-0705'], 'Rialto', 'N2005 LOD 5.2 (13.4%)', '13.4'),
        ('Qfhs.jic-2a', '2A', ['Xgwm515'], 'Spark', 'N2005 LOD 2.9 (7%)', '7'),
        ('Qfhs.jic-3a.1', '3A', ['Xwmc11', 'wPt-7992'], 'Rialto', 'J2006 LOD 2.7 (2.7%, below threshold); across trials LOD 2.8 (4%)', '4'),
        ('Qfhs.jic-3a.2', '3A', ['Xgwm497', 'Xbarc19'], 'Spark', 'N2005 LOD 3.3 (8%); A2005 LOD 2.8 (7.6%)', '8'),
        ('Qfhs.jic-4d.2', '4D', ['Xgwm265'], 'Spark', 'J2006 LOD 3.1 (3.1%); at least 24 cM from Qfhs.jic-4d.1', '3.1'),
        ('Qfhs.jic-5a', '5A', ['Xgwm443'], 'Spark', 'J2006 LOD 3.5 (3.9%)', '3.9'),
        ('Qfhs.jic-6a.2', '6A', ['wPt-8833'], 'Spark', 'J2006 LOD 8.3 (10.5%); across trials LOD 5.2 (7.8%)', '10.5'),
        ('Qfhs.jic-7a', '7A', ['Xpsp3050', 'wPt-7299'], 'Spark', 'J2006 LOD 2.9 (3.7%); N2005 LOD 2.2 (5.2%, below threshold)', '5.2')], 1):
        r = mk(f'S168-{n:03d}', name=nm + ' (Spark x Rialto)', disease='Fusarium head blight', chr=ch, source='Srinivasachary2008TAG', marker=', '.join(ms), pve=pv,
               note=SR + f'Resistance from {org}. {det}.', evidence='Table 3; Fig. 1')
        place(r, ms, ch); N.append(r)
    N.append(mk('S168-009', name='Qfhs.jic-4d.1 at Rht-D1 (Spark x Rialto)', disease='Fusarium head blight', chr='4D', source='Srinivasachary2008TAG', peak=RHTD1, grade='C', marker='XRht-D1b, Xpsp3103', pve='up to 52.7', stat='LOD 6.8-28.2',
        note=SR + 'Major QTL coincident with Rht-D1: resistance from Spark (Rht-D1a), susceptibility associated with Rialto Rht-D1b; up to 52.7% (J2006, LOD 28.2), 50.9% across trials. Not an effect of height per se; near-isogenic lines showed Rht-D1b compromises type I resistance but not type II. Placed at the Rht-D1 gene (TraesCS4D02G040400, gene table).', evidence='Table 3; Abstract'))
    SH = 'Eyespot (O. yallundae) resistance in 189 F5 RILs of Ae. longissima PI 542196 (R) x PI 330486 (S), Sl-genome map of 169 wheat SSRs; four QTL explained 44% (GUS score) and 63% (visual rating) of the variation, all from PI 542196. '
    for n, (q, ch, mk_, nt) in enumerate([('Q.Pch.wsu-1Sl', '1S', 'cfd6, gdm67, gwm642, cfd48', 'Closest markers cfd6, gdm67, gwm642 and cfd48 (3.1-cM cluster).'),
                                         ('Q.Pch.wsu-3Sl', '3S', 'wmc169, gdm72, wmc597, wmc231', 'In a 39.2-cM interval between wmc169 and wmc231, probably on the short arm; 14% (GUS) and 9% (visual), LOD 5.3 and 3.4.'),
                                         ('Q.Pch.wsu-5Sl', '5S', 'wmc415', 'Closest marker wmc415.'),
                                         ('Q.Pch.wsu-7Sl', '7S', 'cfd2', 'Closest marker cfd2.')], 1):
        N.append(mk(f'S169-{n:03d}', name=q + ' (Ae. longissima)', disease='Eyespot', chr=ch, source='Sheng2012TAG', grade='D', ref=f'alien (Ae. longissima {ch}l)', marker=mk_,
            note=SH + nt + ' Ae. longissima chromosome: no wheat RefSeq position.', evidence='Abstract; Results'))
    SH2 = 'Resistance to O. acuformis (eyespot) in the same Ae. longissima RIL population (PI 542196 x PI 330486); three QTL explained 66% (GUS) and 84% (visual) of the variation, all from PI 542196. '
    for n, (q, ch, mk_, nt) in enumerate([('Q.Pch-oa.wsu-1Sl', '1S', 'gdm67, gwm642, cfd48', 'Probably on 1SlL; closest markers within 8.7 cM.'),
                                         ('Q.Pch-oa.wsu-3Sl', '3S', 'cfd79, cfd9, barc71', 'LOD 8.8 (GUS) and 10.0 (visual), 27% and 30%; 40.6-cM interval cfd79-barc71, cfd9 about 6 cM from the peak; probably on 3SlS.'),
                                         ('Q.Pch-oa.wsu-5Sl', '5S', 'wmc415, cfd12, gwm408, gwm271, barc319', 'LOD 6.0 (GUS) and 13.0 (visual), 12% and 23%; wmc415-barc319 (24.9 cM), cfd12, gwm408 and gwm271 closest; probably on 5SlL.')], 1):
        N.append(mk(f'S170-{n:03d}', name=q + ' (Ae. longissima)', disease='Eyespot', chr=ch, source='Sheng2014TAG', grade='D', ref=f'alien (Ae. longissima {ch}l)', marker=mk_,
            note=SH2 + nt + ' Ae. longissima chromosome: no wheat RefSeq position.', evidence='Abstract; Results'))
    N.append(mk('S171-001', name='Fhb1 (TaHRC diagnostic markers)', cls='CG', disease='Fusarium head blight', chr='3B', source='Su2018TAG', peak=FHB1, grade='C', marker='TaHRC deletion markers',
        note='Two functional markers based on the deletion in TaHRC (putative histidine-rich calcium-binding protein) in the Fhb1 region, developed by comparing near-isogenic lines and validated in populations and three diversity panels; diagnostic for Fhb1 and more accurate than earlier markers. Placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Abstract'))
    SL = 'Candidate-gene association mapping (14 KASP markers in nine cross-kingdom RNA trafficking genes) in 356 GABI winter wheat lines, verified in the Apache x Biscay RIL population. '
    N.append(mk('S172-001', name='Dicer1 (DCL1) FHB association', cls='CG', disease='Fusarium head blight', chr='5A', source='Slominska2020TAG', grade='D', marker='DCL1(1), DCL1(2)',
        note=SL + 'Two markers in Dicer1 on 5A were associated with FHB severity (P 1.5e-4 and 9.0e-6) and verified in Apache x Biscay, without association with height or heading. No gene ID given; not placed.', evidence='Table 1; Table 2'))
    N.append(mk('S172-002', name='Ara6 (Rab5-like GTPase) FHB association', cls='CG', disease='Fusarium head blight', chr='4D', source='Slominska2020TAG', grade='D', marker='Ara6',
        note=SL + 'An intronic marker in Ara6 on 4D was associated with FHB severity (P 2.2e-4) and verified in Apache x Biscay, where it was also associated with plant height. No gene ID given; not placed.', evidence='Table 1; Table 2'))
    N.append(mk('S173-001', name='TsrHar (Ethiopian cultivars HAR604, HAR2562)', cls='DG', disease='Tan spot', chr='3B', source='Tadesse2008Euph', grade='D', marker='',
        note='Of 118 Ethiopian bread and durum wheats tested with P. tritici-repentis isolates ASC1a and DW-16, HAR604, HAR2562 and Dashen were resistant to both and carry allelic recessive resistance genes; monosomic analysis placed the gene, temporarily designated TsrHar, on 3B. Its relation to Tsr2 is unknown. Chromosome-level location only.', evidence='Abstract'))
    N.append(mk('S174-001', name='Qfhs.pur-7EL (Th. ponticum 7el2)', disease='Fusarium head blight', chr='7E', source='ShenOhm2007MB', grade='D', ref='alien (Th. ponticum 7el2L)', marker='XBE445653, Xcfa2240', pve='15.1-32.5',
        note='Type II FHB resistance mapped in a population from two wheat-Th. ponticum 7el substitution lines (el1, susceptible; el2, resistant): Qfhs.pur-7EL in the distal long arm of 7el2 between XBE445653 and Xcfa2240 reduced diseased spikelets by 1.5-2.6 and explained 15.1-32.5% in four experiments; 7el2-specific STS/EST markers developed. Alien chromosome: no wheat RefSeq position.', evidence='Abstract'))
    SG = 'Tan spot resistance in F1, F2 and 118 F2:6 RILs of ND 735 (R) x Steele-ND (S), seedlings tested with spore inocula and culture filtrates of P. tritici-repentis races 2 and 5. '
    N.append(mk('S175-001', name='Tsr1 (ND 735; race 2 necrosis insensitivity)', cls='DG', disease='Tan spot', chr='5B', source='Singh2010MB', start=TSN1[0], end=TSN1[1], grade='C', marker='wPt-8285, wPt-3049',
        note=SG + 'A single recessive gene, Tsr1, on 5BL controlled resistance and insensitivity to race 2 necrosis, flanked by wPt-8285 (7.0 cM) and wPt-3049 (2.9 cM). The DArT markers could not be placed; placed at the Tsn1 position in the atlas (Tsr1 is the tsn1 resistance allele).', evidence='Abstract; Fig. 2'))
    N.append(mk('S175-002', name='Tsr6 (ND 735; race 5 chlorosis insensitivity)', cls='DG', disease='Tan spot', chr='2B', source='Singh2010MB', start=TSC2[0], end=TSC2[1], grade='C', marker='wPt-0289, Xwmc25',
        note=SG + 'A second recessive gene, designated Tsr6, on 2BS conferred resistance and insensitivity to chlorosis induced by race 5 spores or culture filtrates, flanked by wPt-0289 (4.6 cM) and Xwmc25 (18.7 cM). The markers could not be placed; placed at the Tsc2 (Ptr ToxB sensitivity, 2BS) position in the atlas as the likely corresponding locus.', evidence='Abstract; Fig. 3'))
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
    json.dump(ops, open('/tmp/ops37.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], (r['asm'].get('RefSeq v2.1') or {}).get('flag'))
    print(len(ops), len(json.dumps(ops)))
