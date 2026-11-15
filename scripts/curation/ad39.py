"""ad39 (Nov 2026). Not added: Walter et al. 2008 (DON gene expression), Yang et al. 2006 (SSR diversity survey), Tessmann and Van Sanford 2019 Euphytica (morphology GWAS), Soltanloo et al. 2011 (diallel), Sucur et al. 2026 (Stb gene survey); one duplicate upload (Zhang et al. 2016)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Yang2010EJPP': ('Yang X, Ma J, Li H, Ma H, Yao J, Liu C (2010) Different genes can be responsible for crown rot resistance at different developmental stages of wheat and barley. Eur J Plant Pathol 128:495-502', 'https://doi.org/10.1007/s10658-010-9680-3', 'Yang et al. (2010)'),
 'Wilde2007MB': ('Wilde F, Korzun V, Ebmeyer E, Geiger HH, Miedaner T (2007) Comparison of phenotypic and marker-based selection for Fusarium head blight resistance and DON content in spring wheat. Mol Breed 19:357-370', 'https://doi.org/10.1007/s11032-006-9067-5', 'Wilde et al. (2007)'),
 'Timonova2013MB': ('Timonova EM, Leonova IN, Roder MS, Salina EA (2013) Marker-assisted development and characterization of a set of Triticum aestivum lines carrying different introgressions from the T. timopheevii genome. Mol Breed 31:123-136', 'https://doi.org/10.1007/s11032-012-9776-x', 'Timonova et al. (2013)'),
 'Xiao2016MB': ('Xiao X, Ohm HW, Hunt GJ, Poland JA, Kong L, Nemacheck JA, Williams CE (2016) Genotyping-by-sequencing to remap QTL for type II Fusarium head blight and leaf rust resistance in a wheat-tall wheatgrass introgression recombinant inbred population. Mol Breed 36:51', 'https://doi.org/10.1007/s11032-016-0472-0', 'Xiao et al. (2016)'),
 'ZhangX2016MB': ('Zhang X, Rouse MN, Nava IC, Jin Y, Anderson JA (2016) Development and verification of wheat germplasm containing both Sr2 and Fhb1. Mol Breed 36:85', 'https://doi.org/10.1007/s11032-016-0502-y', 'Zhang et al. (2016)'),
 'Zheng2017MB': ('Zheng Z, Gao S, Zhou M, Yan G, Liu C (2017) Enhancing Fusarium crown rot resistance by pyramiding large-effect QTL in common wheat (Triticum aestivum L.). Mol Breed 37:107', 'https://doi.org/10.1007/s11032-017-0708-7', 'Zheng et al. (2017)'),
 'WangH2025NG': ('Wang H, Xu B, Wang X, Cheng R, Yuan Y, Li N, Shi J, Kong Z, Zhong J, Zhang X, Luo X, Xie Q, Li G, Jia H, Ma Z (2025) Mutagenesis and analysis of contrasting wheat lines do not support a role for PFT in Fusarium head blight resistance. Nat Genet 57:310-313', 'https://doi.org/10.1038/s41588-024-02056-3', 'Wang et al. (2025)'),
 'Yu2008MolBiol': ('Yu GH, Ma HX, Bai GH, Tang KX (2008) Single-strand conformational polymorphism markers associated with a major QTL for Fusarium head blight resistance in wheat. Mol Biol 42:504-513', 'https://doi.org/10.1134/S0026893308040043', 'Yu et al. (2008)'),
 'ZhangZ2021PJ': ('Zhang Z, Running KLD, Seneviratne S, Peters Haugrud AR, Szabo-Hever A, Shi G, Brueggeman R, Xu SS, Friesen TL, Faris JD (2021) A protein kinase-major sperm protein gene hijacked by a necrotrophic fungal pathogen triggers disease susceptibility in wheat. Plant J 106:720-732', 'https://doi.org/10.1111/tpj.15194', 'Zhang et al. (2021)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad39.py).'
FHB1, FHB5, QCRS3B, SNN3D1 = 8.528, (105.4338, 214.2), (820.501, 822.202), 4.343
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
    keep = [r for r in reports if not re.match(r'S1(8[6-9]|9[0-4])-|S184-002', r['id'])]; N = []
    r = mk('S186-001', name='Qcrs.cpi-4B (Lang; soil-less seedling assay)', disease='Fusarium crown rot', chr='4B', source='Yang2010EJPP', marker='Xbarc199, wPt-7569, wPt-4918', pve='18.9', stat='LOD 4.0',
           note='A soil-less seedling assay for crown rot (inoculated seedlings in moist paper towel) applied to the Lang x CSCR6 RIL population detected only the 4B QTL from Lang (trial 1 Barc199-wPt-7569, LOD 3.0, 14.7%; combined LOD 4.0, 18.9%); the major 3BL QTL (Qcrs.cpi-3B) found with the soil-based assay was not detected, suggesting different genes act at different stages or assays.', evidence='Table 2')
    place(r, ['Xbarc199', 'wPt-7569'], '4B'); N.append(r)
    WI = 'Donor QTL alleles from CM82036 (Sumai 3/Thornbird; 3B, 5A) and Frontana (3A) combined in German spring wheat by a double cross; two- and three-step phenotypic selection compared with marker-based selection (one to three flanking markers per QTL), tested with F. culmorum at four locations in 2004. All selection variants reduced FHB severity and DON; three-step phenotypic selection gave the highest total gain, marker selection for the two CM82036 alleles the highest gain per year. '
    N.append(mk('S187-001', name='Fhb1 (CM82036; MAS vs phenotypic selection)', disease='Fusarium head blight', chr='3B', source='Wilde2007MB', peak=FHB1, grade='C', marker='WMS389, WMS533, BARC133', note=WI + 'Placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Abstract; Methods'))
    N.append(mk('S187-002', name='Qfhs.ifa-5A / Fhb5 (CM82036; MAS vs phenotypic selection)', disease='Fusarium head blight', chr='5A', source='Wilde2007MB', start=FHB5[0], end=FHB5[1], grade='C', marker='WMS156, WMS304', note=WI + 'Placed at the Fhb5 interval in the atlas (105.4-214.2 Mb).', evidence='Abstract; Methods'))
    N.append(mk('S187-003', name='Frontana 3A FHB QTL (MAS vs phenotypic selection)', disease='Fusarium head blight', chr='3A', source='Wilde2007MB', grade='D', marker='WMS720', note=WI + 'The Frontana 3A QTL (Steiner et al. 2004) was tracked with WMS720, which could not be placed.', evidence='Abstract; Methods'))
    TI = 'Common wheat introgression lines carrying T. timopheevii segments, developed by marker-assisted backcrossing of two BC1F20 T. aestivum*2/T. timopheevii lines into Saratovskaya 29 (SSR and in situ hybridisation characterisation), tested with the West Siberian P. triticina population. '
    for n, (q, ch, nt) in enumerate([('QLr.icg-5B', '5B', 'Major leaf rust resistance locus on the introgressed 5GL segment; lines with a single 5GL introgression were completely resistant.'),
                                     ('QLr.icg-2A', '2A', 'Minor locus on a T. timopheevii At-genome segment; suppressed disease development and reduced urediniospore number by up to 25%, without a hypersensitive response.'),
                                     ('QLr.icg-1A', '1A', 'Minor locus on a T. timopheevii At-genome segment; suppressed disease development and reduced urediniospore number by up to 25%, without a hypersensitive response.')], 1):
        N.append(mk(f'S188-{n:03d}', name=q + ' (T. timopheevii introgression)', disease='Leaf rust', chr=ch, source='Timonova2013MB', grade='D', note=TI + nt + ' Originally mapped by Leonova et al. (2008, 2011); no marker positions in this paper, so not placed.', evidence='Abstract; Results'))
    XI = 'GBS (1,700 SNPs) linkage and QTL mapping in a wheat-tall wheatgrass (Th. ponticum) introgression RIL population; the parents proved to be translocation lines in two wheat backgrounds rather than 7E-7D substitution lines in Thatcher. '
    N.append(mk('S189-001', name='Qfhs.pur-7EL (GBS remapping; 7el2 on wheat 7DL)', disease='Fusarium head blight', chr='7E', source='Xiao2016MB', grade='D', ref='alien (Th. ponticum 7el2L, translocated to 7DL)',
        note=XI + 'Type II FHB QTL Qfhs.pur-7EL mapped to a 2.9-cM region (four markers) within a 43.6-cM segment of 7el2 translocated onto 7DL; two linked SNPs converted to KASP markers; five lines carry the resistance alleles of Qfhs.pur-7EL and Lr19 in coupling. Alien chromatin: no wheat RefSeq position.', evidence='Abstract'))
    N.append(mk('S189-002', name='Lr19 (GBS remapping; 7el1 on wheat 7DL)', cls='DG', disease='Leaf rust', chr='7E', source='Xiao2016MB', grade='D', ref='alien (Th. ponticum 7el1L, translocated to 7DL)',
        note=XI + 'Lr19 from 7el1 mapped to a 1.21-cM region containing two markers, in repulsion with Qfhs.pur-7EL. Alien chromatin: no wheat RefSeq position.', evidence='Abstract'))
    N.append(mk('S184-002', name='Lr19 (L. ponticum 7EL)', cls='DG', disease='Leaf rust', chr='7E', source='Zhang2011TAG', grade='D', ref='alien (L. ponticum 7EL)', marker='Xwmc273, XBE404744',
        note='In the 7E map from 237 F7:8 RILs of K11463 (7el1(7D)) x K2620 (7el2(7D)), Lr19 was bracketed by Xwmc273 and XBE404744, 1.54 and 1.43 cM away. Alien chromosome: no wheat RefSeq position. Read from the paper, Nov 2026 (ad39.py; added to the ad38 record set).'.replace(' Read from the paper, Nov 2026 (ad39.py; added to the ad38 record set).', ''), evidence='Abstract'))
    ZX = 'Of 1,949 F2 plants screened with Xgwm533 and csSr2 (Sr2) and UMN10/UMN10v2 (Fhb1), 43 had Fhb1 and Sr2 in coupling; 59 homozygous F3 spring lines carrying both had significantly lower FHB severity and moderate to high stem rust resistance in field nurseries. '
    N.append(mk('S190-001', name='Fhb1 (combined with Sr2)', disease='Fusarium head blight', chr='3B', source='ZhangX2016MB', peak=FHB1, grade='C', marker='UMN10, UMN10v2', note=ZX + 'Placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Abstract'))
    N.append(mk('S190-002', name='Sr2 (combined with Fhb1)', cls='DG', disease='Stem rust', chr='3B', source='ZhangX2016MB', ref='RefSeq v2.1', peak=1.17, grade='C', marker='csSr2, Xgwm533', note=ZX + 'Placed at the Sr2 position in the atlas (1.17 Mb, v2.1).', evidence='Abstract'))
    ZH = 'Pyramiding of three large-effect crown rot QTL in two populations (EGA Wylie//Lang/CSCR6, 178 lines; EGA Wylie/Sumai3 x CSCR6 derived, 181 F8 lines), seedling assays in four trials; lines with two resistant alleles beat those with one or none, and lines with all three were most resistant. '
    N.append(mk('S191-001', name='Qcrs.cpi-3B (CSCR6; pyramiding)', disease='Fusarium crown rot', chr='3B', source='Zheng2017MB', start=QCRS3B[0], end=QCRS3B[1], grade='C', marker='Xgwm181', note=ZH + 'Resistant allele from T. spelta CSCR6 reduced disease index by 25.5-28.5%. Placed at the Qcrs.cpi-3B interval in the atlas.', evidence='Table 1'))
    r = mk('S191-002', name='Qcrs.cpi-5D (EGA Wylie; pyramiding)', disease='Fusarium crown rot', chr='5D', source='Zheng2017MB', marker='Xcfd189 (Xcfd89 not polymorphic)', note=ZH + 'Resistant allele from EGA Wylie (5DS) reduced disease index by 16.6-22.9%.', evidence='Table 1'); place(r, ['Xcfd189'], '5D'); N.append(r)
    r = mk('S191-003', name='Qcrs.cpi-2D (EGA Wylie; pyramiding)', disease='Fusarium crown rot', chr='2D', source='Zheng2017MB', marker='Xcfd73, Xgwm539', note=ZH + 'Resistant allele from EGA Wylie (2DL) reduced disease index by 13.4-19.3%.', evidence='Table 1'); place(r, ['Xcfd73', 'Xgwm539'], '2D'); N.append(r)
    N.append(mk('S192-001', name='Fhb1 = HisR/TaHRC (PFT not required)', cls='CG', disease='Fusarium head blight', chr='3B', source='WangH2025NG', peak=FHB1, grade='C', marker='HisR, PFT',
        note='Lines carrying HisR (TaHRC) with or without PFT, and five HisR-carrying PFT mutants made by gene editing, showed that PFT does not confer FHB resistance and that the HisR effect does not require PFT; PFT mutations were not associated with morphological or phenological traits. This argues against PFT as the Fhb1 gene (Rawat et al. 2016) and supports HisR/TaHRC, about 50 kb from PFT. Placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Abstract'))
    r = mk('S193-001', name='Fhb1 region SSCP markers (Ning 7840)', disease='Fusarium head blight', chr='3B', source='Yu2008MolBiol', marker='Xsscp6, Xsscp20, Xsscp21, Xgwm493, Xgwm533', pve='up to 50',
           note='EST-derived SSCP markers for the 3BS QTL (between Xgwm389 and Xgwm493) from Sumai 3 derivative Ning 7840 in a Ning 7840 x Clark population; Xsscp6, Xsscp20 and Xsscp21 mapped between Xgwm493 and Xgwm533 and explained more variation (R2 50, 47 and 46% in F6:7) than Xgwm533 (45%).', evidence='Abstract; Table 2')
    place(r, ['Xgwm493', 'Xgwm533'], '3B'); N.append(r)
    N.append(mk('S194-001', name='Snn3-D1 (cloned; protein kinase-major sperm protein)', cls='SG', disease='Septoria nodorum blotch', chr='5D', source='ZhangZ2021PJ', peak=SNN3D1, grade='C', marker='',
        note='Map-based cloning of Snn3-D1 from Ae. tauschii (TA2377 x AL8/78 population, BAC contig on AL8/78 scaffold 5289): the gene encodes protein kinase and major sperm protein domains, both essential for SnTox3 sensitivity; the interaction is light-independent and expression is down-regulated by light. Snn3-D1 likely arose in Ae. tauschii through a ~218-kb insertion along the west bank of the Caspian Sea; its homoeolog TaPKMSP-5B lies on 5B. Placed at the Snn3-D1 position in the atlas (Lin et al. 2021).', evidence='Summary; Results'))
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
    json.dump(ops, open('/tmp/ops39.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], (r['asm'].get('RefSeq v2.1') or {}).get('flag'))
    print(len(ops), len(json.dumps(ops)))
