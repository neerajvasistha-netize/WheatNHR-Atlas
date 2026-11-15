"""ad36 (Nov 2026): Rahman 2021, Nishio 2008, Rawat 2016, Mohler 2002, Rosyara 2009, Salameh 2011, Miedaner 2011. Not added: Shen 2004 (substitution lines), Mirdita 2015 (genomic prediction, no loci), Running 2025 (already in atlas), Muhovski 2014 (expression and arm assignment only)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Rahman2021FCR': ('Rahman MM, Davies P, Bansal U, Pasam R, Hayden M, Trethowan R (2021) Relationship between resistance and tolerance of crown rot in bread wheat. Field Crops Res 265:108106', 'https://doi.org/10.1016/j.fcr.2021.108106', 'Rahman et al. (2021)'),
 'Nishio2008BS': ('Nishio Z, Ito M, Tanio M, Tabiki T, Yamauchi H (2008) Microsatellite marker-assisted selection of fusarium head blight resistance in backcrossed Japanese winter wheat lines. Breed Sci 58:183-185', 'https://doi.org/10.1270/jsbbs.58.183', 'Nishio et al. (2008)'),
 'Rawat2016NG': ('Rawat N, Pumphrey MO, Liu S, Zhang X, Tiwari VK, Ando K, Trick HN, Bockus WW, Akhunov E, Anderson JA, Gill BS (2016) Wheat Fhb1 encodes a chimeric lectin with agglutinin domains and a pore-forming toxin-like domain conferring resistance to Fusarium head blight. Nat Genet 48:1576-1580', 'https://doi.org/10.1038/ng.3706', 'Rawat et al. (2016)'),
 'Mohler2002TAG': ('Mohler V, Klahr A, Wenzel G, Schwarz G (2002) A resistance gene analog useful for targeting disease resistance genes against different pathogens on group 1S chromosomes of barley, wheat and rye. Theor Appl Genet 105:364-368', 'https://doi.org/10.1007/s00122-002-0919-z', 'Mohler et al. (2002)'),
 'Rosyara2009TAG': ('Rosyara UR, Gonzalez-Hernandez JL, Glover KD, Gedye KR, Stein JM (2009) Family-based mapping of quantitative trait loci in plant breeding populations with resistance to Fusarium head blight in wheat as an illustration. Theor Appl Genet 118:1617-1631', 'https://doi.org/10.1007/s00122-009-1010-9', 'Rosyara et al. (2009)'),
 'Salameh2011MB': ('Salameh A, Buerstmayr M, Steiner B, Neumayer A, Lemmens M, Buerstmayr H (2011) Effects of introgression of two QTL for fusarium head blight resistance from Asian spring wheat by marker-assisted backcrossing into European winter wheat on fusarium head blight resistance, yield and quality traits. Mol Breed 28:485-494', 'https://doi.org/10.1007/s11032-010-9498-x', 'Salameh et al. (2011)'),
 'Miedaner2011MB': ('Miedaner T, Wurschum T, Maurer HP, Korzun V, Ebmeyer E, Reif JC (2011) Association mapping for Fusarium head blight resistance in European soft winter wheat. Mol Breed 28:647-655', 'https://doi.org/10.1007/s11032-010-9516-z', 'Miedaner et al. (2011)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad36.py).'
FHB1, FHB5, RHTD1 = 8.528, (105.4338, 214.2), 18.781062
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
    keep = [r for r in reports if not re.match(r'S1(5[7-9]|6[0-3])-', r['id'])]; N = []
    RA = 'Crown rot (F. pseudograminearum) GWAS with the 90K SNP array in 196 F4/F5 progeny of the three-way cross AUS29529/2/Syn110/3/Sunco/2*Pastor (207 genotypes with parents and checks), paired inoculated and uninoculated field plots at Narrabri in 2014 and 2015; 43 QTL, of which three for resistance (white stem browning) and 40 for tolerance (yield traits under inoculation; not entered). '
    r = mk('S157-001', name='QTL_7 crown rot resistance 2B (Rahman 2021)', disease='Fusarium crown rot', chr='2B', source='Rahman2021FCR', marker='Tdurum_contig83066_276', pve='6.12', stat='-log10P 3.07',
           note=RA + 'QTL_7 reduced stem browning on adult plants and explained 6.12% of its variation (PopSeq/consensus position 292.15).', evidence='Table 6; Results 3.4'); place(r, ['Tdurum_contig83066_276'], '2B'); N.append(r)
    N.append(mk('S157-002', name='QTL_8 crown rot resistance 2B (Rahman 2021)', disease='Fusarium crown rot', chr='2B', source='Rahman2021FCR', grade='D', pve='6.2', stat='-log10P 3.17',
        note=RA + 'QTL_8 reduced stem browning (R2 6.2%, PopSeq/consensus position 317.99); its SNP is not named in the paper (Supplementary Table 2, not uploaded), so it is unplaced.', evidence='Table 6'))
    r = mk('S157-003', name='QTL_37 crown rot resistance 6A (Rahman 2021)', disease='Fusarium crown rot', chr='6A', source='Rahman2021FCR', marker='Ku_c21490_472, Ra_c8185_676', pve='5.96-6.49',
           note=RA + 'QTL_37, a haplotype of two SNPs at 188.89 (PopSeq/consensus), reduced stem browning and explained 5.96-6.49% of its variation.', evidence='Table 6; Results 3.4'); place(r, ['Ku_c21490_472', 'Ra_c8185_676'], '6A'); N.append(r)
    r = mk('S158-001', name='Fhb1 region (Sumai 3 into Kitamoe backcross lines)', disease='Fusarium head blight', chr='3B', source='Nishio2008BS', marker='Xbarc133, Xbarc147, Xgwm533, Xbarc102',
           note='Six 3BS SSRs tested on BC1F3 families of the FHB-susceptible Japanese winter wheat Kitamoe backcrossed with Sumai 3; families homozygous for the Sumai 3 alleles at barc133, barc147, gwm533.1 and barc102 had lower FHB severity than Kitamoe-type and heterozygous families (unlike barc87), suggesting recessive expression in the winter wheat background.', evidence='Table 1; Fig. 1')
    place(r, ['Xbarc133', 'Xbarc147', 'Xgwm533', 'Xbarc102'], '3B'); N.append(r)
    N.append(mk('S159-001', name='Fhb1 = PFT (pore-forming toxin-like gene, Sumai 3)', cls='CG', disease='Fusarium head blight', chr='3B', source='Rawat2016NG', peak=FHB1, grade='C', marker='STS3B-355, STS3B-334',
        note='Map-based cloning of Fhb1 (Qfhs.ndsu-3BS) from Sumai 3: the 0.08-cM interval between STS3B-355 and STS3B-334 was covered by a Sumai 3 BAC contig, and mutation analysis, gene silencing and transgenic overexpression showed that a pore-forming toxin-like (PFT) gene, encoding a chimeric lectin with two agglutinin domains and an ETX/MTX2 toxin domain, confers FHB resistance. The Chinese Spring Fhb1 region lacks an active copy. Placed at the Fhb1 position in the atlas (TaHRC, Su et al. 2019, which identified a different causal gene at the same locus).', evidence='Abstract; Fig. 1'))
    r = mk('S160-001', name='1DS FHB QTL at RGA XaACT/CAA (Cansas x Ritmo)', disease='Fusarium head blight', chr='1D', source='Mohler2002TAG', marker='XaACT/CAA, Xgwm106', pve='14.5', stat='LOD 3.2',
           note='The barley Mla-linked RGA marker aACT/CAA was assigned to 1DS with Chinese Spring aneuploids and mapped in 96 F5 RILs of Cansas x Ritmo 13.0 cM from Xgwm106; it explained 12% of FHB variance by single-marker regression, and composite interval mapping across four environments gave a LOD of 3.2 (14.5%) for this 1DS QTL. The RGA itself could not be placed.', evidence='Results; Fig. 4B')
    place(r, ['Xgwm106'], '1D'); N.append(r)
    N.append(mk('S161-001', name='Fhb1 (family-pedigree mapping in SDSU breeding populations)', disease='Fusarium head blight', chr='3B', source='Rosyara2009TAG', peak=FHB1, grade='C', marker='Xgwm533.1, Xbarc133, Xgwm493',
        note='Family-based QTL mapping (variance component, pedigree-wide regression and QTDT association) in 82 breeding families (793 individuals) re-detected Fhb1: Xgwm533.1-158, Xbarc133-118 and Xgwm493-192 alleles from resistant parents SD3776, Freyer and SD3851 were associated with FHB severity (QTDT P down to 4e-7 at 21 DAI); markers in the Xgwm389-Xgwm754 interval were linked. Placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Table 2; Results'))
    SA = 'Fhb1 and Qfhs.ifa-5A from CM-82036 transferred by marker-assisted backcrossing into nine European winter wheat lines (15 BC2-derived families); effect on FHB increased in the order Qfhs.ifa-5A < Fhb1 <= both combined, with a large effect of the recipient background and no systematic negative effect on grain yield, thousand grain weight, hectoliter weight or protein content. '
    N.append(mk('S162-001', name='Fhb1 (CM-82036; MABC into European winter wheat)', disease='Fusarium head blight', chr='3B', source='Salameh2011MB', peak=FHB1, grade='C', marker='Gwm533, Gwm493, Barc133',
        note=SA + 'Placed at the cloned Fhb1 (TaHRC) position in the atlas.', evidence='Abstract'))
    N.append(mk('S162-002', name='Qfhs.ifa-5A / Fhb5 (CM-82036; MABC into European winter wheat)', disease='Fusarium head blight', chr='5A', source='Salameh2011MB', start=FHB5[0], end=FHB5[1], grade='C', marker='Gwm293, Gwm156, Gwm1057',
        note=SA + 'Qfhs.ifa-5A alone had a relatively small effect. Placed at the Fhb5 interval in the atlas (105.4-214.2 Mb).', evidence='Abstract'))
    MI = 'Genome-wide association mapping of FHB rating in 455 European soft winter wheat lines (two environments, several ratings) with 115 SSRs and Rht-B1/Rht-D1 markers; nine significant associations on seven chromosomes, three remaining after Bonferroni-Holm correction (1B, 1D, 2D). '
    MM = [('gwm18', '1B', '33', '5.8', 'significant after Bonferroni-Holm'), ('cfd72', '1D', '51', '6.5', 'significant after Bonferroni-Holm; new region'), ('wmc18', '2D', '64', '11.7', 'significant after Bonferroni-Holm'),
          ('wmc336', '1D', '23', '5.0', 'P < 0.05'), ('barc101', '2B', '76', '3.8', 'P < 0.05'), ('wmc11', '3A', '0', '5.6', 'P < 0.05; new region'), ('wmc596', '7A', '73', '7.5', 'P < 0.05')]
    for n, (m, ch, cm, pg, sig) in enumerate(MM, 1):
        r = mk(f'S163-{n:03d}', name=f'FHB association {ch} (X{m}; European winter wheat)', disease='Fusarium head blight', chr=ch, source='Miedaner2011MB', marker='X' + m, pve=pg,
               note=MI + f'X{m} on {ch} ({cm} cM), {pg}% of genotypic variance, {sig}.', evidence='Table 2')
        place(r, ['X' + m], ch); N.append(r)
    N.append(mk('S163-008', name='FHB association 2B/3A (Xbarc45; European winter wheat)', disease='Fusarium head blight', chr='2B', source='Miedaner2011MB', grade='D', marker='Xbarc45', pve='3.4',
        note=MI + 'Xbarc45 maps to 2B or 3A (4/37 cM), 3.4% of genotypic variance, P < 0.05; chromosome ambiguous, so unplaced.', evidence='Table 2'))
    N.append(mk('S163-009', name='Rht-D1 FHB susceptibility (European winter wheat)', disease='Fusarium head blight', chr='4D', source='Miedaner2011MB', peak=RHTD1, grade='C', marker='sw4brht2 (Rht-D1)', pve='4.2',
        note=MI + 'The dwarfing locus Rht-D1 (marker sw4brht2) was associated with FHB rating (4.2% of genotypic variance, P < 0.05), the dwarfing allele increasing susceptibility; an epistatic interaction involved Rht-B1. Placed at the Rht-D1 gene (TraesCS4D02G040400, gene table).', evidence='Table 2; Table 3'))
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
    json.dump(ops, open('/tmp/ops36.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], (r['asm'].get('RefSeq v2.1') or {}).get('flag'))
    print(len(ops), len(json.dumps(ops)))
