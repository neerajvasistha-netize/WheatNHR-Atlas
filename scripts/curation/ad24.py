"""ad24 (Nov 2026): thirteen uploaded papers.
Added: Ban & Suenaga 2000 (Euphytica 113:87), Buerstmayr et al. 1997 (Cereal Res Commun 25:731), Bertucci et al. 2014 (Plant Dis 98:746),
Bokore et al. 2017 (Plant Dis 101:1292), Anh et al. 2015 (Phytopathology 105:1568), Burt & Nicholson 2011 (TAG 123:1387),
Ceoloni et al. 2017 (TAG 130:2005), Bernardo et al. 2012 (Mol Breed 29:477), Cao et al. 2001 (TAG 102:551).
Not added (no locus with a chromosome): Ban 2000 Breed Sci 50:131 (RAPD QTL on an unassigned linkage group), Ban & Suenaga 1997
(seminar report of the B1-linked Sumai 3 gene added from the 2000 paper; its RAPDs are unplaced), Ansari et al. 2007 TAG 114:927
(transcript study), Belhadj Benyahia et al. 2020 (Trichoderma biocontrol)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {
 'Ban2000Euph': ('Ban T, Suenaga K (2000) Genetic analysis of resistance to Fusarium head blight caused by Fusarium graminearum in Chinese wheat cultivar Sumai 3 and the Japanese cultivar Saikai 165. Euphytica 113:87-99', 'https://doi.org/10.1023/A:1003951509797', 'Ban and Suenaga (2000)'),
 'Buerstmayr1997CRC': ('Buerstmayr H, Lemmens M, Ruckenbauer P (1997) Chromosomal location of Fusarium head blight resistance genes in wheat. Cereal Res Commun 25:731-732', 'https://doi.org/10.1007/BF03543826', 'Buerstmayr et al. (1997)'),
 'Bertucci2014': ('Bertucci M, Brown-Guedira G, Murphy JP, Cowger C (2014) Genes conferring sensitivity to Stagonospora nodorum necrotrophic effectors in Stagonospora nodorum blotch-susceptible U.S. wheat cultivars. Plant Dis 98:746-753', 'https://doi.org/10.1094/PDIS-08-13-0820-RE', 'Bertucci et al. (2014)'),
 'Bokore2017': ('Bokore FE, Knox RE, DePauw RM, Clarke F, Cuthbert RD, Campbell HL et al. (2017) Validation of molecular markers for use with adapted sources of Fusarium head blight resistance in wheat. Plant Dis 101:1292-1299', 'https://doi.org/10.1094/PDIS-10-16-1421-RE', 'Bokore et al. (2017)'),
 'Anh2015Rmg8': ('Anh VL, Anh NT, Tagle AG, Vy TTP, Inoue Y, Takumi S et al. (2015) Rmg8, a new gene for resistance to Triticum isolates of Pyricularia oryzae in hexaploid wheat. Phytopathology 105:1568-1572', 'https://doi.org/10.1094/PHYTO-02-15-0034-R', 'Anh et al. (2015)'),
 'Burt2011Pch1': ('Burt C, Nicholson P (2011) Exploiting co-linearity among grass species to map the Aegilops ventricosa-derived Pch1 eyespot resistance in wheat and establish its relationship to Pch2. Theor Appl Genet 123:1387-1400', 'https://doi.org/10.1007/s00122-011-1674-9', 'Burt and Nicholson (2011)'),
 'Ceoloni2017': ('Ceoloni C, Forte P, Kuzmanovic L, Tundo S, Moscetti I, De Vita P et al. (2017) Cytogenetic mapping of a major locus for resistance to Fusarium head blight and crown rot of wheat on Thinopyrum elongatum 7EL and its pyramiding with valuable genes from a Th. ponticum homoeologous arm onto bread wheat 7DL. Theor Appl Genet 130:2005-2024', 'https://doi.org/10.1007/s00122-017-2939-8', 'Ceoloni et al. (2017)'),
 'Bernardo2012': ('Bernardo AN, Ma H, Zhang D, Bai G (2012) Single nucleotide polymorphism in wheat chromosome region harboring Fhb1 for Fusarium head blight resistance. Mol Breed 29:477-488', 'https://doi.org/10.1007/s11032-011-9565-y', 'Bernardo et al. (2012)'),
 'Cao2001SNB': ('Cao W, Hughes GR, Ma H, Dong Z (2001) Identification of molecular markers for resistance to Septoria nodorum blotch in durum wheat. Theor Appl Genet 102:551-554', 'https://doi.org/10.1007/s001220051681', 'Cao et al. (2001)'),
}
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in SRC.items()}
FHB1, B1 = 8.528, 698.5286
T = ' Read from the paper, Nov 2026 (ad24.py).'
def mk(i, **k):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand'); r.update(k)
    r['note'] = r['note'] + T; return r
def build(reports):
    keep = [r for r in reports if not re.match(r'S1(0[4-9]|1[0-2])-', r['id'])]; N = []
    N.append(mk('S104-001', name='FHB resistance gene linked to B1 (Sumai 3)', disease='Fusarium head blight', chr='5A', source='Ban2000Euph', peak=B1, grade='C', marker='B1 (awn inhibitor)',
        note='One of the two major additive FHB resistance genes of Sumai 3 (also in Saikai 165) is linked in repulsion to the awn-inhibitor gene B1 on 5AL: recombination 15.1 +/- 3.3% in Sumai 3/Gamenya DH lines and 21.4 +/- 4.3% in Sumai 3/Emblem DH lines (greenhouse and field tests, F. graminearum). Placed at the cloned B1 gene (TraesCS5A02G542800); the resistance gene lies about 15-21 cM from it.',
        evidence='Table 6: recombination between FHB resistance and awnedness 15.12% and 21.43%'))
    for n, (g, ch, eff, extra) in enumerate([('U-136.1', '5A', '+', ' Affected head symptoms and grain damage.'), ('U-136.1', '1B', '+', ''), ('U-136.1', '3B', '+', ' Significant for grain symptoms only.'),
            ('U-136.1', '4B', '+', ''), ('U-136.1', '6B', '+', ' Affected head symptoms and grain damage.'), ('U-136.1', '6D', '+', ' Affected head symptoms and grain damage.'),
            ('U-136.1', '2B', '-', ''), ('U-226.1', '3A', '+', ''), ('U-226.1', '3B', '+', ''), ('U-226.1', '6B', '+', ''), ('U-226.1', '4D', '+', ' Significant for grain symptoms only.')], 1):
        ped = 'Sagvari/Nobeokabozu-komugi//Mini-Mano/Sumai 3' if g == 'U-136.1' else 'Sumai 3/81.60//Kincso'
        N.append(mk(f'S105-{n:03d}', name=f'FHB chromosome effect {ch} ({g}; monosomic analysis)', disease='Fusarium head blight', chr=ch, source='Buerstmayr1997CRC', grade='D',
            note=f'Back-cross reciprocal monosomic analysis of Hungarian winter wheat line {g} ({ped}) against Hobbit sib, F. culmorum spray inoculation (Tulln). Chromosome {ch} of {g} '
                 + ('enhanced FHB resistance.' if eff == '+' else 'promoted susceptibility (negative effect).') + extra + ' Whole-chromosome effect: no position.', evidence=f'Table 1: {g} chromosome {ch} ({eff})'))
    N.append(mk('S106-001', name='Tsn1 (SNB-susceptible US winter wheats)', cls='DG', disease='Septoria nodorum blotch', chr='5B', source='Bertucci2014', start=546.698, end=546.843, grade='C', marker='Xfcp623',
        note='SnToxA sensitivity (Tsn1) in 32% of the NE-sensitive SNB-susceptible US winter wheat cultivars (infiltration with culture filtrates of 33 southeastern US isolates and SnToxA expressed in Pichia). Xfcp623 (internal to Tsn1) was the only diagnostic marker of ten tested. Survey of a known gene; placed at the Tsn1 position in the atlas.',
        evidence='Table 3; Table 2: Xfcp623 internal to Tsn1 (5BL)'))
    N.append(mk('S106-002', name='Snn3-B1 (SNB-susceptible US winter wheats)', cls='DG', disease='Septoria nodorum blotch', chr='5B', source='Bertucci2014', start=6.648, end=6.975, grade='C', marker='Xcfd20, Xgwm234',
        note='SnTox3 sensitivity (Snn3) in 64% of the NE-sensitive SNB-susceptible US winter wheat cultivars; none were sensitive to SnTox1 (Snn1) and Snn2 markers were not diagnostic. Markers Xcfd20 and Xgwm234 for Snn3-5B were not diagnostic. Survey of a known gene; placed at the Snn3-B1 position in the atlas.',
        evidence='Table 3: SnTox3-sensitive cultivars'))
    N.append(mk('S107-001', name='Fhb1 (Alsen, Carberry, ND3085, ND744)', disease='Fusarium head blight', chr='3B', source='Bokore2017', start=FHB1, end=FHB1, grade='C', marker='Xgwm493, Xgwm533',
        note='Validation in 14 doubled haploid breeding populations from Sumai 3-derived parents (FHB nurseries at Portage la Prairie and Carman, MB): Alsen, Carberry, ND3085 and ND744 carry Fhb1 as Sumai 3 does; no evidence for Glenn. Nursery selection raised the frequency of resistant Fhb1 marker alleles; Xgwm493 and Xgwm533 recommended for selection. Placed at the cloned Fhb1 (TaHRC) position in the atlas.',
        evidence='Tables 3-5: chi-square and t tests on Fhb1 markers'))
    p304, _ = resolve('Xgwm304', '5A')
    N.append(mk('S107-002', name='Fhb5 (Alsen, Carberry, ND3085, ND744)', disease='Fusarium head blight', chr='5A', source='Bokore2017', start=p304, end=214.2, grade='C', marker='Xgwm150, Xgwm304, Xgwm595',
        note='Validation in the same 14 DH populations: Alsen, Carberry, ND3085 and ND744 carry Fhb5 (Qfhs.ifa-5A) as Sumai 3 does; Xgwm150, Xgwm304 and Xgwm595 recommended for selection. Xgwm150 could not be placed and Xgwm595 places on distal 5AL, away from the pericentromeric Fhb5; placed at the Fhb5 interval in the atlas (Xgwm304 105.4 Mb to Xgwm415 214.2 Mb), a wide pericentromeric region.',
        evidence='Tables 3-5: chi-square and t tests on Fhb5 markers'))
    pw, hw = resolve('Xwmc317', '2B')
    N.append(mk('S108-001', name='Rmg8 (S-615; Anh 2015)', disease='Wheat blast', chr='2B', source='Anh2015Rmg8', peak=pw, grade='C', marker='Xwmc317-Xbarc159',
        note='New gene Rmg8 for resistance to Triticum isolate Br48 of Pyricularia oryzae in hexaploid wheat S-615 (seedlings and heads); F2 of S-615 x Hope (93 plants), flanked by Xwmc317 and Xbarc159 on distal 2BL. '
             f'Only Xwmc317 could be placed ({pw:.3f} Mb, {hw}); position is that marker. The cloned Rmg8 is absent from Chinese Spring; its closest homologue is at 796 Mb (atlas record R077).', evidence='Fig. 3: Rmg8 map, Xwmc317 and Xbarc159'))
    N.append(mk('S108-002', name='Rmg7 (St24, tetraploid; Anh 2015)', disease='Wheat blast', chr='2A', source='Anh2015Rmg8', peak=761.8578, grade='C', marker='Xcfd50-Xhbg327',
        note='Rmg7 in tetraploid St24 mapped with 94 F3 lines of St24 x Tat14, flanked by Xcfd50 and Xhbg327 on distal 2AL; Rmg7 and Rmg8 may be homoeologous. Neither marker could be placed; placed at the Rmg7 position in the atlas.',
        evidence='Fig. 3: Rmg7 map, Xcfd50 and Xhbg327'))
    N.append(mk('S109-001', name='Pch1 (Burt and Nicholson 2011)', disease='Eyespot', chr='7D', source='Burt2011Pch1', start=631.678, end=635.512, grade='C', marker='Xwg7S-Xcos7-9',
        note='Ae. ventricosa-derived Pch1 on 7DL located to a 0.26 cM interval between Xwg7S (with co-segregating Xcos7-6, X4CD7A8, Xorw5, Xtr383) and Xcos7-9 in a BC6F2 population (HS x RVPM25), using COS markers from Brachypodium chromosome 1 (364 kb region), rice 6 and sorghum 10. The COS and STS markers cannot be placed on Chinese Spring; placed at the Pch1 interval in the atlas (Xwmc14/Xbarc97/Xcfd175).',
        evidence='Results: Pch1 in 0.26 cM between Xwg7S and Xcos7-9'))
    pc, hc = resolve('Xcfa2040', '7A')
    N.append(mk('S109-002', name='Pch2 (Cappelle Desprez 7A; Burt and Nicholson 2011)', disease='Eyespot', chr='7A', source='Burt2011Pch1', peak=pc, grade='C', marker='Xcfa2040, Xorw1, Xcos7-11',
        note='Homoeologous COS markers extended the Pch2 QTL on 7AL distal to Xcfa2040 (LOD > 15.2, 14 cM region) in the CS x CS/CD7A F2 population; Pch1 and Pch2 are potentially homoeoloci. '
             f'Placed at Xcfa2040 ({pc:.3f} Mb, {hc}).', evidence='Results: Pch2 QTL LOD > 15.2 distal to Xcfa2040'))
    AL = ' The locus is on a Thinopyrum elongatum 7EL segment in wheat 7DL recombinants and has no counterpart position in Chinese Spring.'
    N.append(mk('S110-001', name='Fhb-7EL (FHB)', disease='Fusarium head blight', chr='7D', ref='alien (7EL/7DL recombinant)', source='Ceoloni2017', grade='D',
        note='Major FHB resistance locus Fhb-7EL mapped to the most distal (telomeric) end of Th. elongatum 7EL using 7EL-7el1L bread wheat recombinant lines (Th. ponticum 7el1L of line T4 on 7DL, carrying Lr19 and Sr25); reduced F. graminearum severity and fungal biomass in grains by over 95%. Likely the same locus as Fhb7.' + AL,
        evidence='Abstract; bioassays of 7EL-7el1L recombinants'))
    N.append(mk('S110-002', name='Fhb-7EL (crown rot)', disease='Fusarium crown rot', chr='7D', ref='alien (7EL/7DL recombinant)', source='Ceoloni2017', grade='D',
        note='The same Fhb-7EL locus was effective against F. culmorum and F. pseudograminearum, the main agents of crown rot.' + AL, evidence='Abstract; crown rot bioassays'))
    N.append(mk('S111-001', name='Fhb1 (Ning 7840 x Clark; SNP map)', disease='Fusarium head blight', chr='3B', source='Bernardo2012', peak=FHB1, grade='C', marker='Xsts3B-138, Xsnp3BS-11 (Xgwm533-Xgwm493)',
        note='Seven EST-derived SNP markers mapped between Xgwm533 and Xgwm493 in 71 Ning 7840/Clark BC7F7 RILs; five in a 1.1 cM marker-dense region gave the highest R2 (45-54%) and LOD (9.0-11.6), the most likely Fhb1 location; Xsts3B-138 and Xsnp3BS-11 had the largest R2. Placed at the cloned Fhb1 (TaHRC) position in the atlas.',
        stat='LOD 9.02-11.62', pve='45-54', evidence='Table 1: Xsnp3BS-11 LOD 11.62, R2 54.0%'))
    N.append(mk('S112-001', name='snbTM (durum, S12-1)', cls='DG', disease='Septoria nodorum blotch', chr='3A', source='Cao2001SNB', grade='D', marker='UBC521650, RC37510 (RAPD); SCAR from UBC521650',
        note='Single SNB resistance gene snbTM in durum line S12-1 (RILs of Sceptre/3-6, Sceptre/S9-10, Sceptre/S12-1); RAPD markers UBC521650 and RC37510 flank it at about 15 and 13 cM. RC37510 is also linked (37.1 cM) to the red seed colour gene, consistent with its location on 3A. RAPD markers cannot be placed.',
        evidence='Results: UBC521650 15 cM, RC37510 13 cM'))
    for r in N:
        r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}
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
    json.dump(ops, open('/tmp/ops24.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], r['asm'].get('RefSeq v2.1'))
    print(len(ops))
