"""ad27 (Nov 2026): eight uploads. Added: Friesen & Faris 2004 TAG 109:464 (tsc2, race 5 QTL), Faris et al. 1997 TAG 94:98 (tan spot QTL,
ITMI), Faris et al. 1999 TAG 98:219 (candidate defence genes at tan spot and leaf rust QTL, ITMI), Forte et al. 2014 Mol Breed 34:1701
(Fhb-7el2, alien), Gholizadeh Vazvani et al. 2025 Mol Breed 45:36 (take-all marker association).
Not added: Ferjaoui et al. 2021 preprint (already in the atlas from the published BMC Genomics 2022 version; its Additional Table 4 with
physical positions was not uploaded), Gamba & Lamari 1998 (inheritance only), Faris et al. 2014 MGG 289:641 (spike compactness, not disease).
ITMI RFLP probes cannot be placed on Chinese Spring; QTL marked only by them are grade D."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Friesen2004TAG': ('Friesen TL, Faris JD (2004) Molecular mapping of resistance to Pyrenophora tritici-repentis race 5 and sensitivity to Ptr ToxB in wheat. Theor Appl Genet 109:464-471', 'https://doi.org/10.1007/s00122-004-1678-9', 'Friesen and Faris (2004)'),
 'Faris1997TAG': ('Faris JD, Anderson JA, Francl LJ, Jordahl JG (1997) RFLP mapping of resistance to chlorosis induction by Pyrenophora tritici-repentis in wheat. Theor Appl Genet 94:98-103', 'https://doi.org/10.1007/s001220050387', 'Faris et al. (1997)'),
 'Faris1999TAG': ('Faris JD, Li WL, Liu DJ, Chen PD, Gill BS (1999) Candidate gene analysis of quantitative disease resistance in wheat. Theor Appl Genet 98:219-225', 'https://doi.org/10.1007/s001220051061', 'Faris et al. (1999)'),
 'Forte2014MB': ('Forte P, Virili ME, Kuzmanovic L, Moscetti I, Gennaro A, D\'Ovidio R, Ceoloni C (2014) A novel assembly of Thinopyrum ponticum genes into the durum wheat genome: pyramiding Fusarium head blight resistance onto recombinant lines previously engineered for other beneficial traits from the same alien species. Mol Breed 34:1701-1716', 'https://doi.org/10.1007/s11032-014-0175-3', 'Forte et al. (2014)'),
 'Gholizadeh2025MB': ('Gholizadeh Vazvani M, Dashti H, Saberi Riseh R et al. (2025) Association analysis of response to take-all disease with agronomic traits and molecular markers and selection ideal genotypes in bread wheat (Triticum aestivum L.) genotypes. Mol Breed 45:36', 'https://doi.org/10.1007/s11032-025-01554-4', 'Gholizadeh Vazvani et al. (2025)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad27.py).'
RF = ' The RFLP probe cannot be placed on Chinese Spring, so the locus is unplaced.'
def mk(i, **k):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand'); r.update(k); r['note'] += T; return r
IT = 'W-7984 x Opata 85 (ITMI) recombinant inbred population'
def build(reports):
    keep = [r for r in reports if not re.match(r'S12[2-6]-', r['id'])]; N = []
    FR = f'Resistance to Ptr race 5 (conidial inoculation) and sensitivity to Ptr ToxB in the {IT}. '
    N.append(mk('S122-001', name='tsc2 (Ptr ToxB insensitivity, ITMI)', cls='DG', disease='Tan spot', chr='2B', source='Friesen2004TAG', start=24.09, end=25.02, grade='C', marker='Xfbb274 (distal 2BS)', pve='69',
        note=FR + 'The ToxB-insensitivity gene, designated tsc2, mapped to the distal tip of 2BS, 4.3 cM from Xfbb274, and underlay the major QTL for resistance to race 5 (69% of the variation; resistance from Opata 85). Placed at the Tsc2 position in the atlas (Corsi et al. 2020).',
        evidence='Results; Fig. 4: tsc2 on 2BS'))
    for n, (ch, m, par, r2, ex) in enumerate([('4A', 'Xksu916(Oxo)', 'W-7984', '20.0', ' Oxalate oxidase candidate locus on 4AL; with tsc2 explains 73% of the variation.'),
                                              ('2B', 'Xmwg2025', 'Opata 85', '11.9', ' On 2BL, 44 cM or more from tsc2.'), ('2A', 'Xcdo447', 'Opata 85', '11.7', ' On 2AS.')], 2):
        N.append(mk(f'S122-{n:03d}', name=f'QTs.fcu-{ch} (race 5, ITMI)', disease='Tan spot', chr=ch, source='Friesen2004TAG', grade='D', marker=m, pve=r2,
            note=FR + f'Minor QTL at {m}, R2 {float(r2):.0f}%, resistance from {par}.' + ex + RF, evidence=f'Table 1: {m} R2 {r2}%'))
    FA = f'Seedling resistance to tan spot isolate Pti2 (nec+ chl+) in the {IT} (532-marker RFLP map). '
    N.append(mk('S123-001', name='QTsc.ndsu-1A (ITMI, Pti2)', disease='Tan spot', chr='1A', ref='RefSeq v2.1', start=4.6, end=4.76, grade='C', marker='XGli1', pve='35',
        source='Faris1997TAG', note=FA + 'Major QTL on 1AS at XGli1 (R2 35%, resistance from W-7984), with an epistatic interaction with Xbcd120 on 2DL (11.5%); with Xfba211 (4AL) the model explained 49% of the variation. This 1AS chlorosis-resistance locus of W-7984 is the Tsc1 region; placed at the Tsc1 position in the atlas (Running et al. 2022).',
        evidence='Table 2: XGli1 1AS R2 0.350'))
    for n, (ch, m, par, r2) in enumerate([('4A', 'Xfba211', 'W-7984', '10.7'), ('4A', 'Xfbb1', 'W-7984', '13.7'), ('3B', 'Xpsr903', 'W-7984', '8.6'), ('1B', 'Xcdo346', 'Opata 85', '7.6')], 2):
        arm = {'Xfba211': '4AL', 'Xfbb1': '4AS (probed on 56 lines only)', 'Xpsr903': '3BL (extensive chlorosis; not kept in the multiple regression)', 'Xcdo346': '1BL (extensive chlorosis; not kept in the multiple regression)'}[m]
        N.append(mk(f'S123-{n:03d}', name=f'QTs.ndsu-{ch} at {m} (ITMI, Pti2)', disease='Tan spot', chr=ch, source='Faris1997TAG', grade='D', marker=m, pve=r2,
            note=FA + f'QTL at {m} on {arm}, R2 {r2}%, resistance from {par}.' + RF, evidence=f'Table 2: {m} R2 {r2}%'))
    FC = f'Candidate defence-response gene loci mapped in the {IT} and tested against tan spot (Faris et al. 1997) and leaf rust (Nelson et al. 1997) phenotypes. '
    C = [('S124-001', 'Oxo2/1433a (tan spot, 4AL)', 'Tan spot', '4A', 'Oxo2-4A, 1433a-4A', 'W-7984', 'R2 17.5% (infection type) and 13.4% (leaf area diseased); oxalate oxidase and 14-3-3 genes refine the 4AL QTL first seen at Xfba211', None),
         ('S124-002', 'Oxo2/1433a (leaf rust, 4AL)', 'Leaf rust', '4A', 'Oxo2-4A, 1433a-4A', 'Opata 85', 'R2 6.2-11.2% across seedling, adult and field tests', None),
         ('S124-003', 'Lr34 region with Pr1b/Mpc1 (leaf rust, 7DS)', 'Leaf rust', '7D', 'Xwg834 (Lr34), Pr1b, Mpc1', 'Opata 85', 'Xwg834 R2 up to 42.8% (adult plants and field); Pr1b/Mpc1 R2 8.1%', (47.412062, 47.42449)),
         ('S124-004', 'Cht1b/Tha1/Cat cluster (leaf rust, 7BL)', 'Leaf rust', '7B', 'Cht1b, Tha1, Cat', 'Opata 85', 'R2 11.0-42.2% under natural infection; two closely linked QTL on 7BL', None),
         ('S124-005', 'Per2 (Lr23?) (leaf rust, 2BS)', 'Leaf rust', '2B', 'Per2-2B', 'W-7984', 'R2 33.1% (seedlings, pathotype TBD/TM) and 15.7% (field); peroxidase locus near the reported Lr23 position, 10 cM distal of Xtam72', None),
         ('S124-006', 'Pal near Lr27 (leaf rust, 3BS)', 'Leaf rust', '3B', 'XksuG53 (Lr27), Pal-3B', 'Opata 85', 'XksuG53 R2 16.2%, Pal R2 13.5% (seedlings, pathotype BBB/BN)', None),
         ('S124-007', 'Tha1 (leaf rust, 7DL)', 'Leaf rust', '7D', 'Tha1-7D', 'Opata 85', 'R2 16.8% (seedlings, pathotype BBB/BN); thaumatin may itself contribute', None)]
    for i, nm, dis, ch, m, par, det, pos in C:
        r = mk(i, name=nm + ' (ITMI)', disease=dis, chr=ch, source='Faris1999TAG', marker=m, grade='C' if pos else 'D',
               note=FC + f'Resistance from {par}: {det}.' + (' Placed at the cloned Lr34 gene (TraesCS7D02G080300).' if pos else ' The candidate-gene probes cannot be placed on Chinese Spring, so the locus is unplaced.'), evidence=f'Table 1: {m}')
        if pos: r['start'], r['end'] = pos
        N.append(r)
    N.append(mk('S125-001', name='Fhb-7el2 (Th. ponticum 7el2L in durum)', disease='Fusarium head blight', chr='7A', ref='alien (7el1L+7el2L/7AL recombinant)', source='Forte2014MB', grade='D',
        note='FHB resistance QTL Fhb-7el2 from the Thinopyrum ponticum 7el2L arm, transferred from a bread wheat 7DS.7el2L translocation into durum 7AS.7AL-7el1L recombinant lines (carrying Lr19 and Sr25 from another Th. ponticum accession); the 7el1L+7el2L tetraploid lines had 70-85% lower FHB severity and lower fungal biomass in kernels; similar effects in bread wheat lines with the same assembly on 7DL. The locus is on an alien segment and has no counterpart position in Chinese Spring.',
        evidence='Abstract; FHB tests of 7el1L+7el2L recombinants'))
    p, h = resolve('Xbarc197', '5A')
    N.append(mk('S126-001', name='Xbarc197 association (take-all, 100 bread wheats)', cls='MTA', disease='Take-all', chr='5A', source='Gholizadeh2025MB', peak=p, grade='C', marker='Xbarc197 (allele 7)',
        note='Stepwise regression of take-all (Gaeumannomyces tritici) disease severity on SSR, IRAP and wheat-rye translocation markers in 100 bread wheat genotypes (greenhouse): DS = 67.03 - 18.23 Xbarc56-1 + 14.14 Xbarc197-7 - 25.64 LTR14 (R2 21%); allele 7 of Xbarc197 increased severity. Xbarc56 and the IRAP marker LTR14 (the strongest, lowering severity) cannot be placed and their chromosomes are not given. '
             f'Placed at Xbarc197 ({p:.3f} Mb, {h}).', evidence='Table 6: stepwise regression for DS'))
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
    json.dump(ops, open('/tmp/ops27.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], r['asm'].get('RefSeq v2.1'))
