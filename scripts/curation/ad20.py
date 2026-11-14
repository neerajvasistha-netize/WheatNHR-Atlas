"""ad20 (Nov 2026): Nannuru et al. 2022 TAG 135:2247 (FHB, NMBU Nordic spring wheat panel, 25K and 35K SNP GWAS).
Table 3 prints 13 QTL regions on RefSeq v1.0 (Mbp). The atlas had 9: 8 with the right interval are regraded A with notes (same IDs);
Qfhb.nmbu.6A.1 carried the 6B.1 interval (630-688) and is replaced (new ID, 450-531); 1A.1, 3A.1, 3B.1 and 6A.2 are added."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19'); import lift
T = {  # name: (chr, start, end, -log10p, traits (n locations), panels, published)
 'Qfhb.nmbu.1A.1': ('1A', 9, 46, '3.05-13.27', 'AE (2), FHBcPHDH (3), DONcPHDH (4)', 'complete and European panels', 'previously reported (Jiang et al. 2007; Buerstmayr et al. 2009 and others)'),
 'Qfhb.nmbu.1A.2': ('1A', 520, 590, '3.04-8.16', 'AE (2), FHBcPHDH (2), DONcPHDH (4)', 'complete panel', 'previously reported (Venske et al. 2019; Zheng et al. 2020)'),
 'Qfhb.nmbu.3A.1': ('3A', 683, 738, '3.09-12.27', 'AE (4), FHBcPHDH (3), DONcPHDH (2)', 'complete and European panels', 'previously reported (Venske et al. 2019)'),
 'Qfhb.nmbu.3B.1': ('3B', 7, 9, '3.49-12.86', 'FHBcPHDH (2), DONcPHDH (1)', 'complete panel', 'the Fhb1 region (Anderson et al. 2001; Liu et al. 2006; Buerstmayr et al. 2009)'),
 'Qfhb.nmbu.4A.1': ('4A', 537, 607, '3.20-9.21', 'AE (1), FHBcPHDH (2), DONcPHDH (4)', 'complete panel', 'previously reported (Zheng et al. 2020)'),
 'Qfhb.nmbu.4A.2': ('4A', 721, 743, '3.09-8.40', 'AE (3), DONcPHDH (2)', 'complete and European panels', 'reported as new in this study'),
 'Qfhb.nmbu.4B.1': ('4B', 527, 609, '3.02-6.41', 'DONcPHDH (2)', 'European panel', 'reported as new in this study'),
 'Qfhb.nmbu.5A.1': ('5A', 480, 552, '3.0-12.49', 'AE (4), FHBcPHDH (2), DONcPHDH (2)', 'complete and European panels', 'previously reported (Venske et al. 2019; Zheng et al. 2020); significant haplotype effect on corrected FHB severity'),
 'Qfhb.nmbu.6A.1': ('6A', 450, 531, '3.06-6.66', 'AE (4), FHBcPHDH (1), DONcPHDH (1)', 'complete and European panels', 'reported as new in this study'),
 'Qfhb.nmbu.6A.2': ('6A', 608, 609, '3.09-5.73', 'AE (3), FHBcPHDH (4), DONcPHDH (4)', 'complete and European panels', 'previously reported (Ruan et al. 2020)'),
 'Qfhb.nmbu.6B.1': ('6B', 630, 688, '3.03-6.55', 'AE (4), FHBcPHDH (4), DONcPHDH (4)', 'complete and European panels', 'reported as new in this study'),
 'Qfhb.nmbu.7A.1': ('7A', 120, 129, '3.05-3.86', 'AE (1), FHBcPHDH (3), DONcPHDH (5)', 'complete and European panels', 'reported as new in this study'),
 'Qfhb.nmbu.7A.2': ('7A', 670, 710, '3.50-25.45', 'AE (2), FHBcPHDH (3), DONcPHDH (1)', 'complete and European panels',
                    'previously reported (Semagn et al. 2007; Buerstmayr et al. 2009 and others); the major locus: significant haplotype effects on corrected FHB severity and AE, '
                    'and the only region validated in new breeding lines; peak SNPs BS00098483_51 (698.69 Mb), AX-95248570 (701.36 Mb) and Kukri_c57593_79 (701.63 Mb) converted to KASP'),
}
KEEP = {'Qfhb.nmbu.7A.2': 'R116', 'Qfhb.nmbu.5A.1': 'R150', 'Qfhb.nmbu.1A.2': 'V1341', 'Qfhb.nmbu.4A.1': 'V1342', 'Qfhb.nmbu.4B.1': 'V1343',
        'Qfhb.nmbu.4A.2': 'V1344', 'Qfhb.nmbu.6B.1': 'V1346', 'Qfhb.nmbu.7A.1': 'V1347'}
OLD = ['V1345']
NEWNAMES = ['Qfhb.nmbu.1A.1', 'Qfhb.nmbu.3A.1', 'Qfhb.nmbu.3B.1', 'Qfhb.nmbu.6A.1', 'Qfhb.nmbu.6A.2']
SRC = 'Nannuru et al. 2022 TAG 135:2247 (Supplementary Table 14)'
PR = [{'m': f'{m} (KASP)', 'F': f'{a} / {b}', 'R': c, 'type': 'KASP', 'src': SRC} for m, a, b, c in [
 ('BS00098483_51', 'TCAGATAAGCAGCAGGGACAT', 'TCAGATAAGCAGCAGGGACAG', 'GAAAGGGAATTATACGGTCCAGA'),
 ('AX-95248570', 'TGGGACTGGATGTGGTGAG', 'TGGGACTGGATGTGGTGAA', 'GCAAAGCAATAGGGCTAGGA'),
 ('Kukri_c57593_79', 'CCACAGTAGGCTAAATGGACA', 'CCACAGTAGGCTAAATGGACG', 'ATCGTACACGCTCACTGCTG')]]
def fields(nm):
    ch, a, b, lp, tr, pan, pub = T[nm]
    f = dict(name=nm, cls='QTL', disease='Fusarium head blight', chr=ch, start=float(a), end=float(b), peak=None, ref='RefSeq v1.0',
             stat=f'-log10(p) {lp}', grade='A', curation='hand', source='NMBU2022FHB',
             note=(f'FHB QTL region in the NMBU Nordic spring wheat panel (GWAS with 25K and 35K SNP arrays; 300 lines, spawn-inoculated field trials); traits and number of '
                   f'environments: {tr} (AE = anther extrusion; FHB severity and DON corrected for plant height and days to heading); {pan}; {pub}. '
                   'Physical interval printed in Table 3 on RefSeq v1.0. Corrected from the paper, Nov 2026 (ad20.py).'),
             evidence=f'Table 3: {nm} {a}-{b} Mbp, -log10(p) {lp}, {tr}')
    if nm == 'Qfhb.nmbu.7A.2': f['primers'] = PR; f['marker'] = 'BS00098483_51, AX-95248570, Kukri_c57593_79 (peak SNPs, KASP)'
    return f
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S97-')]
    byid = {r['id']: r for r in keep}; upd = []
    for nm, i in KEEP.items():
        r = byid[i]; r.update(fields(nm)); r['asm'] = lift.lift_record(r); upd.append(r)
    new = []
    for n, nm in enumerate(NEWNAMES, 1):
        r = dict(id=f'S97-{n:03d}', **fields(nm)); r['asm'] = lift.lift_record(r); new.append(r)
    return keep + new, upd, new
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], _, _ = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], upd, new = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    V = {'R116': 2, 'R150': 2}
    ops = [{'op': 'set', 'path': 'reports/V1345', 'data': {'deleted': True}}]
    def doc(r):
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        return d
    for r in upd:
        d = doc(r)
        if r['id'] in V: ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'ifVersion': V[r['id']], 'data': d})
        else: ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    for r in new: ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': doc(r)})
    print(json.dumps(ops, ensure_ascii=False, separators=(',', ':')))
