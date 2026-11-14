"""ad12 (Nov 2026): Li et al. 2024 TAG 137:43 (FCR GWAS, 223 accessions, 90K). Table 4 positions are on CS RefSeq v1.1 (= v1.0
coordinates); the atlas had 35 per-SNP auto records (V0959-V0991, V1309, V1310), many labelled RefSeq v2.1 with values that are
neither the paper's nor v2.1 (e.g. 3A 748.82 for 746.63), two loci missing (3D-1, 3D-2 unplaced). Rebuilt as the 10 Table 4 QTL,
grade A; peak = SNP with the highest BLUP -log10(P). Candidate-gene record Z0256 (TraesCS3D02G524700) kept."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19')
import lift
OLD = ['V%04d' % i for i in range(959, 992)] + ['V1309', 'V1310']
# QTL: [(SNP, Mb, BLUP -log10P, BLUP PVE, trials)]
T = {
 'Qfcr.cau.2B': [('Ex_c55735_1012', 605.00, 3.06, 5.62, 'CRS1, CRS2, BLUP'), ('Kukri_rep_c112675_676', 605.00, 3.26, 3.00, 'CRS1, CRS2, BLUP'), ('wsnp_Ex_c55735_58127324', 605.00, 3.17, 3.00, '')],
 'Qfcr.cau.3A': [('BobWhite_c5246_196', 746.63, 3.33, 5.47, 'CRS1, CRS2, BLUP'), ('RAC875_c4841_753', 747.61, 3.35, 5.12, 'CRS1, CRS2, BLUP')],
 'Qfcr.cau.3B': [('CAP7_c3438_196', 821.08, 3.24, 3.13, 'CRS1, CRS2, BLUP'), ('IAAV686', 822.11, 3.12, 6.11, '')],
 'Qfcr.cau.3D-1': [('RFL_Contig5322_219', 37.11, 4.88, 7.40, 'CRS1, CRS2, CRS3, BLUP')],
 'Qfcr.cau.3D-2': [('BS00062806_51', 517.38, 4.17, 6.55, 'CRS1, CRS2, CRS3, BLUP')],
 'Qfcr.cau.3D-3': [('wsnp_Ex_c12369_19730765', 609.12, 0, 0, ''), ('BS00070060_51', 614.69, 0, 0, ''), ('wsnp_Ex_c12963_20529964', 612.90, 5.67, 4.36, 'CRS1, CRS2, CRS3, BLUP')],
 'Qfcr.cau.4A': [('Tdurum_contig75819_559', 712.86, 3.30, 5.08, 'CRS1, CRS3, BLUP')],
 'Qfcr.cau.5D': [('D_GCE8AKX02IXEFJ_281', 548.86, 3.61, 5.52, ''), ('BS00022688_51', 548.87, 3.70, 5.67, 'CRS1, CRS2, CRS3, BLUP'), ('RFL_Contig4176_605', 548.88, 3.56, 5.55, ''), ('D_GDS7LZN02I3554_251', 548.94, 3.61, 5.67, '')],
 'Qfcr.cau.7A': [('BS00069242_51', 32.62, 3.28, 5.37, 'CRS1, CRS3, BLUP')],
 'Qfcr.cau.7B': [('IAAV5530', 15.30, 3.02, 5.00, 'CRS1, CRS3, BLUP')],
}
EX = {
 'Qfcr.cau.3D-3': (' The major locus on 3DL: 38 SNPs (609.12-614.69 Mb) detected in all three trials and BLUP, -log10(P) up to 5.67, PVE up to '
                   '9.20%. Validated with KASP3D61373 in two RIL populations (Doumai/Shi 4185 and Linmai 2/Zhongmai 892): lines with the '
                   'resistance allele had 21.55% and 21.39% lower disease index. TraesCS3D02G524700 (PR-4) is proposed as a strong candidate gene.'),
 'Qfcr.cau.3A': ' Sequences of its SNPs are collinear with the Qfcr.cau.3D-3 region, so the loci may be homoeologous.',
 'Qfcr.cau.3B': ' Sequences of its SNPs are collinear with the Qfcr.cau.3D-3 region, so the loci may be homoeologous.',
}
PR = {'Qfcr.cau.3D-3': [{'m': 'KASP3D61373', 'F': 'ACCTCGTACTGATCTTGGAATTCG (FAM) / ACCTCGTACTGATCTTGGAATTCA (HEX)', 'R': 'GAGCAGTACGAGGTCGGC',
       'type': 'KASP', 'src': 'Li et al. 2024 TAG 137:43 (Table S2; allele primers shown without the FAM/HEX tails)'}]}
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S89-')]
    new = []
    for n, (nm, ms) in enumerate(T.items(), 1):
        top = max(ms, key=lambda x: (x[2], x[3])); pos = [m[1] for m in ms]; ch = re.search(r'\.(\d[ABD])', nm).group(1)
        multi = len(set(pos)) > 1
        r = dict(id=f'S89-{n:03d}', name=nm, cls='QTL', disease='Fusarium crown rot', chr=ch, ref='RefSeq v1.0',
                 start=min(pos) if multi else None, end=max(pos) if multi else None, peak=top[1],
                 marker=('38 SNPs; peak ' + top[0]) if nm.endswith('3D-3') else ', '.join(m[0] for m in ms),
                 pve=('up to 9.20' if nm.endswith('3D-3') else str(top[3])), stat=f'-log10(P) {top[2]} (BLUP)', source='Z_1244', grade='A', curation='hand')
        r['note'] = (f'GWAS locus for Fusarium crown rot (disease index, F. pseudograminearum) in 223 wheat accessions, 90K SNP array; detected in {top[4]}. '
                     f'Positions printed in Table 4 on Chinese Spring RefSeq v1.1 (same coordinates as v1.0); peak is the SNP with the highest BLUP -log10(P) ({top[0]}).'
                     + EX.get(nm, '') + ' Corrected from the paper, Nov 2026 (ad12.py).')
        r['evidence'] = 'Table 4: ' + ('3D 609.12-614.69 Mb, 38 SNPs' if nm.endswith('3D-3') else '; '.join(f'{m[0]} {ch} {m[1]:.2f} Mb' for m in ms))
        if nm in PR: r['primers'] = PR[nm]
        r['asm'] = lift.lift_record(r); new.append(r)
    return keep + new, new, OLD
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], new, rem = build(M['reports']); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], new2, rem2 = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    ops = [{'op': 'set', 'path': 'reports/' + i, 'data': {'deleted': True}} for i in rem2]
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/home/claude/work/ad12.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['peak'], r['start'], r['end'], r['asm']['RefSeq v2.1'])
