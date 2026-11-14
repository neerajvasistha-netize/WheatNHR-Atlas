"""ad7 (Nov 2026): Hou et al. 2023 TAG 136:101 (Fusarium crown rot GWAS, 361 Chinese landraces).
Positions are on IWGSC RefSeq v2.0 (Methods: markers mapped to RefSeq v2.0); the atlas had them as v1.0/v1.x, unlifted,
two with no position and several taking a wrong locus's position. All 21 loci rebuilt from Table 4 (grade A)."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19')
import lift
OLD = ['R049','V1356','V1357','V1358','V1359','V1360','V1361','V1362']
# (QTL, [(marker, Mb, best -log10P, max PVE, trials)])
T = [
 ('Qfcr.sicau.1A-1', [('SNP59',14.33,3.16,4.86,'FCR2, FCR3, BLUP')]),
 ('Qfcr.sicau.1B-1', [('SNP1110',27.78,3.88,5.68,'FCR1'),('SNP1111',28.00,3.58,5.31,'FCR1'),('SNP1118',28.01,3.46,5.15,'FCR1'),('SNP1115',28.66,3.46,5.15,'FCR1')]),
 ('Qfcr.sicau.1B-2', [('SNP1607',530.24,4.15,6.02,'FCR1, FCR2, BLUP'),('SNP1610',532.25,3.40,5.08,'FCR1')]),
 ('Qfcr.sicau.1B-3', [('SNP1766',607.16,3.02,4.64,'FCR1, FCR3, BLUP')]),
 ('Qfcr.sicau.1B-4', [('SNP1868',641.36,4.41,6.28,''),('SNP1900',641.36,4.50,6.40,''),('SNP1869',641.74,5.98,8.29,'FCR1, FCR2, FCR3, BLUP'),
   ('SNP1901',641.74,5.98,8.29,''),('SNP1870',643.48,5.17,7.25,''),('SNP1879',644.77,3.19,4.82,''),('SNP1913',644.77,3.20,4.83,''),
   ('SNP1882',645.05,4.11,5.97,''),('SNP1917',645.05,3.49,5.19,''),('SNP1883',645.05,3.58,5.30,''),('SNP1884',645.13,3.25,4.84,''),('SNP1918',645.13,3.58,5.24,'')]),
 ('Qfcr.sicau.2B-1', [('SNP5583',794.74,3.34,4.95,'BLUP')]),
 ('Qfcr.sicau.2B-2', [('SNP5655',809.62,4.48,6.44,'FCR1, FCR2, BLUP')]),
 ('Qfcr.sicau.2D-1', [('SNP5957',90.58,3.53,5.24,'FCR1, BLUP')]),
 ('Qfcr.sicau.2D-2', [('SNP6065',577.59,3.08,4.68,'FCR1')]),
 ('Qfcr.sicau.3B-1', [('SNP7310',23.11,3.08,4.78,'FCR2')]),
 ('Qfcr.sicau.3B-2', [('SNP7561',126.92,3.50,5.20,'FCR1, BLUP')]),
 ('Qfcr.sicau.3B-3', [('SNP7718',262.57,3.24,4.88,'FCR1, FCR2, BLUP')]),
 ('Qfcr.sicau.3B-4', [('SNP7936',582.01,3.02,4.61,'FCR1')]),
 ('Qfcr.sicau.3D-1', [('SNP8533',5.30,3.29,4.95,'FCR1')]),
 ('Qfcr.sicau.4B-1', [('SNP9933',50.13,3.05,4.69,'FCR3')]),
 ('Qfcr.sicau.4B-2', [('SNP10193',628.82,3.68,5.61,'FCR2, BLUP'),('SNP10187',628.82,3.44,5.27,'FCR2'),('SNP10195',629.22,3.44,5.27,'FCR2'),('SNP10188',629.22,3.44,5.27,'FCR2'),('SNP10189',629.22,3.44,5.27,'FCR2')]),
 ('Qfcr.sicau.4B-3', [('SNP10261',657.40,3.36,5.14,'FCR3')]),
 ('Qfcr.sicau.5A-1', [('SNP11200',620.15,3.36,5.02,'FCR1')]),
 ('Qfcr.sicau.5B-1', [('SNP11592',76.60,3.18,4.88,'FCR3')]),
 ('Qfcr.sicau.7A-1', [('SNP16572',614.69,3.12,4.83,'FCR2')]),
 ('Qfcr.sicau.7B-1', [('SNP17158',49.03,3.22,4.94,'FCR3'),('SNP17159',49.03,3.30,5.04,'FCR3')]),
]
B4 = (' The major, stable locus: detected in all three repetitions and BLUP. KASP1869 (from SNP1869) was validated in an F2:3 population '
      '(136 lines), where the resistance allele explained up to 39.66% of phenotypic variance. Candidate genes TraesCS1B03G1106600 '
      '(exonic SNP) and TraesCS1B03G1088300 (NB-ARC) were induced at 4 dpi in the resistant landrace Shishoumai.')
PR = [{'m': 'KASP1869', 'F': 'GAAGGTGACCAAGTTCATGCTGCCCATTAAAGTAAAGCTAACGATA / GAAGGTGACCAAGTTCATGCTGCCCATTAAAGTAAAGCTAACGATG',
       'R': 'GCAGCAAGCCAACCAATACA', 'type': 'KASP', 'src': 'Hou et al. 2023 TAG 136:101 (Table S5)'}]
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S85-')]
    new = []
    for n, (nm, ms) in enumerate(T, 1):
        top = max(ms, key=lambda x: (x[2], x[3])); pos = [m[1] for m in ms]
        chr_ = re.search(r'\.(\d[ABD])-', nm).group(1)
        r = dict(id=f'S85-{n:03d}', name=nm, cls='QTL', disease='Fusarium crown rot', chr=chr_, ref='RefSeq v2.0',
                 start=min(pos) if len(set(pos)) > 1 else None, end=max(pos) if len(set(pos)) > 1 else None, peak=top[1],
                 marker=', '.join(m[0] for m in ms), pve=str(top[3]), stat=f'-log10(P) {top[2]}', source='Hou2023TAG', grade='A', curation='hand')
        tr = 'FCR1, FCR2, FCR3, BLUP' if nm.endswith('1B-4') else top[4]
        r['note'] = (f'GWAS locus for Fusarium crown rot (disease index) in 361 Chinese wheat landraces, compressed MLM; detected in {tr}. '
                     f'Positions printed in Table 4 on IWGSC RefSeq v2.0; peak is the most significant SNP ({top[0]}).'
                     + (f' {len(ms)} significant SNPs span the interval.' if len(ms) > 1 else '') + (B4 if nm.endswith('1B-4') else '')
                     + ' Corrected from the paper, Nov 2026 (ad7.py).')
        r['evidence'] = 'Table 4: ' + '; '.join(f'{m[0]} {chr_} {m[1]:.2f} Mb' for m in ms)
        if nm.endswith('1B-4'): r['primers'] = PR
        r['asm'] = lift.lift_record(r); new.append(r)
    return keep + new, new, OLD
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], new, rem = build(M['reports'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</'))
    Sd['reports'], new2, rem2 = build(Sd['reports'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    ops = [{'op': 'set', 'path': 'reports/' + i, 'data': {'deleted': True}} for i in rem2]
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/home/claude/work/ad7.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['peak'], r['start'], r['end'], r['asm']['RefSeq v2.1'])
