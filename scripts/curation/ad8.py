"""ad8 (Nov 2026): Lin et al. 2023 Plant Disease 107:1151 (FCR GWAS, 209 western-China landraces, 55K array).
Positions are on RefSeq v2.0 (Methods); the atlas held 4 of the 10 stable loci as RefSeq v1.x (shifted 0.5-8 Mb on v2.1) and two
auto records (V1349 2A unplaced, V1350 7D at the 7B position) under the wrong source. All 10 Table 2 loci rebuilt (grade A);
peak SNP and detection set from Table S4; KASP primers from Table S6."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19')
import lift
OLD = ['R117', 'R118', 'R119', 'R154', 'V1349', 'V1350']
L = json.load(open('/tmp/pd22/loci.json'))
T = [('Qfcr.sicau.1B', 303.51, 323.69, '3.00-3.61', '5.70-7.19'), ('Qfcr.sicau.2A-1', 68.68, 75.88, '3.00-5.49', '6.01-14.48'),
     ('Qfcr.sicau.2A-2', 107.35, 110.62, '3.08-3.61', '6.24-8.15'), ('Qfcr.sicau.2A-3', 212.98, None, '3.91-5.21', '9.76-13.11'),
     ('Qfcr.sicau.2A-4', 536.72, None, '3.25-3.97', '6.46-8.13'), ('Qfcr.sicau.2A-5', 563.26, None, '3.19-4.00', '6.27-8.37'),
     ('Qfcr.sicau.5A', 478.19, None, '3.37-4.20', '8.19-10.29'), ('Qfcr.sicau.7A', 95.10, None, '3.38-3.77', '6.69-7.63'),
     ('Qfcr.sicau.7B', 685.66, 686.90, '3.22-4.06', '7.13-9.28'), ('Qfcr.sicau.7D', 90.28, 95.48, '3.00-4.41', '5.76-12.21')]
SRC = 'Lin et al. 2023 Plant Dis 107:1151 (Table S6); FAM tail GAAGGTGACCAAGTTCATGCT, HEX tail GAAGGTCGGAGTCAACGGATT'
PR = {'Qfcr.sicau.2A-1': [{'m': 'AX-109889997 (KASP)', 'F': 'agttgtataccaaaatggcatcA (FAM, resistance) / agttgtataccaaaatggcatcG (HEX)', 'R': 'gggtcaatcacgggcatgta', 'type': 'KASP', 'src': SRC}],
      'Qfcr.sicau.7D': [{'m': 'AX-89512658 (KASP)', 'F': 'gctatcggtcgccatctcG (FAM, resistance) / gctatcggtcgccatctcA (HEX)', 'R': 'agagcagagcagggaaactt', 'type': 'KASP', 'src': SRC}]}
EX = {'Qfcr.sicau.2A-1': ' Major locus; KASP marker AX-109889997 validated in an F6 RIL population, where lines with the resistance allele had 17.78% lower FCR severity.',
      'Qfcr.sicau.2A-3': ' Major locus (haplotype analysis).', 'Qfcr.sicau.5A': ' Major locus (haplotype analysis).',
      'Qfcr.sicau.7D': ' Major locus, reported as novel on 7DS; KASP marker AX-89512658 developed.'}
def build(reports):
    keep = [r for r in reports if r['id'] not in OLD and not r['id'].startswith('S86-')]
    new = []
    for n, (nm, a, b, lp, pv) in enumerate(T, 1):
        d = L[nm]; ch = nm.split('.')[-1][:2]
        r = dict(id=f'S86-{n:03d}', name=nm, cls='QTL', disease='Fusarium crown rot', chr=ch, ref='RefSeq v2.0',
                 start=a if b else None, end=b, peak=d['pos'], marker=d['pk'] + (f' (peak of {d["n"]} significant SNPs)' if d['n'] > 1 else ''),
                 pve=pv, stat=f'-log10(P) {lp}', source='PD2023FCR', grade='A', curation='hand')
        r['note'] = (f'Stable GWAS locus for Fusarium crown rot (disease index) in 209 accessions of four western-China landrace subspecies '
                     f'(Tibetan semi-wild, Yunnan hulled, Xinjiang rice and Sichuan white wheat), 55K SNP array; detected in {", ".join(d["tr"])}. '
                     f'Positions printed in Table 2 on RefSeq v2.0; peak is the most significant SNP in Table S4 ({d["pk"]}, {d["pos"]} Mb, '
                     f'-log10(P) {d["lp"]}, PVE {d["pve"]}%).' + EX.get(nm, '') + ' Corrected from the paper, Nov 2026 (ad8.py).')
        r['evidence'] = f'Table 2: {nm} {ch} ' + (f'{a}-{b}' if b else f'{a}') + f' Mb; -log10(P) {lp}; PVE {pv}%'
        if nm in PR: r['primers'] = PR[nm]
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
    V = {'R117': 2, 'R118': 2, 'R119': 2, 'R154': 2}
    ops = [dict({'op': 'set', 'path': 'reports/' + i, 'data': {'deleted': True}}, **({'ifVersion': V[i]} if i in V else {})) for i in rem2]
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/home/claude/work/ad8.json', 'w'), ensure_ascii=False)
    for r in new: print(r['id'], r['name'], r['peak'], r['start'], r['end'], r['asm']['RefSeq v2.1'])
