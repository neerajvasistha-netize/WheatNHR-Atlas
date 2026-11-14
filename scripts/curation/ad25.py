"""ad25 (Nov 2026): eight uploaded papers. Added: Chaubey et al. 2023 preprint (Tsn1, spot blotch), Choupannejad et al. 2023 preprint
(Tsn1, Snn1, Snn3 survey, SNB), Chartrain et al. 2005 Plant Pathol 54:134 (Stb6 survey), Chapman et al. 2009 TAG 118:1045 (Pch2
candidates), Chen et al. 2007 Mol Breed 19:113 (Qfhs.ndsu-3AS). Not added: Chu et al. 2006 TAG 112:1374 (hybrid necrosis genes Ne1/Ne2,
not disease resistance), Mol Breed 2017 37:141 (GWAS of agronomic traits only), Trop Plant Pathol 2025 50:10 (cultivar screening with
SSRs, no locus mapped)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Chaubey2023pre': ('Chaubey RK, Thakur D, Navathe S, Sharma S, Mishra VK, Singh PK, Chand R (2023) Heterologous expression and characterization of ToxA1 haplotype from India and its interaction with Tsn1 for spot blotch susceptibility in spring wheat. Research Square preprint', 'https://doi.org/10.21203/rs.3.rs-2796717/v1', 'Chaubey et al. (2023)'),
 'Choupannejad2023pre': ('Choupannejad R, Sharifnabi B, Collemare J, Massah A, Gholami J, Mehrabi R (2023) Host susceptibility and fungal virulence in the wheat-Parastagonospora nodorum interactions. Research Square preprint', 'https://doi.org/10.21203/rs.3.rs-3198290/v1', 'Choupannejad et al. (2023)'),
 'Chartrain2005PP': ('Chartrain L, Brading PA, Brown JKM (2005) Presence of the Stb6 gene for resistance to septoria tritici blotch (Mycosphaerella graminicola) in cultivars used in wheat-breeding programmes worldwide. Plant Pathol 54:134-143', 'https://doi.org/10.1111/j.1365-3059.2005.01164.x', 'Chartrain et al. (2005)'),
 'Chapman2009TAG': ('Chapman NH, Burt C, Nicholson P (2009) The identification of candidate genes associated with Pch2 eyespot resistance in wheat using cDNA-AFLP. Theor Appl Genet 118:1045-1057', 'https://doi.org/10.1007/s00122-009-0961-1', 'Chapman et al. (2009)'),
 'Chen2007MB': ('Chen X, Faris JD, Hu J, Stack RW, Adhikari T, Elias EM, Kianian SF, Cai X (2007) Saturation and comparative mapping of a major Fusarium head blight resistance QTL in tetraploid wheat. Mol Breed 19:113-124', 'https://doi.org/10.1007/s11032-006-9049-7', 'Chen et al. (2007)'),
}.items()}
TSN1, SNN1, SNN3, STB6 = (546.698, 546.843), (2.356, 2.357), (6.648, 6.975), 26.2035
T = ' Read from the paper, Nov 2026 (ad25.py).'
def mk(i, **k):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand', grade='C'); r.update(k); r['note'] += T; return r
def build(reports):
    keep = [r for r in reports if not re.match(r'S11[3-7]-', r['id'])]; N = []
    N.append(mk('S113-001', name='Tsn1 (ToxA1 of Bipolaris sorokiniana)', cls='DG', disease='Spot blotch', chr='5B', source='Chaubey2023pre', start=TSN1[0], end=TSN1[1],
        note='Preprint (not peer reviewed). ToxA1 cloned from Indian B. sorokiniana and expressed; infiltration of 18 Tsn1/tsn1 genotypes at seedling and adult stage confirmed the inverse gene-for-gene ToxA-Tsn1 interaction for spot blotch susceptibility (largest necrosis on K68, Tsn1; smallest on WAMI280, tsn1). Survey of a known gene; placed at the Tsn1 position in the atlas.',
        evidence='Results: necrosis area on Tsn1 and tsn1 genotypes'))
    for n, (g, pos, txt) in enumerate([('Tsn1', TSN1, 'present in 70% of highly susceptible cultivars and absent from all resistant ones; SnToxA infiltration caused necrosis only on Tsn1 cultivars'),
                                       ('Snn1', SNN1, 'SnTox1-Snn1 was the most frequent effector-gene pair in Iranian isolates and cultivars (100%); cultivars with Tsn1, Snn1 and Snn3 were susceptible to all isolates'),
                                       ('Snn3-B1', SNN3, 'Snn3 was detected in all Iranian cultivars and in 62.5% of world cultivars')], 1):
        ch = '5B' if g != 'Snn1' else '1B'
        N.append(mk(f'S114-{n:03d}', name=f'{g} (Iranian and world cultivars)', cls='DG', disease='Septoria nodorum blotch', chr=ch, source='Choupannejad2023pre', start=pos[0], end=pos[1],
            note=f'Preprint (not peer reviewed). Reaction of 40 wheat genotypes (19 Iranian cultivars, 16 world cultivars, five differentials) to 33 P. nodorum isolates and occurrence of necrotrophic effector sensitivity genes: {g} {txt}. Survey of a known gene; placed at the {g} position in the atlas.',
            evidence='Results: sensitivity gene frequencies'))
    N.append(mk('S115-001', name='Stb6 (cultivars worldwide)', disease='Septoria tritici blotch', chr='3A', source='Chartrain2005PP', peak=STB6, marker='Xgwm369',
        note='Resistance to M. graminicola isolate IPO323 in many cultivars used in wheat breeding worldwide was attributed to Stb6; eleven Xgwm369 alleles (closely linked to Stb6) were found, and pedigrees combined with Xgwm369 alleles indicate several independent entries of Stb6 into breeding programmes. Survey of a known gene; placed at the Stb6 position in the atlas.',
        evidence='Table 1; Xgwm369 allele analysis'))
    N.append(mk('S116-001', name='Pch2 candidate region (cDNA-AFLP)', disease='Eyespot', chr='7A', source='Chapman2009TAG', start=713.852, end=717.079, marker='4CD7A8, 19CD7A4',
        note='cDNA-AFLP of Chinese Spring versus CS/Cappelle Desprez 7A (Pch2) with Oculimacula acuformis: 34 of about 4,700 fragments differed; 14 clones mapped to 7A and three in the Pch2 region, notably 4CD7A8 (homologous to a rice callose synthase) and 19CD7A4 (homologous to the cereal cyst nematode NBS-LRR protein RCCN). Most differences appear to be allelic rather than expression differences. Placed at the Pch2 interval in the atlas.',
        evidence='Results: fragments mapped to the Pch2 region'))
    p, h = resolve('Xgwm2', '3A')
    N.append(mk('S117-001', name='Qfhs.ndsu-3AS (T. dicoccoides; saturation map)', disease='Fusarium head blight', chr='3A', source='Chen2007MB', peak=p, marker='Xgwm2 (Xfcp401-Xfcp397.2)', pve='42',
        note='Type II FHB resistance QTL from T. dicoccoides chromosome 3A in 83 Langdon x LDN(DIC-3A) RICLs: Xgwm2 most significant (P < 0.00001, 38%); CIM explained 42% within an 11.5 cM interval flanked by TRAP loci Xfcp401 and Xfcp397.2; not homoeologous to Qfhs.ndsu-3BS. '
             f'TRAP markers cannot be placed; placed at the peak marker Xgwm2 ({p:.3f} Mb, {h}).', evidence='Results: CIM, Xgwm2 R2 38%, QTL 42%'))
    for r in N: r['asm'] = lift.lift_record(r)
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
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1',) if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/tmp/ops25.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['asm'].get('RefSeq v2.1'))
