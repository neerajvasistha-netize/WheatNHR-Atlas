"""ad31 (Nov 2026): five uploads. Added: Kharina et al. 2025 Agric Sci Euro-North-East 26:1298 (Stb2, Stb9, Stb11, Stb13 survey),
Kwiatek et al. 2012 J Plant Prot Res 52:248 (Pch1 survey), Leonard et al. 2008 TAG 116:261 (Ep-D1 candidate gene and STS markers for Pch1).
Not added: Lendenmann et al. 2016 Heredity 116:384 (QTL in the pathogen Z. tritici, not wheat), Khodarahmi et al. 2020 J Crop Breed 12:62
(inheritance of FHB resistance by generation means, no mapping)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Kharina2025': ('Kharina AV, Novoselova NV, Pakholkova EV (2025) Search for donors of Stb genes for resistance to septoria in collection samples of spring soft wheat. Agricultural Science Euro-North-East 26:1298-1308', 'https://doi.org/10.30766/2072-9081.2025.26.6.1298-1308', 'Kharina et al. (2025)'),
 'Kwiatek2012JPPR': ('Kwiatek M, Pankiewicz K, Wisniewska H, Korbas M, Danielewicz J (2012) Identification of Pch1 eyespot resistance gene in the collection of wheat lines (Triticum aestivum L.). J Plant Prot Res 52:254-258', 'https://doi.org/10.2478/v10045-012-0040-5', 'Kwiatek et al. (2012)'),
 'Leonard2008TAG': ('Leonard JM, Watson CJW, Carter AH, Hansen JL, Zemetra RS, Santra DK, Campbell KG, Riera-Lizarazu O (2008) Identification of a candidate gene for the wheat endopeptidase Ep-D1 locus and two other STS markers linked to the eyespot resistance gene Pch1. Theor Appl Genet 116:261-270', 'https://doi.org/10.1007/s00122-007-0664-4', 'Leonard et al. (2008)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad31.py).'
def mk(i, **k):
    r = dict(id=i, cls='DG', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand', grade='C'); r.update(k); r['note'] += T; return r
def build(reports):
    keep = [r for r in reports if not re.match(r'S14[2-4]-', r['id'])]; N = []
    KH = '146 spring bread wheat cultivars of the VIR world collection screened with SSR markers for Stb genes and tested with Z. tritici strains (seedlings, laboratory): Stb genes were found in 57.7% of cultivars; all four in Tepsey and Eminent; Stb13 may contribute most to resistance against the strains used. '
    for n, (g, ch, m, frq) in enumerate([('Stb2', '3B', 'Xgwm389', '14.3'), ('Stb9', '2B', 'Xwmc317', '55.9'), ('Stb11', '1B', 'Xbarc008', '67.9'), ('Stb13', '7B', 'Xwmc396', '79.8')], 1):
        p, h = resolve(m, ch)
        extra = ' Stb2 was screened with the 3BS markers Xgwm389 and Xgwm533 (Xgwm533 cannot be placed).' if g == 'Stb2' else ''
        N.append(mk(f'S142-{n:03d}', name=f'{g} (VIR spring wheat collection)', disease='Septoria tritici blotch', chr=ch, source='Kharina2025', peak=p, marker=m,
            note=KH + f'{g} occurred in {frq}% of the genotypes.' + extra + f' Survey of a known gene; placed at the diagnostic marker {m} ({p:.3f} Mb, {h}).', evidence=f'Abstract and Table 1: {g} {frq}%, marker {m}'))
    N.append(mk('S143-001', name='Pch1 (134 Polish breeding lines)', disease='Eyespot', chr='7D', source='Kwiatek2012JPPR', start=631.678, end=635.512, marker='EpD1b, XustSSR2001-7DL',
        note='134 wheat breeding lines (Randevous as resistant control) tested with the endopeptidase marker EpD1b, the STS marker XustSSR2001-7DL and field phenotyping: eight completely or highly resistant Pch1 lines were selected; marker results agreed 100% with the phenotype. Survey of a known gene; placed at the Pch1 interval in the atlas.',
        evidence='Abstract: eight resistant lines; markers 100% concordant'))
    pl = [(m,) + resolve(m, '7D') for m in ('Xwmc14', 'Xbarc97', 'Xcfd175')]
    N.append(mk('S144-001', name='Pch1 / Ep-D1 (Coda x Brundage)', disease='Eyespot', chr='7D', source='Leonard2008TAG', start=min(x[1] for x in pl), end=max(x[1] for x in pl), grade='B', marker='Xorw1 (Ep-D1 candidate), Xorw5, Xorw6, Xwmc14, Xbarc97, Xcfd175',
        note='In 254 Coda x Brundage RILs the STS marker Xorw1, from a wheat EST encoding oligopeptidase B (candidate for the endopeptidase locus Ep-D1), and STS markers Xorw5 and Xorw6, plus SSRs Xwmc14, Xbarc97 and Xcfd175, were completely linked to Ep-D1 and Pch1 on 7DL; the SSRs recombined in W7984 x Opata, suggesting reduced recombination near Pch1 in Coda/Brundage. In a variety panel Xorw1, Xorw5 and Xorw6 predicted Pch1 with 100% accuracy. '
             'Placed in silico across the co-segregating SSRs: ' + ', '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in pl) + '.', evidence='Results: complete linkage of Xorw1, Xorw5, Xorw6, Xwmc14, Xbarc97, Xcfd175 with Pch1'))
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
    json.dump(ops, open('/tmp/ops31.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    print(open('/tmp/ops31.json').read())
