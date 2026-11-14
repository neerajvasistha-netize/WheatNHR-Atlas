"""ad26 (Nov 2026): seven uploads. Added: Cruppe et al. 2023 Plant Dis 107:2407 (wheat blast QTL), Esvelt Klos et al. 2014 Plant Pathol
63:400 (Pch2), Dinglasan et al. 2019 TAG 132:149 (yellow spot GWAS), Dutta et al. 2021 Euphytica 217:108 (Stb6 GWAS; uploaded as the
preprint and the correction 217:175). Not added: Duguid & Brule-Babel 2001 (inheritance only, no mapping), Du et al. 2015 Plant Mol
Biol Rep 33:414 (Psathyrostachys huashanica 2Ns(2D) substitution line: alien chromosome, no wheat-genome position)."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {k: {'citation': c, 'url': u, 'short': s} for k, (c, u, s) in {
 'Cruppe2023': ('Cruppe G, Lemes da Silva C, Lollato RP, Fritz AK, Kuhnem P, Cruz CD, Calderon L, Valent B (2023) QTL pyramiding provides marginal improvement in 2NS-based wheat blast resistance. Plant Dis 107:2407-2416', 'https://doi.org/10.1094/PDIS-09-22-2030-RE', 'Cruppe et al. (2023)'),
 'EsveltKlos2014': ('Esvelt Klos KL, Wetzel HC III, Murray TD (2014) Resistance to Oculimacula yallundae and Oculimacula acuformis is conferred by Pch2 in wheat. Plant Pathol 63:400-404', 'https://doi.org/10.1111/ppa.12107', 'Esvelt Klos et al. (2014)'),
 'Dinglasan2019': ('Dinglasan EG, Singh D, Shankar M, Afanasenko O, Platz G, Godwin ID, Voss-Fels KP, Hickey LT (2019) Discovering new alleles for yellow spot resistance in the Vavilov wheat collection. Theor Appl Genet 132:149-162', 'https://doi.org/10.1007/s00122-018-3204-5', 'Dinglasan et al. (2019)'),
 'Dutta2021': ('Dutta A, Croll D, McDonald BA, Krattinger SG (2021) Genome-wide association study for septoria tritici blotch resistance reveals the occurrence and distribution of Stb6 in a historic Swiss landrace collection. Euphytica 217:108', 'https://doi.org/10.1007/s10681-021-02843-7', 'Dutta et al. (2021)'),
}.items()}
T = ' Read from the paper, Nov 2026 (ad26.py).'
def mk(i, **k):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand'); r.update(k); r['note'] += T; return r
def place(r, ms, ch):
    pl = [(m,) + resolve(m, ch) for m in ms]; pl = [x for x in pl if x[1]]
    if len(pl) >= 2:
        lo, hi = min(x[1] for x in pl), max(x[1] for x in pl); r['start'], r['end'] = lo, hi; r['grade'] = 'B' if hi - lo <= 100 else 'C'
        t = ' Placed in silico: ' + ', '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in pl) + '.'
        if r['grade'] == 'C': t += f' The markers spread over {hi - lo:.0f} Mb: a chromosome region, not a locus.'
    elif pl: m, p, h = pl[0]; r['peak'], r['grade'] = p, 'C'; t = f' Only {m} could be placed ({p:.3f} Mb, {h}); position is that marker.'
    else: r['grade'] = 'D'; t = ' None of the markers could be placed.'
    r['note'] = r['note'][:-len(T)] + t + T
CR = 'Wheat blast QTL in the TBIO Alvorada x TBIO Sossego population (4,169 SNPs; field Quirusillas, Bolivia 2017 and 2019; greenhouse Bolivia and Brazil). '
CQ = [('2AS (2NvS)', '2A', ['IWB9863', 'IWB36110'], 'Sossego', 'field 2017 severity LOD 17.5 (33.5%), incidence 23.0 (40.0%); field 2019 severity 44.7 (67.1%), incidence 52.6 (73.0%); greenhouse Bolivia severity 60.2 (75.8%), AUDPC 30.2 (52.7%); greenhouse Brazil severity 77.2 (82.4%), AUDPC 56.5 (74.1%)', '82.4',
       ' The locus corresponds to the Ae. ventricosa 2NvS translocation on 2AS; the markers are placed at their Chinese Spring positions on 2AS.'),
      ('5AL', '5A', ['IWB65128', 'IWB6827', 'IWA7789'], 'Sossego', 'field 2017 severity LOD 3.4 (5.5%), incidence 4.4 (6.0%); greenhouse Brazil AUDPC 3.1 (1.9%)', '6.0', ''),
      ('1DS', '1D', ['IWA4645'], 'Alvorada', 'field 2017 incidence LOD 3.3 (4.5%)', '4.5', ''),
      ('4DS', '4D', ['Xwmc285'], 'Sossego', 'greenhouse Brazil severity LOD 3.5 (1.3%)', '1.3', ''),
      ('5DL', '5D', ['IWB74430'], 'Sossego', 'greenhouse Bolivia severity LOD 2.8 (1.6%)', '1.6', ''),
      ('6DL', '6D', ['IWB6427'], 'Sossego', 'greenhouse Brazil severity LOD 3.9 (1.5%)', '1.5', '')]
DQ = [('1A.1', '1A', ['1150277', '1251432'], 'seedling; -log10 P up to 3.2, R2 up to 3.6%; aligns with tsc1', '3.6'),
      ('1A.2', '1A', ['995402'], 'adult plant (Field_Warwick 2016, API 2015); -log10 P up to 3.5, R2 up to 4.9%', '4.9'),
      ('2B.1', '2B', ['1167654', '1134162'], 'adult plant (Field_Warwick 2016); -log10 P up to 4.0, R2 up to 5.8%', '5.8'),
      ('2B.2', '2B', ['1109219', '1151114', '3021559'], 'seedling and adult; -log10 P up to 3.6, R2 up to 4.2%', '4.2'),
      ('4A', '4A', ['985921', '1107332', '1152450', '2288919'], 'seedling; -log10 P up to 3.5, R2 up to 4.1%', '4.1'),
      ('4B.1', '4B', ['1126690', '1129901'], 'seedling and adult; -log10 P up to 3.2, R2 up to 4.0%', '4.0'),
      ('4B.2', '4B', ['3534297', '978340'], 'seedling; -log10 P up to 3.3, R2 up to 3.8%', '3.8'),
      ('5B', '5B', ['4408847', '983884', '3955588', '3028160'], 'seedling; -log10 P up to 6.2, R2 up to 7.4%; aligns with tsn1 (Tsn1)', '7.4'),
      ('6B', '6B', ['1222803', '4542162', '4542467'], 'seedling; -log10 P up to 3.2, R2 up to 3.6%', '3.6'),
      ('6D', '6D', ['1223301', '1081184'], 'adult plant (Field_Warwick 2015, 2016); -log10 P up to 3.1, R2 up to 4.2%', '4.2'),
      ('7B', '7B', ['1130858', '1686218', '1091402', '1100309', '1140579'], 'seedling and adult; -log10 P up to 3.6, R2 up to 4.6%', '4.6')]
def build(reports):
    keep = [r for r in reports if not re.match(r'S1(1[89]|2[01])-', r['id'])]; N = []
    for n, (nm, ch, ms, par, det, mx, ex) in enumerate(CQ, 1):
        r = mk(f'S118-{n:03d}', name=f'QWB {nm} (Alvorada x Sossego)', disease='Wheat blast', chr=ch, source='Cruppe2023', marker=', '.join(ms), pve=f'up to {mx}',
               note=CR + f'Resistance allele from TBIO {par}. Detections: {det}.' + ex, evidence=f'Table 1: {nm}, peak markers {", ".join(ms)}')
        place(r, ms, ch); N.append(r)
    N.append(mk('S119-001', name='Pch2 (O. yallundae and O. acuformis)', disease='Eyespot', chr='7A', source='EsveltKlos2014', start=713.852, end=717.079, grade='C', marker='Xcfa2040 region',
        note='102 recombinant substitution lines of CS x CS/Cappelle Desprez 7A, seedling tests with GUS-transformed isolates: a QTL for resistance to O. acuformis (LOD 27 and 67) at the same 7AL position as resistance to O. yallundae (LOD 132 and 114), confirming that Pch2 confers resistance to both eyespot species. Placed at the Pch2 interval in the atlas.',
        evidence='Abstract and Table 1: LOD 27/67 (O. acuformis), 132/114 (O. yallundae)'))
    for n, (nm, ch, ms, det, mx) in enumerate(DQ, 1):
        r = mk(f'S120-{n:03d}', name=f'qNV.YS-{nm}', cls='MTA', disease='Tan spot', chr=ch, source='Dinglasan2019', marker=', '.join(ms), pve=f'up to {mx}',
               note=f'Yellow spot (tan spot) resistance QTL from GWAS of the N. I. Vavilov wheat collection (DArTseq markers; seedling tests with Ptr race 1 and field trials at Warwick, Qld): {det}.',
               evidence=f'Table 1: qNV.YS-{nm} DArTseq {", ".join(ms)}')
        place(r, ms, ch); N.append(r)
    N.append(mk('S121-001', name='Stb6 region (Swiss landraces)', cls='MTA', disease='Septoria tritici blotch', chr='3A', source='Dutta2021', peak=26.461994, grade='A', marker='Excalibur_c34649556', stat='p below Bonferroni 7.8e-09, MAF 0.12',
        note='GWAS of a historic Swiss landrace collection with four Swiss Z. tritici isolates: one SNP (Excalibur_c34649556) on 3A at 26,461,994 bp associated with necrosis and pycnidia for isolate 1E4, 255 kb from Stb6; 66% of landraces carry the functional Stb6-1 allele and seven new Stb6 alleles were found. Responses to the other isolates were quantitative with no significant SNP. Position as printed (RefSeq v1.0). Uploaded as the Research Square preprint; the published article and its correction (Euphytica 217:175) were consulted for the citation.',
        evidence='Results: Excalibur_c34649556 3A 26,461,994 bp'))
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
    json.dump(ops, open('/tmp/ops26.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    for r in new2: print(r['id'], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], r['asm'].get('RefSeq v2.1'))
