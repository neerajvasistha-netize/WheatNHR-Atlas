"""ad23 (Nov 2026): six uploaded papers.
- Ballini et al. 2020 Euphytica 216:92 (durum EPO panel GWAS; SNB, wheat blast, STB, leaf rust). Positions in the paper are on the
  wild emmer Zavitan v2.0 assembly. Loci the paper ties to a Chinese Spring gene are placed at that gene (Snn1, Rmg7, Rmg8,
  TraesCS3A01G041800); the rest use the Zavitan position as an approximation of the Chinese Spring position (grade C).
  LR mock- and chitin-treated detections at the same SNP are one record each.
- Arraiano et al. 2007 TAG 116:113 (Hobbit sib 5BS susceptibility; Stb6 in Bezostaya 1).
- Agostinelli et al. 2012 Euphytica 186:115 (Fhb1 and QFhs.nau-2DL effects in KY93C-1238-17-2 x VA01W-476).
- Arruda et al. 2016 Mol Breed 36:84 reprints the GWAS MTAs of Arruda et al. 2016 Plant Genome (A_FHB_024); the significant ones
  are added to that original study. CSS scaffold SNP names cannot be placed; MTAs in LD with Fhb1 are placed at Fhb1.
- Not added: An et al. 2019 TAG (rye 4R addition line) and Anderson et al. 2010 Euphytica (Lophopyrum substitution lines):
  resistance is on alien chromosomes with no position on the wheat genome."""
import json, re, sys
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
SRC = {
 'Ballini2020': {'citation': 'Ballini E, Tavaud M, Ducasse A, Sanchez D, Paillard S, Morel JB et al. (2020) Genome wide association mapping for resistance to multiple fungal pathogens in a panel issued from a broad composite cross-population of tetraploid wheat Triticum turgidum. Euphytica 216:92',
                 'url': 'https://doi.org/10.1007/s10681-020-02631-9', 'short': 'Ballini et al. (2020)'},
 'Arraiano2007TAG': {'citation': 'Arraiano LS, Kirby J, Brown JKM (2007) Cytogenetic analysis of the susceptibility of the wheat line Hobbit sib (Dwarf A) to Septoria tritici blotch. Theor Appl Genet 116:113-122',
                     'url': 'https://doi.org/10.1007/s00122-007-0651-9', 'short': 'Arraiano et al. (2007)'},
 'Agostinelli2012': {'citation': 'Agostinelli AM, Clark AJ, Brown-Guedira G, Van Sanford DA (2012) Optimizing phenotypic and genotypic selection for Fusarium head blight resistance in wheat. Euphytica 186:115-126',
                     'url': 'https://doi.org/10.1007/s10681-011-0499-6', 'short': 'Agostinelli et al. (2012)'},
}
SNN1, RMG7, RMG8, FHB1, STB6, TRP = 2.356, 761.8578, 796.0231, 8.528, 26.2035, 22.528
ZAV = ' Position in the paper is on the wild emmer Zavitan v2.0 assembly; no Chinese Spring anchor is given, so the Zavitan position is used directly as an approximate RefSeq v2.1 position (grade C).'
BA = 'Locus from GWAS of the durum EPO panel (composite cross-population, TaBW280K array, 168,725 SNPs), seedling tests in the glasshouse.'
B = [  # id, name, disease, chr, zav Mb, size Mb, stats, placement(peak, grade, text)
 ('SNB_Tox1B', 'Septoria nodorum blotch', '1B', 1.03, 1.8, 'Tox1 infiltration; -log10 P 9.2, R2 0.22', (SNN1, 'B', ' The peak is 66 kb from Snn1 and the paper concludes it is Snn1; placed at the Snn1 position in the atlas (Chinese Spring).')),
 ('SNB_Tox7A', 'Septoria nodorum blotch', '7A', 48.3, 1.9, 'Tox1 infiltration; -log10 P 8.0, R2 0.19', None),
 ('SNB_Nec4B', 'Septoria nodorum blotch', '4B', 641.3, 1.8, 'spore inoculation (necrosis %); -log10 P 6.63, R2 0.18', None),
 ('SNB_Nec7B', 'Septoria nodorum blotch', '7B', 684.2, 2.2, 'spore inoculation (necrosis %); -log10 P 4.73, R2 0.09', None),
 ('WB_BR43_2A', 'Wheat blast', '2A', 771.0, 0.5, 'isolate BR43; -log10 P 11.2, R2 0.23', (RMG7, 'C', ' Co-localises with the 2 Mb Rmg7 locus (Zavitan 771.35-773.37 Mb); placed at the Rmg7 position in the atlas (Chinese Spring).')),
 ('WB_BR43_2B', 'Wheat blast', '2B', 795.0, 5.0, 'isolate BR43; -log10 P 8.8, R2 0.19', (RMG8, 'C', ' Maps close to the region of Rmg8; placed at the Rmg8 position in the atlas (Chinese Spring).')),
 ('WB_BR43_5B', 'Wheat blast', '5B', 711.0, 0.1, 'isolate BR43; -log10 P 6.7, R2 0.16', None),
 ('WB_BR32_7A', 'Wheat blast', '7A', 722.0, 17.0, 'isolate BR32; -log10 P 6.3, R2 0.18; top SNP 75 kb from NLR TRIDC7AG077360', None),
 ('WB_BR32_7B', 'Wheat blast', '7B', 744.0, 24.0, 'isolate BR32; -log10 P 5.1, R2 0.14', None),
 ('STB_Nec3A', 'Septoria tritici blotch', '3A', 25.1, 0.3, 'isolate P1a (necrosis %); -log10 P 13.0, R2 0.45; 327 kb region with six WAK genes, one homologous to Stb6', (TRP + 3.3, 'C', ' 3.3 Mb from STB_Pyc3A in the paper; placed 3.3 Mb distal of the Chinese Spring anchor of STB_Pyc3A (TraesCS3A01G041800).')),
 ('STB_Pyc3A', 'Septoria tritici blotch', '3A', 21.8, 0.3, 'isolate P1a (pycnidia %); -log10 P 6.7, R2 0.11; complete inhibition of pycnidia confirmed with isolates P1a and M1a in elite lines Coussur, Miradoux, Nefer and Plussur', (TRP, 'B', ' The 326 kb region contains the ortholog of Chinese Spring TraesCS3A01G041800 (tryptophan synthase); placed at that gene (22.53 Mb).')),
 ('STB_Nec5A', 'Septoria tritici blotch', '5A', 649.3, 4.7, 'isolate P1a (necrosis %); -log10 P 5.0, R2 0.11', None),
 ('LR_2B', 'Leaf rust', '2B', 773.2, 0.18, 'isolate AD2015; mock -log10 P 16.47 (R2 0.22), chitin 17.33 (0.19)', None),
 ('LR_6A', 'Leaf rust', '6A', 0.89, 0.01, 'isolate AD2015; mock -log10 P 9.11 (R2 0.21), chitin 7.92 (0.17)', None),
 ('LR_6B11', 'Leaf rust', '6B', 12.3, 3.24, 'isolate AD2015; mock -log10 P 7.36 (R2 0.18), chitin 5.39 (0.13); strongest LR QTL, NLR and LRR-kinase clusters; relation to Lr36 and Lr61 unresolved', None),
 ('LR_6B135', 'Leaf rust', '6B', 135.2, 2.4, 'isolate AD2015; mock -log10 P 8.78 (R2 0.14), chitin 7.82 (0.17)', None),
]
def mk(i, **k):
    r = dict(id=i, cls='QTL', ref='RefSeq v1.0', start=None, end=None, peak=None, marker='', pve='', stat='', curation='hand'); r.update(k); return r
def build(reports):
    old = ('S100-', 'S101-', 'S102-', 'S103-')
    keep = [r for r in reports if not r['id'].startswith(old)]; new = []
    for n, (nm, dis, ch, z, sz, st, pl) in enumerate(B, 1):
        r = mk(f'S100-{n:03d}', name=nm + ' (EPO durum panel)', cls='MTA', disease=dis, chr=ch, source='Ballini2020', stat=st)
        if pl: r['peak'], r['grade'] = round(pl[0], 4), pl[1]; tail = pl[2] + f' Zavitan v2.0 peak {z} Mb, support interval {sz} Mb.'
        else:
            r['start'], r['end'], r['peak'], r['grade'], r['ref'] = round(max(0, z - sz / 2), 3), round(z + sz / 2, 3), z, 'C', 'RefSeq v2.1'; tail = ZAV + f' Zavitan v2.0 peak {z} Mb, support interval {sz} Mb.'
        r['note'] = BA + f' {st}.' + tail + ' Read from the paper (Table 1), Nov 2026 (ad23.py).'; r['evidence'] = f'Table 1: {nm} {ch} {z} Mb (Zavitan v2.0)'; new.append(r)
    p1, h1 = resolve('Xgwm540', '5B'); p2, h2 = resolve('Xgwm67', '5B')
    new.append(mk('S101-001', name='Septoria susceptibility, Hobbit sib 5BS', disease='Septoria tritici blotch', chr='5B', source='Arraiano2007TAG', start=p1, end=p2, grade='C',
        marker='Xwmc386-Xgwm67', note='Hobbit sib (Dwarf A) chromosome arm 5BS carries genes that promote susceptibility to septoria or suppress resistance: the substitution line in which the proximal '
        'part of Hobbit sib 5BS (Xwmc386-Xgwm67, 29-49 cM on the Somers consensus map) had been replaced by Chinese Spring 5BS was much less diseased; confirmed with substitution and nullitetrasomic lines. '
        f'Placed between Xgwm540 ({p1:.3f} Mb, {h1}) and Xgwm67 ({p2:.3f} Mb, {h2}): a chromosome-arm region, not a locus. Read from the paper (Table 1, Results), Nov 2026 (ad23.py).',
        evidence='Table 1: SSR alleles of Hs (Bez 5BS-7BS) match Chinese Spring from Xwmc386 to Xgwm67 on 5BS'))
    new.append(mk('S101-002', name='Stb6 (Bezostaya 1)', cls='QTL', disease='Septoria tritici blotch', chr='3A', source='Arraiano2007TAG', peak=STB6, grade='C', marker='Xgwm369',
        note='Bezostaya 1 had isolate-specific resistance to M. graminicola IPO323, attributed to Stb6 on 3AS (CD (Bez3A) substitution line resistant; Xgwm369 allele shared with Flame, the Stb6 source). '
        'Placed at the Stb6 position in the atlas. Bezostaya 1 also had partial resistance spread over several chromosomes (not mapped). Read from the paper, Nov 2026 (ad23.py).',
        evidence='Results: IPO323 resistance of CD (Bez3A); Xgwm369 allele as Flame'))
    pc, hc = resolve('Xcfd233', '2D')
    AG = 'Effect in F2-derived lines of KY93C-1238-17-2 x VA01W-476 (soft red winter wheat; scab nurseries at Lexington 2007-2008 and Princeton 2008, KY); resistance allele from VA01W-476.'
    new.append(mk('S102-001', name='Fhb1 (VA01W-476)', disease='Fusarium head blight', chr='3B', source='Agostinelli2012', start=FHB1, end=FHB1, grade='C', marker='Xsts3B-256, Xgwm533, Xbarc147',
        note=AG + ' Fhb1 alone reduced DON by 25.4% and FDK by 31.7%; with QFhs.nau-2DL, 69.8% and 63.5%. Stable across environments. Placed at the cloned Fhb1 (TaHRC) position in the atlas. Read from the paper (Table 3), Nov 2026 (ad23.py).',
        evidence='Table 3: RR+RS (Fhb1) DON effect 25.4%, FDK 31.7%'))
    new.append(mk('S102-002', name='QFhs.nau-2DL (VA01W-476)', disease='Fusarium head blight', chr='2D', source='Agostinelli2012', peak=pc, grade='C', marker='Xcfd233, Xgwm608',
        note=AG + ' QFhs.nau-2DL alone reduced DON by 54.9% and FDK by 40.4%, more than Fhb1 in this background; significant QTL x environment interaction. '
        f'Only Xcfd233 could be placed ({pc:.3f} Mb, {hc}); position is that marker. Read from the paper (Table 3), Nov 2026 (ad23.py).', evidence='Table 3: RR+SR (QFhs.nau-2DL) DON effect 54.9%, FDK 40.4%'))
    AR = ('MTA from GWAS of 273 soft red winter wheat lines (19,992 GBS SNPs, compressed MLM), as reprinted in Arruda et al. 2016 Mol Breed 36:84, Table 2 (significant at FDR 0.10). ')
    A = [('Fhb1 region (3B, GBS SNPs in LD with umn10/gwm533)', '3B', 'IWGSC_CSS_3B_scaff_10676713_7175, _10352272_5482, _10699215_3620', 'SEV, FHB index, ISK and DON; 3.4-18.3 cM; r2 up to 0.08', FHB1,
          ' The SNPs are in LD with the Fhb1 markers umn10 and gwm533; placed at the cloned Fhb1 position in the atlas.'),
         ('MTA 3B 6.9 cM (severity)', '3B', 'IWGSC_CSS_3B_scaff_10698462_2332', 'SEV; 6.86 cM; r2 0.04', None, ''),
         ('MTA 3B 73.7 cM (DON)', '3B', 'IWGSC_CSS_3B_scaff_10413672_4839, IWGSC_CSS_3B_scaff_10764618_2168', 'DON; 73.67 cM; r2 0.06-0.07', None, ''),
         ('MTA 7D 70.8 cM (incidence)', '7D', 'IWGSC_CSS_7DS_scaff_3876750_2023', 'INC and ISK; 70.84 cM; r2 0.16', None, ''),
         ('MTA 6A 134.2 cM (incidence)', '6A', 'IWGSC_CSS_6AL_scaff_5780077_12152', 'INC; 134.15 cM; r2 0.06', None, ''),
         ('MTA 4D 0 cM (incidence)', '4D', 'IWGSC_CSS_4DS_scaff_2300354_4482', 'INC; 0 cM; r2 0.06', None, ''),
         ('MTA 4A 78.4 cM (incidence)', '4A', 'IWGSC_CSS_4AL_scaff_7146617_11335', 'INC; 78.35 cM; r2 0.06', None, ''),
         ('MTA 7A 22.8 cM (incidence)', '7A', 'IWGSC_CSS_7AS_scaff_4132011_1400', 'INC; 22.82 cM; r2 0.05', None, ''),
         ('MTA 1D 19.0 cM (DON)', '1D', 'IWGSC_CSS_1DS_scaff_1879930_3352', 'DON; 19.04 cM; r2 0.06', None, '')]
    for n, (nm, ch, mk_, st, p, ex) in enumerate(A, 1):
        r = mk(f'S103-{n:03d}', name=nm + ' (Arruda 2016)', cls='MTA', disease='Fusarium head blight', chr=ch, source='A_FHB_024', marker=mk_, stat=st, peak=p, grade='C' if p else 'D')
        r['note'] = AR + st + '.' + (ex if p else ' The SNP is named by its IWGSC chromosome survey sequence scaffold and only a cM position is given, so it cannot be placed.') + ' Read from the paper, Nov 2026 (ad23.py).'
        r['evidence'] = 'Arruda et al. 2016 Mol Breed Table 2: ' + mk_ + '; ' + st; new.append(r)
    for r in new:
        r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}
    return keep + new, new, []
if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json')); M['reports'], new, _ = build(M['reports']); M['sources'].update(SRC); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'], new2, _ = build(Sd['reports'])
    if 'sources' in Sd: Sd['sources'].update(SRC)
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    ops = [{'op': 'set', 'path': 'sources/' + k, 'data': v} for k, v in SRC.items()]
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + r['id'], 'data': d})
    json.dump(ops, open('/tmp/ops23.json', 'w'), ensure_ascii=False)
    for o in ops[3:]: d = o['data']; print(d['id'], d['chr'], d.get('start'), d.get('end'), d.get('peak'), d['grade'], d['asm'].get('RefSeq v2.1'))
