"""ad57 (Nov 2026): journal versions uploaded for papers whose preprints were removed in ad55, plus one new paper.
- Peters Haugrud et al. (2023) Mol Breed 43:54 (A_SNB_001): Table 4, ten loci (RefSeq v2 bp), plus the 5BS SnTox3 SNP pair at 65.76 Mb.
  R0633 (Tsc1) corrected: tan spot, not SNB.
- Gardner et al. (2025) TAG 138:133 (Z_1123): Table 2, 34 replicated loci (brown rust 12, powdery mildew 2, yellow rust 20; RefSeq v1.0).
  The 24 existing records were all labelled Septoria tritici blotch (no STB locus replicated); corrected in place, ten loci added.
- Liu et al. (2021) TAG 134:3279 (V_Liu2021b): QBp.caas-3BL fine-mapped to 1.7 Mb between TraesCS3B02G404200 (M8) and TraesCS3B02G404800 (M10).
  V1442/V1443 removed: they restate QTL of another study cited in the introduction.
- Winn et al. (2022) TAG 135:3177 (Z_1463): Table 1 KASP positions (RefSeq v1.0) and Table 3 SUWWSN effects. Z0137 had accuracy values (0.96-0.98)
  read as Mb; corrected. V0353/V0354 removed: positions came from other papers; Qfhb.vt-1B is covered by Z0137.
- Chaudhary et al. (2024) G3 14:jkae240 (new, CHAUD2024G3): fusarium foot rot (F. graminearum) disease-score MTAs, placed in silico."""
import json, re, sys, gzip
exec(open('scripts/curation/lib/insilico_res.py').read())
sys.path.insert(0, 'scripts/curation/sb19'); import lift
T = ' Read from the journal PDF, Nov 2026 (ad57.py).'
SRC = {'CHAUD2024G3': {'citation': 'Chaudhary S, Ricardo RMN, Dubey M, Jensen DF, Grenville-Briggs L, Karlsson M (2024) Genotypic variation in winter wheat for fusarium foot rot and its biocontrol using Clonostachys rosea. G3 Genes|Genomes|Genetics 14:jkae240',
                       'url': 'https://doi.org/10.1093/g3journal/jkae240', 'short': 'Chaudhary et al. (2024)'}}
UPD, NEW, DEL = [], [], ['V1442', 'V1443', 'V0353', 'V0354']
def rec(i, **k):
    r = dict(id=i, cls='MTA', curation='hand'); r.update(k); return r
# ---- Peters Haugrud 2023 (RefSeq v2 bp; positions printed -> A)
PH = 'Global winter wheat panel (264 lines), 90K SNPs, GWAS (MLM/BLINK, FDR-adjusted). '
ph = [('S200-001', 'Tsn1 (SnToxA sensitivity; SNB Sn2000, NOR4, SnIr05H71a)', 'Septoria nodorum blotch', '5B', 549.788309, 549.933008, 'IACX9261, Excalibur_rep_c105815_305, tplb0027f13_452, Tdurum_contig25513_123', '14 SNPs for SnToxA; also the top SNPs for disease with isolates Sn2000, SnIr05H71a (tplb0027f13_452, 549,914,572 bp) and NOR4 (Tdurum_contig25513_123).'),
      ('S200-002', 'Snn1 (SnTox1 sensitivity)', 'Septoria nodorum blotch', '1B', 2.358687, 2.381165, 'BS00064465_51, RAC875_rep_c97766_145', 'Snn1 gene at 2,383,732 bp (RefSeq v2.0). Table 4 prints the end as "2,381-165".'),
      ('S200-003', 'Snn3-B1 (SnTox3 sensitivity; SNB Sn4)', 'Septoria nodorum blotch', '5B', 6.704549, 7.091681, 'GENE-3324_338, Kukri_c912_1423, BS00064298_51', '11 SNPs for SnTox3; the four Sn4 SNPs span 6,704,549-7,031,414 bp. Snn3-D1 BLAST hit on 5B at 7,041,691 bp.'),
      ('S200-004', 'Snn5 (SnTox5 sensitivity)', 'Septoria nodorum blotch', '4B', None, None, 'BS00023035_51', 'BS00023035_51 at 635,843,775 bp, within the Snn5 region.', 635.843775),
      ('S200-005', 'Snn7 region (SNB AR2-1)', 'Septoria nodorum blotch', '2D', None, None, 'IAAV5315', 'Table 4 gives 652,955,221 bp; the text gives 652,955,177 bp for IAAV5315.', 652.955221),
      ('S200-006', 'SnTox5 sensitivity 5A (RAC875_c44214_322)', 'Septoria nodorum blotch', '5A', None, None, 'RAC875_c44214_322', 'Within TraesCS5A03G0882600 (development and cell death domain protein); possibly a novel SnTox5 sensitivity locus.', 569.482195),
      ('S200-007', 'SnTox5 sensitivity 5B (BS00002208_51)', 'Septoria nodorum blotch', '5B', None, None, 'BS00002208_51', 'Within TraesCS5B03G0555300; possibly a novel SnTox5 sensitivity locus.', 381.83262),
      ('S200-008', 'SNB 2AL (Sn2000, SnIr05H71a)', 'Septoria nodorum blotch', '2A', 747.578919, 784.753007, 'Excalibur_rep_c105106_665, wsnp_Ex_c63909_62932437', 'Two single SNPs, one per isolate: Excalibur_rep_c105106_665 (747,578,919 bp, SnIr05H71a; beta-glucosidase TraesCS2A03G1210000) and wsnp_Ex_c63909_62932437 (784,753,007 bp, Sn2000; TraesCS2A03G1360000). Table 4 lists them as one locus; they lie 37 Mb apart.'),
      ('S200-009', 'Tan spot 7DS (DW5, race 5)', 'Tan spot', '7D', None, None, 'Excalibur_c34115_226', 'Single SNP; predicted disease resistance genes within 2 Mb.', 1.306537),
      ('S200-010', 'SnTox3-associated 5BS proximal (wsnp_Ku_c439_913400)', 'Septoria nodorum blotch', '5B', 65.759133, 65.759716, 'wsnp_Ku_c439_913400, Ex_c1846_1818', 'Two further SnTox3 SNPs proximal to Snn3-B1, within TraesCS5B03G0147600 (PK-NLR).')]
for x in ph:
    i, nm, dis, ch, s, e, mk, nt = x[:8]; pk = x[8] if len(x) > 8 else None
    NEW.append(rec(i, name=nm, disease=dis, chr=ch, start=s, end=e, peak=pk, ref='RefSeq v2.1', marker=mk, source='A_SNB_001', grade='A',
                   evidence='Table 4 and Results', note=PH + nt + ' Positions as printed (IWGSC RefSeq v2).' + T))
UPD.append(rec('R0633', name='Tsc1 region (tan spot Pti2, AR CrossB10)', disease='Tan spot', chr='1A', start=5.048054, end=8.249438, peak=None, ref='RefSeq v2.1',
               marker='Kukri_c22349_754, BobWhite_c20902_177', source='A_SNB_001', grade='A', evidence='Table 4: Tsc1 1A 5,048,054-8,249,438',
               note=PH + 'Kukri_c22349_754 (8,249,438 bp) for Pti2 (race 1) and BobWhite_c20902_177 (5,048,054 bp) for AR CrossB10. Earlier record wrongly labelled SNB; corrected.' + T))
# ---- Gardner 2025 Table 2 (RefSeq v1.0)
G = '''BR_1B117||1B|1.20|11.31|BS00009715_51|112.4|6.44|-0.50|5|V0852
BR_2A015|MT25Mb_1|2A|2.23|24.28|BS00004089_51|14.8|13.95|-0.67|2|V0853
BR_2A704||2A|703.98|717.17|Excalibur_c40617_983|704.0|3.58|0.83|2|Z0373
BR_2B026|MT25Mb_2|2B|26.30|33.84|Kukri_c40764_367|26.3|5.80|-0.48|2|Z0374
BR_2B640||2B|214.59|640.79|BS00069685_51|640.8|3.83|-0.36|2|V0854
BR_2B777|MT25Mb_3|2B|||RAC875_c19685_944|777.3|3.57|0.47|2|V0855
BR_2D014|MT25Mb_4|2D|||BobWhite_c15073_502|14.4|3.34|-0.27|2|V0856
BR_3A134||3A|134.25|471.36|BobWhite_c35303_192|134.2|4.04|-0.44|2|V0857
BR_3A741|MT25Mb_5|3A|736.39|741.24|Tdurum_contig5009_349b|741.2|3.61|-0.28|2|Z0375
BR_3B826||3B|810.33|826.1|Tdurum_contig42131_1300|826.1|3.47|1.18|2|Z0376
BR_4A714|MT25Mb_6|4A|698.04|743.98|BobWhite_c20306_111|713.5|16.27|-0.75|2|Z0377
BR_6A016||6A|11.41|23.72|Excalibur_rep_c105463_330|15.7|4.30|-0.38|3|Z0378
PM_4A734|MT25Mb_6|4A|610.49|742.09|BS00110758_51|734.0|5.22|-0.62|2|S201-001
PM_6B664||6B|663.68|681.92|RAC875_c5129_280|663.7|5.04|0.68|2|Z0379
YR_1B164||1B|14.85|163.78|Excalibur_rep_c92475_275|163.8|4.82|0.26|2|Z0380
YR_2A010|MT25Mb_1|2A|0.40|24.28|wsnp_Ra_c8771_14786376|9.5|12.48|-0.77|5|S201-002
YR_2B034|MT25Mb_2|2B|24.91|33.84|BobWhite_rep_c64429_660|33.8|10.61|1.96|4|Z0381
YR_2B047||2B|||BS00041587_51|47.4|5.95|0.81|6|S201-003
YR_2B155||2B|133.7|154.99|Kukri_c36783_91|155.0|4.17|-0.57|2|Z0382
YR_2B763|MT25Mb_3|2B|0.00|763.09|BS00070301_51|763.1|4.54|-0.30|2|S201-004
YR_2D014|MT25Mb_4|2D|0.00|13.99|RAC875_c90426_151|14.0|5.96|1.24|3|S201-005
YR_3A008||3A|7.43|8.87|BS00037189_51|7.8|4.86|0.73|2|S201-006
YR_3A746|MT25Mb_5|3A|||wsnp_Ex_c60462_60905848|746.2|5.69|-0.68|2|S201-007
YR_3B739||3B|738.75|750.36|wsnp_Ex_rep_c101457_86817938|739.1|4.85|-0.48|2|Z0383
YR_4A738|MT25Mb_6|4A|734.00|738.78|Excalibur_c65272_341|737.5|6.70|-0.62|8|Z0384
YR_5A030||5A|15.85|29.51|Excalibur_rep_c90275_262|29.5|6.07|0.77|2|Z0385
YR_5A560||5A|559.52|568.27|BS00021955_51|559.5|4.21|-0.60|2|Z0386
YR_5A677||5A|677.13|684.94|BS00022867_51|677.1|5.07|-0.46|4|Z0387
YR_6A002||6A|0.29|13.07|Excalibur_c50323_215|1.9|5.83|-0.44|5|S201-008
YR_6A610||6A|596.58|611.13|GENE-4021_496|610.3|12.91|0.54|6|S201-009
YR_7A090||7A|89.84|90.65|Kukri_c5757_530|90.7|4.17|1.23|2|Z0388
YR_7A626||7A|625.74|632.59|BS00105558_51|625.7|5.73|-0.48|3|Z0389
YR_7A730||7A|||Kukri_c11451_1882|730.4|4.97|-0.63|2|S201-010
YR_7B701||7B|700.57|705.73|Excalibur_c7338_563|700.8|6.04|0.51|3|Z0390'''
GD = {'BR': 'Leaf rust', 'PM': 'Powdery mildew', 'YR': 'Stripe rust'}
GB = 'North-western European winter wheat panel (480 cultivars), 90K SNP array, 31 field trials in the UK, Germany, Denmark, France and Sweden; loci significant in two or more trials. '
FIX = {'BR_1B117': ' Table 2 prints the interval as 1.20-11.31 Mbp while the peak SNP is at 112.4 Mbp; the interval is likely misprinted, so only the peak is used.',
       'YR_3A746': ' Table 2 prints the peak as 7746.2 Mbp; read as 746.2 Mbp (interval 746.17).',
       'YR_6A610': ' Table 2 prints the interval as 596.58-11.13 Mbp; read as 596.58-611.13 Mbp (peak 610.3).',
       'YR_2A010': ' Resistance from the Ae. ventricosa 2NvS introgression (Yr17/Lr37 segment) carries effective yellow and brown rust resistance beyond Yr17 and Lr37; validated in bi-parental populations.',
       'BR_2A015': ' Same 2NvS (Ae. ventricosa) segment as YR_2A010.', 'YR_6A610': ' Table 2 prints the interval as 596.58-11.13 Mbp; read as 596.58-611.13 Mbp (peak 610.3). Validated in bi-parental populations.'}
for l in G.splitlines():
    loc, mt, ch, s, e, mk, pk, lp, eff, nt, i = l.split('|')
    s = float(s) if s else None; e = float(e) if e else None; pk = float(pk)
    if loc == 'BR_1B117': s = e = None
    if s is not None and s == e: e = None
    wide = s is not None and e - s > 100
    note = GB + f'Peak SNP {mk} at {pk} Mbp; effect of the reference allele {eff.replace("-", "minus ")} % leaf infection; significant in {nt} trials.' + (f' Part of multi-resistance locus {mt}.' if mt else '') + FIX.get(loc, '') + \
        (f' Interval width {e - s:.0f} Mb: a chromosome region as the source reports it, not a locus.' if wide else '') + ' Positions as printed (IWGSC RefSeq v1.0).' + T
    r = rec(i, name=f'{loc} ({mk})', disease=GD[loc[:2]], chr=ch, start=s, end=e, peak=pk, ref='RefSeq v1.0', marker=mk, stat=f'-log10P {lp}',
            source='Z_1123', grade='A', evidence=f'Table 2: {loc} {ch} ' + (f'{s}-{e}' if s is not None and e else '') + f' peak {pk}', note=note, pve='')
    (NEW if i.startswith('S2') else UPD).append(r)
# ---- Liu 2021 QBp.caas-3BL (gene models, RefSeq v1.0/v1.1)
UPD.append(rec('S202-001', cls='QTL', name='QBp.caas-3BL', disease='Spot blotch', chr='3B', start=638.960471, end=640.068645, peak=None, ref='RefSeq v1.0',
               marker='M8 (TraesCS3B02G404200), M10 (TraesCS3B02G404800); KASP from TraesCS3B02G404700', pve='6.0-12.2', source='V_Liu2021b', grade='A',
               evidence='Fine mapping between M8 (TraesCS3B02G404200) and M10 (TraesCS3B02G404800), 1.7 Mb',
               note='Black point resistance (the paper notes Bipolaris sorokiniana, Alternaria alternata and Fusarium proliferatum among causes). Linmai 2 x Zhong 892 RILs, seven environments; four NIL sets narrowed the QTL to 1.7 Mb with five genes (TraesCS3B02G404300 to 404700). Candidates: TraesCS3B02G404300, TraesCS3B02G404600 (ABI5 homologues) and TraesCS3B02G404700 (protein phosphatase; nonsynonymous SNP converted to a KASP marker, validated in 166 cultivars). Interval from the RefSeq v1.0 coordinates of the flanking gene models.' + T))
NEW.append(UPD.pop())
# ---- Winn 2022
W = 'SunGrains breeding panels genotyped with KASP assays and GBS; haplotype calls predicted by machine learning; effects estimated in the 2011-2020 Southern Uniform Winter Wheat Scab Nursery (Table 3). '
w = [('S203-001', 'Qfhb.nc-1A', '1A', 541.9, 569.75, 'IWA3805, IWA886', 'From NC-Neuse. Resistant haplotype: DON -1.23 ppm (P 0.0011), FDK -2.56% (P 0.016), severity -2.37% (P 0.014).', 'DON -1.23 ppm, P = 0.0011'),
     ('S203-002', 'Qfhb.nc-4A', '4A', 543.87, 575.63, 'IWA2900, IWA402, IWA2793', 'From NC-Neuse. No significant effect in the scab nursery (DON +0.65 ppm, P 0.15; severity -0.17%, P 0.49).', 'not significant (P > 0.14)'),
     ('Z0137', 'Qfhb.vt-1B (Jamestown haplotype, QFHB.vt-1B.1 and QFHB.vt-1B.2)', '1B', 336.46, 477.41, 'IWB43992, IWA6259, IWA7594', 'From Jamestown. Resistant haplotype: DON -1.66 ppm (P 0.0001), severity -3.20% (P 0.020). Earlier record held accuracy values (0.96-0.98) read as Mb; corrected. Interval width 141 Mb: a chromosome region spanned by the haplotyping markers, not a locus.', 'DON -1.66 ppm, P = 0.0001'),
     ('S203-003', 'Fhb1 (KASP haplotype)', '3B', 13.64, 13.96, 'TaHRC, snp3BS-8', 'Resistant haplotype: DON -2.19 ppm, FDK -6.72%, severity -8.82% (all P < 0.001); also delays heading by 0.93 d. Positions as printed in Table 1 (TaHRC at 13.64 Mbp); other atlas records place TaHRC near 8.5 Mb.', 'severity -8.82%, P < 0.001')]
for i, nm, ch, s, e, mk, nt, stt in w:
    r = rec(i, cls='QTL', name=nm, disease='Fusarium head blight', chr=ch, start=s, end=e, peak=None, ref='RefSeq v1.0', marker=mk, stat=stt, source='Z_1463', grade='A' if e - s <= 100 else 'C',
            evidence='Table 1 KASP positions (Mbp); Table 3', note=W + nt + ' Interval spans the KASP markers (RefSeq v1.0, Table 1).' + T, pve='')
    (UPD if i == 'Z0137' else NEW).append(r)
# ---- Chaudhary 2024 G3 (in silico)
C = 'Fusarium foot and root rot (F. graminearum) seedling assay, 190 northern European winter wheat genotypes, 20K SNP array (7,360 SNPs), GWAS on disease score in treatment Fg (P <= 0.00014). '
c = [('S204-001', 'Foot rot 1A (53 cM)', '1A', ['Excalibur_c7026_2635'], 'No difference in disease score between alleles.'),
     ('S204-002', 'Foot rot 2A (115-116 cM)', '2A', ['BS00089497_51', 'Kukri_c40121_373'], 'Minor alleles GG (BS00089497_51) and AA (Kukri_c40121_373) reduce disease score (P < 0.05).'),
     ('S204-003', 'Foot rot 4B (71-73 cM)', '4B', ['BS00096604_51', 'RFL_Contig2459_2314', 'Ku_c33858_325'], 'Minor alleles TT, GG and CC reduce disease score (P < 0.05). The candidate-gene paragraph puts BS00096604_51 on 4D; the association results put all three on 4B.')]
for i, nm, ch, ms, nt in c:
    ok = [(m,) + resolve(m, ch) for m in ms]; ok = [x for x in ok if x[1]]
    r = rec(i, name=nm, disease='Fusarium crown rot', chr=ch, ref='RefSeq v1.0', marker=', '.join(ms), source='CHAUD2024G3', evidence='Results: genome-wide marker-trait associations', pve='')
    if len(ok) >= 2:
        lo, hi = min(x[1] for x in ok), max(x[1] for x in ok); r.update(start=lo, end=hi, grade='B' if hi - lo <= 100 else 'C')
    elif ok: r.update(peak=ok[0][1], grade='C')
    else: r['grade'] = 'D'
    pl = ('; '.join(f'{m} {p:.3f} Mb ({h})' for m, p, h in ok) + '.') if ok else 'no marker placed.'
    r['note'] = C + nt + ' Placed in silico (RefSeq v1.0): ' + pl + (f' Span {hi-lo:.0f} Mb: a chromosome region, not a locus.' if len(ok) >= 2 and hi - lo > 100 else '') + T
    NEW.append(r)
for r in UPD + NEW: r['asm'] = lift.lift_record(r) if r['grade'] != 'D' else {'_note': 'no physical position in source'}
def apply(rs):
    U = {r['id']: r for r in UPD}; out = []
    for r in rs:
        if r['id'] in DEL: continue
        if r['id'] in U:
            n = dict(r); n.update(U[r['id']])
            for k in ('snpseq', 'primers'): 
                if k in r: n[k] = r[k]
            out.append(n)
        else: out.append(r)
    have = {r['id'] for r in out}
    return out + [r for r in NEW if r['id'] not in have]
if __name__ == '__main__':
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    M = json.load(open('data/merged_v13f.json')); M['reports'] = apply(M['reports']); M['sources'].update(SRC); json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read(); m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</')); Sd['reports'] = apply(Sd['reports']); Sd['sources'].update(SRC)
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    R = {r['id']: r for r in M['reports']}
    ops = [{'op': 'set', 'path': 'sources/' + k, 'data': v} for k, v in SRC.items()]
    for i in [r['id'] for r in UPD + NEW]:
        r = R[i]; d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}; d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        ops.append({'op': 'set', 'path': 'reports/' + i, 'data': d})
    for i in DEL: ops.append({'op': 'set', 'path': 'reports/' + i, 'data': {'deleted': True, 'deleted_by': 'ad57', 'deleted_at': '2026-11-17', 'name': i}})
    json.dump(ops, open('/tmp/ops57.json', 'w'), ensure_ascii=False, separators=(',', ':'))
    print(len(UPD), 'updated', len(NEW), 'new', len(DEL), 'deleted', len(ops), 'ops', len(json.dumps(ops)))
    for r in UPD + NEW: print(r['id'], r['disease'][:14], r['chr'], r.get('start'), r.get('end'), r.get('peak'), r['grade'], (r['asm'].get('RefSeq v2.1') or r['asm'].get('RefSeq v1.0') or {}).get('flag'))
