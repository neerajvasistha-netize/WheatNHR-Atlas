"""Batch 6 (Nov 2026): spot blotch / black point / leaf spot GWAS papers (GWAS_SB_1.zip, GWAS_BPR.zip, GWAS_LS.zip).
- new papers: Bainsla 2020, Aggarwal 2026, Lv 2020, Tang 2022, Iqbal 2022 (all positions printed in the papers, grade A)
- Chen 2025 (A_spot_blotch_033): 5 auto-extracted records replaced by the 6 QTL of Table 2 (printed intervals)
- Juliana 2022: marker 5B_68806 was a line-wrap truncation of 5B_688063050; 5B_687326022 (named in the text) added
- regrades only (positions unchanged): see REGRADE
- Thirunarayanan 2024 citation re-keyed (short citation read "P et al.")
Writes data/merged_v13f.json, web/index.html (seed) and /home/claude/work/gw60_db.json."""
import json, re, sys
sys.path.insert(0, 'scripts/curation/sb19')
import lift

BP, TS, CB, LR, LS, YR = 'Spot blotch', 'Tan spot', 'Common bunt', 'Leaf rust', 'Leaf spot complex', 'Stripe rust'
def S(cit, url, short): return dict(citation=cit, url=url, short=short)
NEW_SOURCES = {
 'Bainsla2020PP': S('Bainsla NK, Phuke RM, He X, Gupta V, Bishnoi SK, Sharma RK, Ataei N, Dreisigacker S, Juliana P, Singh PK (2020) Genome-wide association study for spot blotch resistance in Afghan wheat germplasm. Plant Pathology 69:1161-1171', 'https://doi.org/10.1111/ppa.13191', 'Bainsla et al. (2020)'),
 'Aggarwal2026TPG': S('Aggarwal N, He X, Rathore M, Jan F, Gupta V, Mir RR, Singh PK (2026) Identification of resistance sources and genomic regions regulating spot blotch resistance in Asian bread wheat (Triticum aestivum L.) via genome-wide association study. The Plant Genome 19:e70228', 'https://doi.org/10.1002/tpg2.70228', 'Aggarwal et al. (2026)'),
 'Lv2020PD': S('Lv G, Dong Z, Wang Y, Geng J, Li J, Lv X, Sun C, Ren Y, Zhang J, Chen F (2020) Identification of genetic loci of black point in Chinese common wheat by genome-wide association study and linkage mapping. Plant Disease 104:2005-2013', 'https://doi.org/10.1094/PDIS-12-19-2733-RE', 'Lv et al. (2020)'),
 'Tang2022PD': S('Tang H, Tan Z, Wang X, Yang L, Chen G, Yu H, Pu Z, Jiang Q, Li M, Cheng M, Qi P, Li W, Liu Y, Wang J (2022) Genome-wide association study of kernel black point resistance in Chinese wheat landraces. Plant Disease 106:1428-1433', 'https://doi.org/10.1094/PDIS-09-21-1898-RE', 'Tang et al. (2022)'),
 'Iqbal2022Plants': S('Iqbal M, Semagn K, Jarquin D, Randhawa H, McCallum BD, Howard R, Aboukhaddour R, Ciechanowska I, Strenzke K, Crossa J, Ceron-Rojas JJ, N\'Diaye A, Pozniak C, Spaner D (2022) Identification of disease resistance parents and genome-wide association mapping of resistance in spring wheat. Plants 11:2905', 'https://doi.org/10.3390/plants11212905', 'Iqbal et al. (2022)'),
 'Thirunarayanan2024PMPP': S('Thirunarayanan P, Korra T, Chandra R (2024) Genome-wide association study of lignin trait in elite spring wheat against spot blotch disease. Physiological and Molecular Plant Pathology 133:102388', 'https://doi.org/10.1016/j.pmpp.2024.102388', 'Thirunarayanan et al. (2024)'),
}
RESOURCE = {'A_spot_blotch_022': 'Thirunarayanan2024PMPP'}
REMOVE = ['R0627', 'R0628', 'R0630', 'R0631', 'V0482', 'R0413']
REGRADE = {}
G = json.load(open('/home/claude/work/b6_regradeA.json'))
for k in ('A_spot_blotch_022', 'A_spot_blotch_021', 'A_spot_blotch_018', 'Liu2017BMC'):
    for i in G[k]: REGRADE[i] = ('A', 'position printed in the paper')
REGRADE_SRC_A = {'A_spot_blotch_017': 'RefSeq v1.0 position encoded in the GBS marker name', 'A_spot_blotch_013': 'RefSeq v1.0 position encoded in the GBS marker name'}
REGRADE_SRC_B = {'A_multi_002': 'the paper gives cM only; position is that of the 90K SNP from another source'}
REGRADE_IDS_B = {'R0353': 'the paper gives cM only; position is that of the 90K SNP from another source',
                 'R0356': 'the paper gives cM only; position is that of the 90K SNP from another source'}

R = []
def rec(src, name, ch, dis, cls, note, ev, ref='RefSeq v1.0', s=None, e=None, p=None, pve='', stat='', marker=''):
    R.append(dict(name=name, chr=ch, source=src, pve=pve, stat=stat, note=note, evidence=ev, ref=ref, curation='hand', cls=cls,
                  disease=dis, marker=marker, start=s, end=e, peak=p, grade='A'))
mb = lambda x: round(int(str(x).replace(',', '')) / 1e6, 6)

# Chen et al. 2025 JIA (Zhongmai 578 x Jimai 22 RIL), Table 2
CH = 'Zhongmai 578 x Jimai 22 RIL; black point (incl. Bipolaris sorokiniana) in five environments; interval printed in Table 2 (RefSeq v1.0).'
for q, c, a, b, s, e, r2, ld, al in [('QBp.caas-2A', '2A', 'AX-95629616', 'AX-179560327', 752.1, 753.6, '3.4-13.5', '3.3-21.5', 'Jimai 22'),
        ('QBp.caas-2B1', '2B', 'AX-111057916', 'AX-95659274', 159.4, 167.5, '19.7-28.8', '12.5-15.9', 'Zhongmai 578'),
        ('QBp.caas-2B2', '2B', 'AX-179560603', 'AX-86165074', 558.3, 563.9, '4.1-15.2', '7.3-10.4', 'Zhongmai 578'),
        ('QBp.caas-2D', '2D', 'AX-94641475', 'AX-94731323', 74.9, 84.8, '2.1-5.7', '3.7-4.3', 'Jimai 22'),
        ('QBp.caas-3A', '3A', 'AX-111258206', 'AX-95631300', 682.3, 688.9, '3.4-6.0', '2.9-3.7', 'Jimai 22'),
        ('QBp.caas-5B', '5B', 'AX-179388377', 'AX-110910964', 673.0, 679.7, '3.6-8.7', '3.4-6.4', 'Jimai 22')]:
    rec('A_spot_blotch_033', q + ' (black point)', c, BP, 'QTL', CH + ' Resistance allele from %s.' % al, 'Table 2: %s %s-%s %s-%s Mb' % (q, a, b, s, e),
        s=s, e=e, pve=r2, stat='LOD ' + ld, marker='%s, %s' % (a, b))
# Juliana et al. 2022
JU = 'GWAS in CIMMYT breeding panels; RefSeq v1.0 position encoded in the GBS marker name.'
rec('A_spot_blotch_017', '5B_688063050', '5B', BP, 'MTA', JU + ' The atlas had this marker as 5B_68806 (name cut at a line break).', 'Results text: 5B_688063050', p=688.06305, marker='5B_688063050')
rec('A_spot_blotch_017', '5B_687326022', '5B', BP, 'MTA', JU + ' Named in the results among markers explaining 2.5-3% of the variation.', 'Results text: 5B_687326022', p=687.326022, marker='5B_687326022')
# Bainsla et al. 2020: Table 3 (FarmCPU) and Table 2 (MLMM), DArTseq markers, RefSeq v1.0 bp
BA = 'Panel of 301 Afghan wheat genotypes; spot blotch AUDPC, 2017 and 2018; DArTseq marker position (RefSeq v1.0) printed in the table.'
for m, c, pos, pv, env in [('1193782', '1A', '44,038,064', '4.41', '2018'), ('1113207', '1B', '641,140,613', '2.98', '2017, 2018, mean'), ('1203688', '1B', '587,848,104', '2.71', '2018'),
        ('4909825', '1B', '557,997,411', '4.01 (FarmCPU); 8.81 (MLMM)', 'mean (FarmCPU); 2018 (MLMM)'), ('100181038', '1D', '6,202,603', '3.69', '2018, mean'), ('1215338', '2B', '653,788,494', '2.75', '2017'),
        ('1012395', '2D', '648,467,518', '2.52', '2018'), ('3029570', '3A', '142,889,812', '3.37', '2017'), ('3064675', '3B', '17,109,515', '2.59', '2017'), ('3953802', '3B', '770,969,478', '2.92', '2018'),
        ('5411703', '3B', '33,238,672', '3.86', '2017'), ('1091178', '4A', '544,628,512', '4.78', '2017'), ('100177527', '5A', '3,319,047', '4.48 (FarmCPU); 17.62 (MLMM)', '2018'),
        ('5411867', '5A', '586,600,348', '6.20 (FarmCPU); 17.66 (MLMM)', '2017, 2018, mean'), ('100120844', '5B', '43,437,349', '2.03', 'mean'), ('6036625', '5B', '407,623,507', '5.94', '2017'),
        ('3021829', '6A', '184,371,714', '3.39', '2018'), ('1128142', '7A', '586,239,026', '2.79', 'mean'), ('1001447', '7D', '607,181,532', '3.41', '2017')]:
    rec('Bainsla2020PP', 'DArT ' + m, c, BP, 'MTA', BA + ' FarmCPU; significant in %s.' % env, 'Table 3: %s %s %s' % (m, c, pos), p=mb(pos), pve=pv, marker=m)
for m, c, pos, pv, env in [('1085203', '3A', '595,935,042', '17.66', '2017, mean'), ('1128070', '3A', '502,404,423', '2.20', '2018'), ('1220348', '3A', '598,916,422', '13.24', '2017, mean'),
        ('7354241', '3B', '12,103,697', '8.83', 'mean'), ('991620', '4A', '658,343,324', '12.33', '2018'), ('998276', '5A', '569,660,176', '10.59', '2017, 2018, mean')]:
    rec('Bainsla2020PP', 'DArT ' + m, c, BP, 'MTA', BA + ' MLMM; significant in %s.%s' % (env, ' Position printed as "56,966,0,176" (read as 569,660,176).' if m == '998276' else ''),
        'Table 2: %s %s %s' % (m, c, pos), p=mb(pos), pve=pv, marker=m)
# Aggarwal et al. 2026, Table 3
AG = 'Panel of 187 spring wheat lines from India and Bangladesh; DArTseq; stable MTA (both years, at least two models); RefSeq v1.0 position printed in Table 3.'
for m, c, pos, pv, note in [('1695432', '3D', '561,275,268', '3.94-15.81', ''), ('994016', '5A', '705,576,477', '3.96-18.17', ''), ('1091615', '5A', '569,892,433', '3.86-6.09', ' Chromosome cell merged with the row above (5A).'),
        ('5411850', '3B', '684,334,097', '3.95-15.02', ''), ('1037234', '3B', '681,714,448', '4.55-17.71', ' Chromosome cell merged with the row above (3B).'),
        ('1003923', '1B', '1,295,732', '13.73-15.72', ' Chromosome cell merged with the row below (1B).'), ('1090465', '1B', '55,315,950', '15.56-18.11', '')]:
    rec('Aggarwal2026TPG', 'DArT ' + m, c, BP, 'MTA', AG + note, 'Table 3: %s %s' % (m, pos), p=mb(pos), pve=pv, marker=m)
# Lv et al. 2020, Table 2 (GWAS, 406 cultivars) and linkage QTL (UC1110 x PI610750 RIL)
LV = '406 Chinese wheat cultivars, 90K/660K SNPs; black point; multi-environment SNP (Table 2); positions on the Chinese Spring genome, taken as RefSeq v1.0.'
for m, c, pos, r2 in [('RAC875_rep_c71149_738', '2B', '782,536,180', '12.29'), ('Excalibur_c10055_221', '6A', '28,017,180', '18.20'), ('GENE-4167_145', '6A', '611,328,876', '18.20'),
        ('BS00099401_51', '6A', '613,256,522', '18.20'), ('GENE-1337_298', '2A', '107,383,862', '16.21'), ('Kukri_c17_1587', '2A', '774,850,845', '17.80'), ('BS00097147_51', '2B', '786,033,209', '16.68'),
        ('Tdurum_contig51605_1130', '3B', '822,680,599', '13.44'), ('Tdurum_contig51605_194', '3B', '822,682,075', '13.44'), ('BS00087695_51', '3B', '822,826,547', '13.44'),
        ('Tdurum_contig11114_335', '3B', '822,828,808', '13.44'), ('Tdurum_contig50364_1051', '6B', '711,702,135', '13.91'), ('GENE-4021_496', '6D', '463,743,572', '15.49'),
        ('Excalibur_rep_c77404_218', '1A', '33,205,346', '18.78'), ('D_contig25392_201', '1D', '31,425,108', '18.78'), ('RAC875_rep_c71149_798', '2B', '782,536,120', '18.78'),
        ('Kukri_c3831_2136', '4A', '634,841,828', '18.78'), ('Kukri_c22504_459', '4A', '636,574,055', '18.78'), ('Kukri_c2757_1535', '4A', '638,174,291', '18.78'),
        ('Excalibur_rep_c102533_116', '5B', '677,154,676', '18.78'), ('wsnp_CAP11_c3570_1730776', '6A', '61,409,322', '18.78'), ('CAP7_c9162_265', '6A', '61,409,514', '18.78'),
        ('Excalibur_c35871_596', '6A', '615,467,141', '14.79'), ('Tdurum_contig25789_260', '6B', '702,751,203', '18.78'), ('RAC875_c50572_106', '6B', '720,474,930', '20.55'),
        ('RAC875_c57219_1058', '6B', '720,507,929', '18.78'), ('AX-109939573', '2A', '733,030,778', '18.54')]:
    rec('Lv2020PD', m + ' (black point)', c, BP, 'MTA', LV, 'Table 2: %s %s %s' % (m, c, pos), p=mb(pos), pve=r2, marker=m)
LQ = 'UC1110 x PI610750 RIL; black point; QTL position printed in the text.'
rec('Lv2020PD', 'QBP.hau-3A (black point)', '3A', BP, 'QTL', LQ + ' Interval wms2-rPt-9057; resistance from UC1110.', 'Results: QBP.hau-3A at 60,200,972 bp', p=60.200972, pve='6.76')
rec('Lv2020PD', 'QBP.hau-6D (black point)', '6D', BP, 'QTL', LQ + ' Interval wPt-6661-wPt-1695; resistance from PI610750; covers six GWAS SNPs.', 'Results: QBP.hau-6D 5,806,474-470,627,654 bp', s=5.806474, e=470.627654, pve='7.79')
rec('Lv2020PD', 'QBP.hau-7D (black point)', '7D', BP, 'QTL', LQ + ' Resistance from UC1110.', 'Results: QBP.hau-7D 41,403,730-53,641,214 bp', s=41.40373, e=53.641214, pve='8.84')
# Tang et al. 2022, Table 1
TA = '272 Chinese wheat landraces, 55K SNP array; kernel black point; interval printed in Table 1 (reference not stated; 55K annotation on RefSeq v1.0 assumed).'
for q, c, s, e, lp, r2, mk in [('q.bp_1', '1B', 8.66, 15.51, '5.26', '9.85', 'AX-89582405, AX-110442065'), ('q.bp_2', '2A', 14.6, None, '4.84', '8.85', 'AX-109948875'),
        ('q.bp_3', '2B', 94.83, None, '4.70', '8.22', 'AX-111476445'), ('q.bp_4', '3A', 66.67, 73.87, '5.15', '9.94', 'AX-109863037, AX-110402379'),
        ('q.bp_5', '3A', 150.49, 154.39, '5.21', '10.64', 'AX-109837196, AX-110667599'), ('q.bp_6', '4A', 172.33, 181.65, '4.96', '9.37', 'AX-110416533, AX-109548239'),
        ('q.bp_7', '6A', 3.22, 6.35, '5.00', '9.19', 'AX-108877487, AX-111651180'), ('q.bp_8', '6B', 1.20, 16.60, '5.61', '10.42', 'AX-108908450, AX-110037894'),
        ('q.bp_9', '6B', 146.06, 150.72, '5.42', '9.95', 'AX-110472844, AX-108736461'), ('q.bp_10', '6B', 185.75, 187.10, '4.79', '8.75', 'AX-108921722, AX-110961220'),
        ('q.bp_11', '7A', 27.53, 38.46, '5.73', '10.99', 'AX-110667702, AX-110956682')]:
    if e is None: rec('Tang2022PD', q + ' (black point)', c, BP, 'QTL', TA, 'Table 1: %s %s %s Mb' % (q, c, s), p=s, pve=r2, stat='-log10P ' + lp, marker=mk)
    else: rec('Tang2022PD', q + ' (black point)', c, BP, 'QTL', TA, 'Table 1: %s %s %s-%s Mb' % (q, c, s, e), s=s, e=e, pve=r2, stat='-log10P ' + lp, marker=mk)
# Iqbal et al. 2022, Table 2 (RefSeq v2.0 bp)
IQ = '196 Canadian spring wheat lines; GWAS across combined environments; genomic region printed in Table 2 on IWGSC RefSeq v2.0.'
DIS = {'QCbt': CB, 'QLr': LR, 'QLs': LS, 'QYr': YR}
T2 = """QCbt.dms-1A.1 1 1A 4,383,342 4,383,342 6.9|QCbt.dms-1A.2 2 1A 13,371,025 14,027,876 7.2|QCbt.dms-1A.3 1 1A 556,873,095 556,873,095 7.4|QLs.dms-1B 1 1B 18,203,933 18,203,933 6.6|QYr.dms-1B 1 1B 540,583,488 540,583,488 7.0|QLr.dms-1B 4 1B 547,332,148 549,647,134 7.5|QCbt.dms-1B 1 1B 645,251,547 645,251,547 6.6|QCbt.dms-1D 1 1D 10,669,243 10,669,243 6.6|QLs.dms-2A.1 1 2A 19,019,241 19,019,241 7.9|QLs.dms-2A.2 6 2A 45,928,136 45,932,241 6.7|QLr.dms-2A 1 2A 758,316,082 758,316,082 7.9|QLs.dms-2B 1 2B 19,530,481 19,530,481 7.4|QLr.dms-2B.1 1 2B 690,898,021 690,898,021 6.7|QLr.dms-2B.2 10 2B 771,850,360 778,230,186 8.0|QCbt.dms-2B 1 2B 811,019,075 811,019,075 7.1|QYr.dms-2B 1 2B 812,244,914 812,244,914 7.4|QLr.dms-2D 5 2D 624,625,220 624,952,858 8.3|QCbt.dms-3A.1 1 3A 10,279,544 10,279,544 10.4|QCbt.dms-3A.2 1 3A 671,290,927 671,290,927 8.5|QLr.dms-3A 1 3A 51,644,908 51,644,908 7.1|QLs.dms-3B 1 3B 557,835,101 557,835,101 6.6|QLr.dms-3B.1 1 3B 616,148,037 616,148,037 6.6|QLr.dms-3B.2 4 3B 743,632,624 743,922,021 7.0|QLr.dms-3D 1 3D 550,280,985 550,280,985 7.0|QLs.dms-4A 1 4A 580,845,356 580,845,356 7.3|QYr.dms-4B 3 4B 638,809,518 638,813,440 8.9|QLr.dms-5A.1 2 5A 331,454,313 331,884,318 7.2|QLr.dms-5A.2 1 5A 338,666,459 338,666,459 6.8|QYr.dms-5A 1 5A 547,615,657 547,615,657 7.4|QLr.dms-5B 2 5B 281,176,164 284,717,655 7.2|QCbt.dms-5D.1 1 5D 244,100,829 244,100,829 7.0|QCbt.dms-5D.2 1 5D 565,867,455 565,867,455 6.6|QCbt.dms-6D 1 6D 7,431,984 7,431,984 16.9|QCbt.dms-7A 1 7A 15,774,259 15,774,259 6.7|QLr.dms-7B 1 7B 36,162,027 36,162,027 6.8|QYr.dms-7D 1 7D 48,955,909 48,955,909 8.6|QLs.dms-7D 1 7D 266,720,824 266,720,824 6.6"""
for row in T2.split('|'):
    q, n, c, a, b, pv = row.split()
    kw = dict(p=mb(a)) if a == b else dict(s=mb(a), e=mb(b))
    rec('Iqbal2022Plants', q, c, DIS[q.split('.')[0]], 'QTL', IQ + ' %s significant SNP(s).' % n, 'Table 2: %s %s %s-%s bp' % (q, c, a, b), ref='RefSeq v2.0', pve=pv, **kw)

SID = lambda i: 'S60-%03d' % i

def build(reports, sources):
    new = []
    for i, x in enumerate(R, 1):
        r = {k: v for k, v in x.items()}; r['id'] = SID(i)
        r['note'] = r['note'] + ' Read from the paper, Nov 2026 (gw60.py).'
        r['asm'] = lift.lift_record(r)
        new.append(r)
    keep, changed = [], []
    for r in reports:
        if r['id'] in REMOVE or r['id'].startswith('S60-'): continue
        upd = {}
        if r['source'] in RESOURCE: upd['source'] = RESOURCE[r['source']]
        if r['id'] in REGRADE and r.get('grade') != 'A': upd['grade'] = 'A'
        if r['source'] in REGRADE_SRC_A and r.get('grade') == 'B' and r.get('peak') is not None: upd['grade'] = 'A'
        if (r['source'] in REGRADE_SRC_B and r.get('grade') == 'A') or r['id'] in REGRADE_IDS_B: upd['grade'] = 'B'
        if upd:
            r = dict(r); r.update(upd); changed.append((r['id'], upd))
        keep.append(r)
    src = dict(sources); src.update(NEW_SOURCES)
    return keep + new, src, new, [i for i in REMOVE if any(x['id'] == i for x in reports)], changed

if __name__ == '__main__':
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], M['sources'], new, rem, ch = build(M['reports'], M['sources'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</'))
    Sd['reports'], Sd['sources'], new2, rem2, ch2 = build(Sd['reports'], Sd['sources'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs, 'removed': rem2, 'changed': ch2, 'sources': NEW_SOURCES}, open('/home/claude/work/gw60_db.json', 'w'), ensure_ascii=False)
    import collections
    print('new', len(new), 'removed', rem2, 'changed', len(ch2), collections.Counter(tuple(sorted(u.items())) for _, u in ch2),
          collections.Counter((r['asm'].get('RefSeq v2.1') or {}).get('flag', 'none') for r in new))
