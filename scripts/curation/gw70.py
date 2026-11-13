"""Batch 7 (Nov 2026): STB, SNB and tan spot GWAS papers from the user's GWAS/STB, GWAS/TS, GWAS/GWAS_SNB and GWAS/SGB folders.
Papers already in the atlas were left unchanged; this adds the papers that were missing.
- grade A: position printed in the paper (Chaudhary 2025 on RefSeq v2.1, Berisha 2024 on v2.0, the rest v1.0)
- grade B: 90K/35K/DArT marker placed with the sb26 resolver (paper gives cM only)
- grade C: locus tagged in the paper to a known gene (Phan 2018); grade D: marker not placeable
Writes data/merged_v13f.json, web/index.html (seed) and /home/claude/work/gw70_db.json."""
import json, re, sys, io, contextlib
sys.path.insert(0, 'scripts/curation'); sys.path.insert(0, 'scripts/curation/sb19')
with contextlib.redirect_stdout(io.StringIO()):
    import sb26_fix as F, stb44
import lift

STB, SNB, TS = 'Septoria tritici blotch', 'Septoria nodorum blotch', 'Tan spot'
def S(cit, doi, short): return dict(citation=cit, url=('https://doi.org/' + doi) if doi else '', short=short)
NEW_SOURCES = {
 'Alemu2021FG': S('Alemu A, Brazauskas G, Gaikpa DS, Henriksson T, Islamov B, Jorgensen LN, Koppel M, Koppel R, Liatukas Z, Svensson JT, Chawade A (2021) Genome-wide association analysis and genomic prediction for adult-plant resistance to Septoria tritici blotch and powdery mildew in winter wheat. Frontiers in Genetics 12:661742', '10.3389/fgene.2021.661742', 'Alemu et al. (2021)'),
 'Aouini2025bR': S('Aouini L, Gaire R, Scofield S, Brown-Guedira G, Mohammadi M, Goodwin SB (2025) Genome-wide association study in a US soft winter wheat population reveals novel and known sources of resistance to the Septoria tritici blotch pathogen Zymoseptoria tritici. bioRxiv 2025.09.30.679614', '10.1101/2025.09.30.679614', 'Aouini et al. (2025)'),
 'Chaudhary2025BMC': S('Chaudhary S, Zakieh M, Dubey M, Jensen DF, Grenville-Briggs L, Chawade A, Karlsson M (2025) Plant genotype-specific modulation of Clonostachys rosea-mediated biocontrol of septoria tritici blotch disease in wheat. BMC Plant Biology 25:576', '10.1186/s12870-025-06620-9', 'Chaudhary et al. (2025)'),
 'Kassie2025PO': S('Kassie MM, Abebe TD et al. (2025) Genome-wide association analysis of Septoria tritici blotch for adult plant resistance in elite bread wheat (Triticum aestivum L) genotypes. PLoS ONE 20(3):e0317603', '10.1371/journal.pone.0317603', 'Kassie et al. (2025)'),
 'Kidane2017FPS': S('Kidane YG, Hailemariam BN, Mengistu DK, Fadda C, Pe ME, Dell\'Acqua M (2017) Genome-wide association study of Septoria tritici blotch resistance in Ethiopian durum wheat landraces. Frontiers in Plant Science 8:1586', '10.3389/fpls.2017.01586', 'Kidane et al. (2017)'),
 'Mahboubi2022JAG': S('Mahboubi M, Talebi R, Mehrabi R, Naji AM, Maccaferri M, Kema GHJ (2022) Genetic analysis of novel resistance sources and genome-wide association mapping identified novel QTLs for resistance to Zymoseptoria tritici, the causal agent of septoria tritici blotch in wheat. Journal of Applied Genetics 63:429-445', '10.1007/s13353-022-00696-x', 'Mahboubi et al. (2022)'),
 'Mikaberidze2025MPMI': S('Mikaberidze A, McDonald BA, Kronenberg L (2025) A genome-wide association study identifies markers and candidate genes affecting tolerance to the wheat pathogen Zymoseptoria tritici. Molecular Plant-Microbe Interactions 38:265-274', '10.1094/MPMI-08-24-0085-FI', 'Mikaberidze et al. (2025)'),
 'Odilbekov2019FG': S('Odilbekov F, Armoniene R, Koc A, Svensson J, Chawade A (2019) GWAS-assisted genomic prediction to predict resistance to Septoria tritici blotch in Nordic winter wheat at seedling stage. Frontiers in Genetics 10:1224', '10.3389/fgene.2019.01224', 'Odilbekov et al. (2019)'),
 'Yates2019PP': S('Yates S, Mikaberidze A, Krattinger SG, Abrouk M, Hund A, Yu K, Studer B, Fouche S, Meile L, Pereira D, Karisto P, McDonald BA (2019) Precision phenotyping reveals novel loci for quantitative resistance to Septoria tritici blotch. Plant Phenomics 2019:3285904', '10.34133/2019/3285904', 'Yates et al. (2019)'),
 'Miedaner2013BMC': S('Miedaner T, Zhao Y, Gowda M, Longin CFH, Korzun V, Ebmeyer E, Kazman E, Reif JC (2013) Genetic architecture of resistance to Septoria tritici blotch in European wheat. BMC Genomics 14:858', '10.1186/1471-2164-14-858', 'Miedaner et al. (2013)'),
 'Seed2023Thesis': S('Seed PM (2023) Mapping and characterisation of resistance to Septoria tritici blotch in winter wheat. PhD thesis, University of Nottingham, UK', '', 'Seed (2023)'),
 'Liu2015PD': S('Liu Z, El-Basyoni I, Kariyawasam G, Zhang G, Fritz A, Hansen J, Marais F, Friskop A, Chao S, Akhunov E, Baenziger PS (2015) Evaluation and association mapping of resistance to tan spot and Stagonospora nodorum blotch in adapted winter wheat germplasm. Plant Disease 99:1333-1341', '10.1094/PDIS-11-14-1131-RE', 'Liu et al. (2015)'),
 'Phan2018TAG': S('Phan HTT, Rybak K, Bertazzoni S, Furuki E, Dinglasan E, Hickey LT, Oliver RP, Tan KC (2018) Novel sources of resistance to Septoria nodorum blotch in the Vavilov wheat collection identified by genome-wide association studies. Theoretical and Applied Genetics 131:1223-1238', '10.1007/s00122-018-3073-y', 'Phan et al. (2018)'),
 'Berisha2024PB': S('Berisha P, Michel S, Loschenberger F, Ametz C, Bistrich H, Burstmayr H (2024) Genome-wide association mapping and genomic prediction of Septoria nodorum blotch resistance in Central European winter wheat (Triticum aestivum L.). Plant Breeding', '10.1111/pbr.13233', 'Berisha et al. (2024)'),
 'Mandal2024Genes': S('Mandal R, He X, Singh G, Kabir MR, Joshi AK, Singh PK (2024) Screening of CIMMYT and South Asian bread wheat germplasm reveals marker-trait associations for seedling resistance to Septoria nodorum blotch. Genes 15:890', '10.3390/genes15070890', 'Mandal et al. (2024)'),
}

R = []
def recA(src, name, ch, dis, cls, note, ev, ref='RefSeq v1.0', s=None, e=None, p=None, pve='', stat='', marker=''):
    R.append(dict(name=name, chr=ch, source=src, pve=pve, stat=stat, note=note, evidence=ev, ref=ref, curation='hand', cls=cls,
                  disease=dis, marker=marker, start=s, end=e, peak=p, grade='A'))
def recM(src, name, ch, dis, cls, note, ev, markers=(), peak=None, pve='', stat=''):
    r = dict(name=name, chr=ch, source=src, pve=pve, stat=stat, note=note, evidence=ev, ref='RefSeq v1.0', curation='hand', cls=cls,
             disease=dis, marker=', '.join([m for m in markers] + ([peak] if peak else [])), start=None, end=None, peak=None, grade=None)
    F.pos = stb44.P
    try: F.place(r, list(markers), ch, peak_marker=peak)
    finally: F.pos = stb44.POS0
    if r['peak'] is None and r['start'] is None: r['grade'] = 'D'
    elif peak and r['peak'] is not None: r['grade'] = 'B'
    for k in ('start', 'end', 'peak'):
        if r[k] is not None: r[k] = round(float(r[k]), 4)
    R.append(r)
mb = lambda x: round(int(str(x).replace(',', '')) / 1e6, 6)

# ---- STB ----
AL = 'Nordic-Baltic winter wheat panel; adult-plant STB in eight environments; 90K SNP, the paper gives cM only.'
for m, c, cm, env in [('wsnp_Ex_c33012_41567026', '4A', 153, 'E1, E6, E7, E8'), ('BS00021714_51', '1A', 78, 'E2, E6, E7'), ('RAC875_c68120_285', '5A', 63, 'E2, E6, E7'),
        ('wsnp_Ku_c38451_47086066', '6A', 79, 'E4, E8'), ('wsnp_CAP11_c59_99769', '3B', 115, 'E2, E6'), ('D_F5XZDLF01A85DT_301', '1D', 61, 'E7, E8'), ('Ra_c69221_1167', '5A', 42, 'E3'),
        ('RAC875_rep_c116515_181', '3B', 71, 'E8'), ('RAC875_c14309_317', '6B', 76, 'E1'), ('RAC875_c23654_214', '6B', 113, 'E4'), ('Excalibur_rep_c69522_83', '1B', 171, 'E3'),
        ('Kukri_rep_c103893_875', '2B', 65, 'E5')]:
    recM('Alemu2021FG', m, c, STB, 'MTA', AL + ' Significant in %s.' % env, 'Table 2: %s %s %s cM' % (m, c, cm), peak=m)
AO = 'US soft red winter wheat panel; seedling STB with isolates T2 and T46; GBS SNP named by its RefSeq v1.0 position.'
for m, tr in [('S1D_8183618', 'T2 pycnidia; Stb19 region'), ('S1D_8183661', 'T2 pycnidia; Stb19 region'), ('S6A_611714334', 'T46 necrosis'), ('S1A_500074551', 'T46 pycnidia'),
        ('S1A_577377232', 'T46 pycnidia'), ('S4A_605028654', 'T46 pycnidia'), ('S4A_729925867', 'T46 pycnidia; Stb7 region'), ('S6A_48048307', 'T46 pycnidia'),
        ('S6A_532657432', 'T46 pycnidia; Stb15 region (paper)'), ('S7A_719723711', 'T46 pycnidia'), ('S3B_726510291', 'T46 pycnidia'), ('S5B_519034813', 'T46 pycnidia; Stb1 region (paper)')]:
    c, p = m[1:3], int(m.split('_')[1])
    recA('Aouini2025bR', m, c, STB, 'MTA', AO + ' %s.' % tr, 'Table 6: %s' % m, p=mb(p), marker=m)
CU = 'Swedish winter wheat panel; STB severity with Z. tritici alone (Zt) or with the biocontrol fungus Clonostachys rosea (ZtCr); position (RefSeq v2.1) printed in Table 2.'
for m, c, pos, pv, tr in [('Excalibur_c49875_479', '2B', '788,524,048', '4.0', 'Zt'), ('IAAV4876', '3B', '3,492,598', '6.6', 'Zt'), ('Kukri_rep_c70198_1436', '3B', '495,024,661', '5.1', 'Zt'),
        ('Excalibur_c29625_222', '3B', '493,110,618', '7.6', 'Zt'), ('RAC875_rep_c83245_239', '3B', '493,776,130', '7.2', 'Zt'), ('BS00022902_51', '1B', '73,310,732', '3.8', 'ZtCr'),
        ('BS00070856_51', '6D', '494,584,503', '0.7', 'ZtCr')]:
    recA('Chaudhary2025BMC', m, c, STB, 'MTA', CU + ' Trait %s.' % tr, 'Table 2: %s %s %s bp' % (m, c, pos), ref='RefSeq v2.1', p=mb(pos), pve=pv, marker=m)
KA = 'Ethiopian elite bread wheat panel (186 genotypes), 20K SNP array; adult-plant STB severity at heading (H) and maturity (M); position printed in Table 4.'
T4 = """qtSTB.01 1B 21473796 21473796 1 Adet 2023 H|qtSTB.02 1B 563254727 598339967 4 Adet 2022 M, Adet 2023 H|qtSTB.03 1B 649005218 669948745 3 Dabat 2023 M, Koga 2023 M|qtSTB.04 2A 4620089 21095554 4 Koga 2023 M|qtSTB.05 2A 102685812 102685812 1 Koga 2023 M|qtSTB.06 2A 719566734 719566734 1 Koga 2023 H|qtSTB.07 2B 1832830 1832830 1 Kulumsa 2023 H|qtSTB.08 2B 28138249 42619355 2 Koga 2023 M|qtSTB.09 2B 57876972 72456453 2 Dabat 2023 M|qtSTB.10 2B 150464672 175028367 7 Koga 2023 M, Dabat 2023 M, Dabat 2023 H|qtSTB.11 2B 200385458 200385458 1 Koga 2023 M|qtSTB.12 2B 565080256 579821341 2 Kulumsa 2023 H|qtSTB.13 2D 8297177 8297177 1 Kulumsa 2023 H|qtSTB.14 3B 252595780 252595780 1 Adet 2023 M|qtSTB.15 3B 623076245 644321527 3 Kulumsa 2023 M, Kulumsa 2023 H|qtSTB.16 3D 563861601 563861601 1 Koga 2023 M|qtSTB.17 4A 17265696 17265696 1 Dabat 2023 M|qtSTB.18 4A 616874298 616874298 1 Adet 2023 H|qtSTB.19 4D 403572723 403572723 1 Koga 2023 H|qtSTB.20 5A 205132950 205132950 2 Dabat 2023 H, Koga 2023 M|qtSTB.21 5A 570714642 570714642 1 Adet 2022 H|qtSTB.22 5A 577942781 577942781 1 Adet 2022 H|qtSTB.23 5A 581620919 597890151 41 Adet 2022 M and H, Adet 2023 H, Dabat 2023 H and M (stable)|qtSTB.24 5B 50469684 50469684 1 Dabat 2023 H|qtSTB.25 5B 307408948 307408948 1 Adet 2023 H|qtSTB.26 5B 409876449 409876449 1 Adet 2023 H|qtSTB.27 5B 424925845 441798211 2 Koga 2023 H|qtSTB.28 5D 451674901 451674901 1 Adet 2022 M|qtSTB.29 5D 470473562 480049293 2 Dabat 2023 M, Adet 2022 H|qtSTB.30 6A 560812338 567805780 3 Adet 2023 H|qtSTB.31 6A 594818372 594818372 1 Dabat 2023 M|qtSTB.32 6A 598636889 605651331 2 Adet 2022 M|qtSTB.33 6D 463743278 463743278 2 Kulumsa 2023 H and M|qtSTB.34 7A 977272 8132084 4 Adet 2023 M|qtSTB.35 7A 22834490 22834490 1 Koga 2023 M|qtSTB.36 7A 704792521 704792521 1 Adet 2023 H|qtSTB.37 7B 491397671 491397671 1 Adet 2022 M|qtSTB.38 7B 658449259 709463743 5 Dabat 2023 M, Koga 2023 H, Kulumsa 2023 H (stable)|qtSTB.39 7D 386646259 386646259 2 Adet 2022 M"""
for row in T4.split('|'):
    q, c, a, b, n, env = row.split(' ', 5)
    kw = dict(p=mb(a)) if a == b else dict(s=mb(a), e=mb(b))
    recA('Kassie2025PO', q, c, STB, 'QTL', KA + ' %s MTA(s); %s.' % (n, env), 'Table 4: %s %s %s-%s bp' % (q, c, a, b), **kw)
KI = 'Ethiopian durum wheat landraces; STB severity (SDS) and progress coefficient (SPC); 90K SNP, the paper gives cM only.'
recM('Kidane2017FPS', 'qSTB.1', '1A', STB, 'QTL', KI + ' SDS at heading.', 'Table 2: qSTB.1 67.7-69.2 cM IACX3496/Kukri_c10239_2186', markers=('IACX3496', 'Kukri_c10239_2186'))
for q, c, m, cm, tr in [('qSTB.2', '2B', 'Jagger_c9472_305', '85.8', 'SDS at heading'), ('qSTB.4', '4A', 'Excalibur_c113341_139', '167.5', 'SPC'), ('qSTB.5', '5A', 'Tdurum_contig10210_425', '135.2', 'SPC, 2013')]:
    recM('Kidane2017FPS', q, c, STB, 'QTL', KI + ' %s.' % tr, 'Table 2: %s %s cM %s' % (q, cm, m), peak=m)
recM('Kidane2017FPS', 'qSTB.3', '3A', STB, 'QTL', KI + ' SPC; strongest QTL (9 MTAs).', 'Table 2: qSTB.3 71.6-72.5 cM wsnp_Ex_c3478_6369892/wsnp_RFL_Contig4404_5157920',
     markers=('wsnp_Ex_c3478_6369892', 'wsnp_RFL_Contig4404_5157920'))
MA = '185 globally diverse bread wheat genotypes; seedling STB with 10 Z. tritici isolates; DArTseq; RefSeq v1.0 position printed in Table 5.'
T5 = """Qstb.iau-1 5332931 RM6 1A 2,187,418|Qstb.iau-2 4544165 IPO99031 2A 748,128,615|Qstb.iau-2 3533473 IPO99031 2A 776,037,491|Qstb.iau-3 3064517 IPO323 2B 11,132|Qstb.iau-3 1864355 IPO323,IPO86013 2B 11,066|Qstb.iau-4 1093912 IPO98022 2B 775,155,639|Qstb.iau-5 1265720 RM230 2D 12,587,173|Qstb.iau-5 3937084 RM6 2D 22,776,052|Qstb.iau-6 1205624 IPO98022 3B 310,392,089|Qstb.iau-6 1124803 IPO323 3B 557,176,126|Qstb.iau-7 5582113 IPO98022,IPO99031,RM6,RM183 4A 140,700,686|Qstb.iau-8 2257551 RM22,RM183 4A 641,506,835|Qstb.iau-9 1238557 RM230 4A 722,708,344|Qstb.iau-10 2258488 RM183 5A 568,491,318|Qstb.iau-11 5970385 IPO02166 5B 488,112,132|Qstb.iau-11 1088825 IPO02166 5B 490,053,583|Qstb.iau-12 2354562 IPO92034 5B 557,353,675|Qstb.iau-13 6038202 IPO86013 5D 541,681,037|Qstb.iau-13 7157166 IPO86013 5D 541,902,658|Qstb.iau-14 1110173 IPO92034 6A 5,125,140|Qstb.iau-15 1009838 IPO02166 6B 129,858,768|Qstb.iau-15 1266810 IPO98022 6B 648,849,374|Qstb.iau-15 5577074 IPO98022 6B 655,271,759|Qstb.iau-16 2275399 IPO98022 6D 433,578,091|Qstb.iau-17 4008741 IPO86013 7A 670,929,737|Qstb.iau-18 3020733 RM22 7D 15,219,082|Qstb.iau-19 2242097 IPO323 7D 556,246,143"""
for row in T5.split('|'):
    q, m, iso, c, pos = row.split()
    recA('Mahboubi2022JAG', '%s (DArT %s)' % (q, m), c, STB, 'MTA', MA + ' Isolate(s) %s.' % iso.replace(',', ', '), 'Table 5: %s %s %s %s bp' % (q, m, c, pos), p=mb(pos), marker=m)
recA('Mikaberidze2025MPMI', 'Ra_c7322_2294 (resistance)', '5A', STB, 'MTA', 'Swiss elite winter wheat panel (GABI), automated leaf image analysis; resistance MTA (the tolerance MTAs are not resistance loci and are not listed); RefSeq v1.0 position printed in Table 1. Same region as Yates et al. (2019) interval 2.',
     'Table 1: Resistance 5A Ra_c7322_2294 20,816,495', p=mb('20,816,495'), pve='12.72', stat='LOD 8.12', marker='Ra_c7322_2294')
OD = 'Nordic winter wheat panel (175 genotypes); seedling STB; RefSeq v1.0 position printed in Table 1.'
for m, c, pos, r2 in [('BobWhite_c1361_1187', '1A', 1525253, ''), ('BobWhite_c42716_71', '1B', 623712765, '0.11'), ('wsnp_Ex_rep_c66255_64400455', '1B', 623729791, '0.09'),
        ('RFL_Contig5937_1677', '1B', 623730512, '0.08'), ('RAC875_c47427_75', '1B', 623731255, '0.10'), ('Excalibur_rep_c72368_68', '1B', 623770763, '0.09'), ('BS00110231_51', '1B', 623989423, '0.09'),
        ('wsnp_Ex_c22423_31615798', '2B', 215593752, ''), ('wsnp_Ex_c5929_10402147', '3A', 481018206, ''), ('Excalibur_c17553_84', '5A', 375375809, '')]:
    recA('Odilbekov2019FG', m, c, STB, 'MTA', OD, 'Table 1: %s %s %d' % (m, c, pos), p=mb(pos), pve=r2, marker=m)
YT = 'Swiss elite winter wheat panel (GABI); STB measured by automated leaf image analysis; interval (RefSeq v1.0) printed in Table 1.'
T1 = """PLACL 1 5A 644275 1225606 5.51 0.065|PLACL 2 5A 20816445 20998554 5.95 0.068|PLACL 3 5B 671238689 671360289 5.93 0.089|PLACL 4 5D 1612786 2664314 9.16 0.103|PDL 5 2B 648912764 651932324 5.93 0.078|PDL 6 2B 781584349 781830297 5.59 0.085|PDL 7 4A 648893451 658970236 6.37 0.021|PDL 8 5D 52615122 60899123 6.13 0.044|PDL 9 5D 403075850 413073376 6.19 0.066|PDL 10 6A 411095227 415905355 6.24 0.019|PDL 11 6A 421527350 425277806 6.47 0.091|PDL 12 6B 175719208 176403432 5.79 0.073|PDL 13 6B 453931609 455449000 6.03 0.062|PDL 14 6B 462130236 468736016 6.99 0.052|PDL 15 6B 472609890 474177953 7.10 0.093|PDL 16 6D 394936301 406340008 5.70 0.074|PDL 17 7B 265441725 325980844 5.46 0.034|MPA 18 2B 245888219 247935045 5.50 0.059|PGV 19 1A 12369332 12506454 6.64 0.083|PGV 20 1A 472140874 472168930 5.53 0.063|PGV 21 2A 635581134 639988522 6.12 0.078|PGV 22 3B 30319666 32286228 5.60 0.075|PGV 23 4D 485524935 502402946 8.71 0.106|PGV 24 5A 685438196 685537900 5.67 0.075|PGV 25 5B 442374822 455735533 6.40 0.096|PGV 26 7B 65661891 83030127 6.34 0.078"""
TR = {'PLACL': 'leaf area covered by lesions', 'PDL': 'pycnidia density within lesions', 'MPA': 'mean pycnidia area', 'PGV': 'pycnidia grey value'}
for row in T1.split('|'):
    t, i, c, a, b, lod, r2 = row.split()
    recA('Yates2019PP', 'Interval %s (%s)' % (i, t), c, STB, 'QTL', YT + ' Trait: %s.' % TR[t], 'Table 1: interval %s %s %s-%s bp' % (i, c, a, b),
         s=mb(a), e=mb(b), pve=r2, stat='LOD ' + lod)
for m, c, cm in [('we00819', '1B', '18'), ('we02521', '6A', '14.2'), ('we03018', '5B', '53'), ('we03636', '2B', '74.2')]:
    R.append(dict(name=m, chr=c, source='Miedaner2013BMC', pve='', stat='', ref='RefSeq v1.0', curation='hand', cls='MTA', disease=STB, marker=m,
                  note='Hybrid and inbred winter wheat (GABI-Wheat), field STB; proprietary SNP with a cM position only (%s cM); not placeable on the genome, left unplaced.' % cm,
                  evidence='Table 2: %s %s (%s cM)' % (m, c, cm), start=None, end=None, peak=None, grade='D'))
SE = 'Cougar x Scout DH population; field STB assessment; the thesis gives the QTL interval in cM; position is that of the Axiom peak marker.'
for c, yr, iv, m, lod, pv in [('4D', 2018, '15.8-17.5', 'AX-94733088', '10.86', '44.1'), ('1B', 2019, '39.6-42.6', 'AX-95652713', '5.93', '16'), ('2A', 2019, '56.3-58.4', 'AX-95171275', '3.6', '17.6'),
        ('2B', 2019, '66.3-68.5', 'AX-95230317', '4.22', '10.9'), ('2B', 2019, '228.2-228.5', 'AX-94814715', '3.12', '15.4'), ('3B', 2019, '240.1', 'AX-94701598', '3.22', '15.8')]:
    recM('Seed2023Thesis', 'QStb.%s (Cougar x Scout, %s cM)' % (c, iv.split('-')[-1]), c, STB, 'QTL', SE + ' %d; interval %s cM.%s' % (yr, iv, ' Major 4D resistance from Cougar.' if c == '4D' else ''),
         'Chapter 2 Table 5: %s %s %s LOD %s' % (c, iv, m, lod), peak=m, pve=pv, stat='LOD ' + lod)

# ---- SNB and tan spot ----
LI = 'US hard winter wheat cultivars; seedling inoculation; 90K SNP (IWB names), the paper gives cM only.'
for iso, dis, rows in [('Sn4', SNB, [('IWB67424', '5A'), ('IWB36366', '5B'), ('IWB38178', '5B'), ('IWB45668', '5D')]),
        ('Pti2 (race 1)', TS, [('IWB20112', '2D'), ('IWB11339', '4A'), ('IWB47813', '4A'), ('IWB45442', '5D'), ('IWB68544', '7A')]),
        ('331-9 (race 3)', TS, [('IWB64643', '2A'), ('IWB74851', '2A'), ('IWB9575', '3B'), ('IWB51094', '3B'), ('IWB64735', '3D'), ('IWB35182', '4A'), ('IWB51002', '7A')]),
        ('DW5 (race 5)', TS, [('IWB64607', '3B'), ('IWB6788', '5D'), ('IWB70366', '6B'), ('IWB45141', '7B')]),
        ('AR CrossB10', TS, [('IWB64643', '2A'), ('IWB35182', '4A'), ('IWB3978', '6A'), ('IWB62163', '6B')])]:
    for m, c in rows:
        recM('Liu2015PD', '%s (%s)' % (m, iso.split()[0]), c, dis, 'MTA', LI + ' Isolate %s.%s' % (iso, ' Same marker as the ToxA (Tsn1) association.' if m == 'IWB36366' else ''),
             'Table 4: %s %s %s' % (iso, m, c), peak=m)
PH = 'Vavilov wheat collection; seedling SNB with P. nodorum SN15 or the ToxA-deleted strain toxa13; DArT marker, the paper gives cM only.'
for m, c, gene, cm, iso in [('1268227', '1B', 'Snn1', '7.025', 'SN15'), ('1195188', '1B', 'Snn1', '7.025', 'SN15'), ('2248448', '2D', 'Snn2', '195.48', 'SN15'),
        ('1254900', '3A', '', '209.7', 'SN15'), ('1676515', '4A', '', '240.629', 'SN15'), ('1104887', '4B', '', '5.495', 'SN15'), ('3939782', '5A', '', '40.05', 'SN15'),
        ('1126921', '5B', 'Snn3-B1', '4.119', 'SN15'), ('1151694', '5B', 'Snn3-B1', '5.452', 'SN15'), ('1232949', '5B', 'Snn3-B1', '5.452', 'SN15'), ('2276671', '5B', 'Snn3-B1', '6.153', 'SN15'),
        ('2374474', '7A', '', '110.881', 'SN15'), ('1231149', '2D', 'Snn7', '140.67', 'toxa13'), ('1126530', '5A', '', '72.798', 'toxa13'), ('1118447', '7D', '', '231.408', 'toxa13')]:
    nt = PH + ' Isolate %s.' % iso + (' The paper assigns this locus to %s.' % gene if gene else '')
    recM('Phan2018TAG', 'DArT %s (%s)' % (m, iso), c, SNB, 'MTA', nt, 'Table 1: %s %s %s cM' % (m, c, cm), peak=m)
    r = R[-1]
    if r['grade'] == 'D' and gene in F.GENE and F.GENE[gene][0] == c:
        r['peak'] = round(F.GENE[gene][1], 4); r['grade'] = 'C'
        r['note'] = r['note'].replace(' No marker of this locus could be placed on %s; left unplaced.' % c, '') + ' Marker not placeable; position is that of %s, to which the paper assigns the locus.' % gene
BE = 'Central European winter wheat breeding lines (Saatzucht Donau), field SNB; DArTseq; RefSeq v2.0 position printed in the results text. The other MTAs are in Table S1.'
recA('Berisha2024PB', 'DArT 1102861', '2A', SNB, 'MTA', BE + ' Most important peak marker across data sets (-log10 P 5.9).', 'Results: 1102861 at 19,730,738 bp on 2A', ref='RefSeq v2.0', p=mb('19,730,738'), marker='1102861')
recA('Berisha2024PB', 'DArT 5411524', '2B', SNB, 'MTA', BE + ' Peak marker with the highest probability.', 'Results: 5411524 at 24,829,137 bp on 2B', ref='RefSeq v2.0', p=mb('24,829,137'), marker='5411524')
MD = 'CIMMYT and South Asian bread wheat lines; seedling SNB; position printed in Table 2 (RefSeq v1.0: the printed positions match the v1.0 positions of these array SNPs).'
for m, c, pos, r2 in [('AX-94677640', '5B', '6,838,168', '8.48-8.96'), ('AX-94864577', '5B', '6,974,754', '8.63-9.68'), ('AX-94978939', '3D', '52,806,848', '7.35-7.65'),
        ('BS00091519_51', '5B', '6,648,517', '8.44-8.80'), ('Excalibur_c47452_183', '5B', '6,654,116', '7.83-8.22'), ('GENE-3324_338', '5B', '3,521,225', '7.40-7.76'),
        ('IACX9261', '5B', '546,703,936', '10.31-13.30'), ('Tdurum_contig25513_123', '5B', '565,753,503', '8.59-15.23'), ('Tdurum_contig31131_198', '5B', '417,443,545', '7.86-11.34'),
        ('tplb0027f13_1493', '5B', '565,754,215', '7.68-15.67'), ('wsnp_Ex_c210_411604', '5A', '584,054,461', '7.50'), ('wsnp_Ex_c29304_38355434', '5A', '584,451,511', '7.58')]:
    recA('Mandal2024Genes', m, c, SNB, 'MTA', MD, 'Table 2: %s %s %s' % (m, c, pos), p=mb(pos), pve=r2, marker=m)

SID = lambda i: 'S70-%03d' % i

def build(reports, sources):
    new = []
    for i, x in enumerate(R, 1):
        r = dict(x); r['id'] = SID(i)
        r['note'] = (r['note'] + ' Read from the paper, Nov 2026 (gw70.py).').strip()
        r['asm'] = lift.lift_record(r)
        if (r['asm'].get('RefSeq v2.1') or {}).get('flag') == 'wide_gap':
            import ts52
            v = r['asm']['RefSeq v2.1']; v.update(ts52.wide_fix(r))
        new.append(r)
    keep = [r for r in reports if not r['id'].startswith('S70-')]
    src = dict(sources); src.update(NEW_SOURCES)
    return keep + new, src, new

if __name__ == '__main__':
    import collections
    M = json.load(open('data/merged_v13f.json'))
    M['reports'], M['sources'], new = build(M['reports'], M['sources'])
    json.dump(M, open('data/merged_v13f.json', 'w'))
    H = open('web/index.html').read()
    m = re.search(r'(<script type="application/json" id="seed">)(.*?)(</script>)', H, re.S)
    Sd = json.loads(m.group(2).replace('<\\/', '</'))
    Sd['reports'], Sd['sources'], new2 = build(Sd['reports'], Sd['sources'])
    ser = lambda o: json.dumps(o, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
    open('web/index.html', 'w').write(H[:m.start(2)] + ser(Sd) + H[m.end(2):])
    docs = []
    for r in new2:
        d = {k: v for k, v in r.items() if v not in (None, '') and k != 'asm'}
        d['asm'] = {k: r['asm'][k] for k in ('RefSeq v2.1', 'RefSeq v1.0', '_note') if k in r['asm']}
        docs.append(d)
    json.dump({'docs': docs, 'sources': NEW_SOURCES}, open('/home/claude/work/gw70_db.json', 'w'), ensure_ascii=False)
    print('new', len(new), collections.Counter(r['grade'] for r in new), collections.Counter(r['source'] for r in new),
          collections.Counter((r['asm'].get('RefSeq v2.1') or {}).get('flag', 'none') for r in new))
