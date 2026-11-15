# Correction log: SNB and tan spot papers (Archive.zip) and Röder 1998 GWM primers

November 2026. Positions are on IWGSC RefSeq v1.0 (atlas frame) with RefSeq v2.1 shown in the atlas. Scripts: `scripts/curation/sb26_fix.py`, `sb26_apply.py`, `r98_fix.py`, `r98_apply.py`.

## Part 1. Septoria nodorum blotch and tan spot papers

The zip held 29 PDFs, which are 26 unique papers (Abeysekara 2011/2012 appear three times, Gao 2015 twice, and Shi 2015.1 is the same file as Shi 2015b). Every QTL table was re-read. The 116 old atlas records from these papers were removed and replaced by 182 records (S26-001 to S26-182): grade A 34, B 44, C 58, D 46. All placed records have a RefSeq v2.1 position with flag ok, rechecked, reoriented_checked or v2.0_as_reported.

### What was wrong in the old records

| Paper | Problem found | Now |
|---|---|---|
| Lin et al. (2020), NIAB MAGIC | 19 "MTA" records stored cM values as Mb; wrong ends for 2A.3, 2A.4, 3A, 5B.1, 7D.3; 6A.2 peak was a cM value; QSnb.niab-2D.1 missing | 19 QTL from Table 3 bp positions (grade A). QSnb.niab-7D.3: the printed end, 724.13 Mb, lies beyond the end of 7D and is the 7A position of RAC875_c10022_23, so the locus sits at the 7D flank JD_c2708_1512 (grade B) |
| Lin et al. (2021), BMWpop MAGIC | 3 records were Lin 2020 loci; point positions wrong; 5 QTL missing; citation year shown as 2020 | 10 QTL from Table 2 bp (grade A), new citation key Lin2021TAG |
| Singh et al. (2019) | positions not from the paper | 7 tan spot and 5 SNB QTL; XFCP393 flank placed at Tsn1 |
| Schnurbusch et al. (2003) | incomplete | 13 QSng loci from the SSR flanks |
| Shankar et al. (2008) | incomplete | 28 records, each with its peak marker |
| Rivera-Burgos et al. (2022) | positions not from the paper | 12 records; Qncb.snl-5Ac at the printed v2.0 position (579.7-579.8 Mb); Ppd-A1 and Ppd-D1 loci at the genes; the rest unplaced (markers only in Fig. 2) |
| Czembor et al. (2019) | positions not checked | 5 QTL |
| Abeysekara et al. (2012) | 5 duplicate records for 2 QTL | QSnb.fcu-1AS at Snn4; QSnb.fcu-7AS at XBF293121 (primers from the paper, 7A 61.43 Mb) |
| Chu et al. (2010) | incomplete | Snn1 and Tsn1 loci plus 5 QSnb.fcu loci |
| Francki et al. (2011) | 2A loci merged | 5 records, 2A in 2004 and 2005 kept apart |
| Lu and Lillemo (2014) | positions not checked | 8 closest-marker QTL |
| Li et al. (2021) | IWB27379 record not in the paper's table | 5 R-gene MTAs at printed bp plus Xfcp623 at Tsn1 |
| Thapa et al. (2016) | its only record was Schnurbusch's QSng.sfr-3BS | 1 record at Xgwm389 (0.81 Mb) |
| Hu et al. (2019) | only record was another study's QTL named in the text | 20 tan spot and SNB QTL from Tables 2 and 3 |
| Francki et al. (2018) | positions not checked | 2D glume QTL IWB31450-IWB21124 |
| Ruud et al. (2017) | incomplete | 16 QTL (5 seedling, 11 field) |
| Arseniuk et al. (2004) | wrong PVE | QSnl.ihar-6A at Xgwm570 (14% and 36%); putative 6DL QTL at barc175 |
| Liu et al. (2004), Gonzalez-Hernandez et al. (2009), Reszka et al. (2007), Zhang et al. (2009), Shi et al. (2016) | not in the atlas | added with new citation keys |

Kept unchanged after checking: Shi et al. 2015 (Snn7), Gao et al. 2015 (Snn6), Friesen et al. 2012 (Snn5).

## Part 2. Röder et al. 1998 primers used across the atlas

Sources: Röder et al. (1998) Genetics 149:2007-2023 (Appendix, 229 GWM primer sets, 279 loci) and Röder et al. (1998) Genome 41:278-283 (Table 2, group 2 deletion-bin markers; it adds 2A for GWM71). Primers are in `scripts/curation/sb19/roder1998_primers.json`; the exact-match in-silico PCR hits on RefSeq v1.0 are in `ispcr_roder.json`.

- 149 of the 229 primer sets amplify on every chromosome Röder mapped them to, 19 on some, 61 on none (exact matches fail where Chinese Spring differs from the mapping parents).
- For 8 markers the Röder sequence differs by one base from the primer used before; the Röder version gives the product for gwm140, gwm153, gwm294, gwm295 and gwm361.
- Rule: a GWM product is used only on a chromosome that Röder et al. mapped that primer set to. A Röder-mapped homoeologous locus without an exact product (for example gwm165 on 4B) keeps its marker-table position instead of being rejected.
- The resolvers for both batches now read these primers. The 19 spot blotch papers gave the same positions as before. In the SNB batch the only change was that QSng.sfr-4BL keeps its gwm165 position.
- QSb.bhu-3B and QSb.bhu-7B (Chirya 3) stay unplaced. Their markers (gwm533, gwm263, gwm1043) are not in Röder 1998.

### Atlas-wide check

Every atlas record with a GWM marker was checked against the positions of its own markers. Meta-QTL projections (Liu 2020, Sharma 2024, Z_1228, Pundir 2023) were left out. A record was corrected when its position lay more than 5 Mb outside the span of its placed markers, or more than 20 Mb from its only placed marker. That gave 35 corrections, all in the Tong et al. (2024) rust compilation. Its v2.1 anchor points often sat on the wrong arm (for example QLr.hebau-2DS at 383 Mb, whereas Xgwm261 and Xgwm484 lie at 19.6 and 48.2 Mb). Each corrected record has a new id (R98-001 to R98-035), its note names the record it replaces, and the old record was deleted. Yr46 was placed at the cloned Lr67/Yr46 gene (TraesCS4D02G364400, 509.28 Mb).

| New id | Locus | Chr | Old (v1.0 Mb) | New (v1.0 Mb) | Grade |
|---|---|---|---|---|---|
| R98-001 | QLr.B22-1A | 1A | 0.000 | 570.05 | C |
| R98-002 | QLr.hbu-2BS.1 | 2B | 369.446 | 88.83-95.34 | B |
| R98-003 | QLr.hebau-2DS | 2D | 383.283 | 19.62-48.17 | B |
| R98-004 | Lr66 | 3A | 164.016 | 214.14 | C |
| R98-005 | QLr.cim-3DC | 3D | 203.617 | 109.55 | C |
| R98-006 | (unnamed, Schnurbusch et al) | 4B | 115.169 | 59.61 | C |
| R98-007 | (unnamed, William et al) | 4B | 235.223 | 59.61 | C |
| R98-008 | QLr.hebau-4B | 4B | 456.525 | 544.65-568.56 | B |
| R98-009 | Lr25 | 4B | 456.525 | 568.56 | C |
| R98-010 | QLr.ubo-7B.1 | 7B | 66.084 | 97.94 | C |
| R98-011 | (unnamed, Schnurbusch et al) | 7B | 131.896 | 97.94 | C |
| R98-012 | (unnamed, Christiansen et al) | 2A | 136.766 | 2.32-31.83 | B |
| R98-013 | QYr.uga-2AS | 2A | 246.770 | 28.20 | C |
| R98-014 | QYr.inra-2AL | 2A | 686.103 | 716.55 | C |
| R98-015 | QYr.ucw-2BS | 2B | 190.765 | 73.62 | C |
| R98-016 | QYrlu.cau-2BS2 | 2B | 312.608 | 100.85 | C |
| R98-017 | QYr.caas-2BS.1 | 2B | 408.946 | 100.85 | C |
| R98-018 | Yr41 | 2B | 652.880 | 95.34 | C |
| R98-019 | QYr.caas-2DS | 2D | 364.089 | 12.36-19.62 | B |
| R98-020 | QYr.cim-3DC | 3D | 203.617 | 109.55 | C |
| R98-021 | QYrid.ui-4A | 4A | 413.373 | 51.14 | C |
| R98-022 | Yr62 | 4B | 253.688 | 568.56 | C |
| R98-023 | QYr.ufs-4B | 4B | 448.266 | 482.82-509.04 | B |
| R98-024 | QYr.caas-4BL | 4B | 478.995 | 509.04-544.65 | B |
| R98-025 | QYrus.vt-4BL | 4B | 478.995 | 509.04-544.65 | B |
| R98-026 | QYrhm.nwafu-4B | 4B | 493.934 | 531.34-568.56 | B |
| R98-027 | Yr46 | 4D | 231.652 | 509.28-509.28 | B |
| R98-028 | Yr-tem-5B.1 | 5B | 398.893 | 315.53 | C |
| R98-029 | QYr.inra-5BL.1 | 5B | 416.968 | 504.30 | C |
| R98-030 | QYrco.wpg-5B | 5B | 416.968 | 504.30 | C |
| R98-031 | (unnamed, Yuan et al) | 6A | 530.501 | 579.13 | C |
| R98-032 | (unnamed, Christiansen et al) | 6B | 300.720 | 89.00 | C |
| R98-033 | QYrpi.vt-6BL | 6B | 344.798 | 430.08 | C |
| R98-034 | QYrst.wgp.6BS.1 | 6B | 716.144 | 44.50 | C |
| R98-035 | QYr.sgi-7A.2 | 7A | 24.440 | 681.39 | C |


## Part 3. Five further SNB sources (11 Nov 2026, `scripts/curation/sb31.py`)

| Source | In the atlas before | Now |
|---|---|---|
| Cockram et al. (2015) G3, file Cockram2010.pdf | absent | Snn1 MAGIC QTL, 1B 0.17-2.10 Mb, peak Excalibur_c21898_1423 at 1.42 Mb (grade B); QSnn.niab-5A.1, 5A 679.49-702.46 Mb, peak BS00068108_51 (grade B) |
| Reddy et al. (2008) Plant Genome | absent | Snn1 high-resolution map (0.46 cM, Xksud14.2-XBF293222), placed at the cloned Snn1 gene (grade B) |
| Virdi et al. (2016) G3 | absent | 7 loci: QTs.fcu-2B at Tsc2 (grade C), QTs.fcu-4B (Xiwa2126-Xgwm113, wide), QTs.fcu-5B at Xgwm234, QTs.fcu-6B (Xiwa5148-Xiwa967), QSnb.fcu-5B at Tsn1 in both populations, LD5B QTs.fcu-5B unplaced (RFLP flanks) |
| Running (2022) PhD dissertation, front matter and chapter 2 only | Snn5-B2 placed on 4B "by analogy" (V1644) | Snn5-B2 corrected to 2BS (the abstract says so); left unplaced because chapter 5 with the markers was not supplied |
| Saini et al. (2022) Mol Breed | meta-QTL source, no records | not added: the atlas does not keep meta-QTL |

The Xfcp618 primers in Reddy et al. 2008 give no exact product on RefSeq v1.0 (the marker is dominant for the T. dicoccoides allele); Xpsp3000 amplifies at 1B 5.14 Mb.

## Part 4. Septoria tritici blotch (STB) batch, Nov 2026

Script: scripts/curation/stb44.py (records in stb44_recs.py; Riaz/Teagasc Table 5 in stb44_riaz_table5.json).

- 174 new records S44-001 to S44-174 (grade A 62, B 28, C 45, D 39) replace 93 older STB records from 22 sources, which were deleted.
- 19 sources added or re-keyed with corrected citations (Arraiano2007PP, Ghaffary2011TAG, Ghaffary2012TAG, Liu2013Eu, Ferjaoui2022BMC, LanglandsPerry2022 and 13 others).
- Main corrections:
  - Teagasc MAGIC (Riaz et al.) records had cM stored as Mb; now the printed Mbp intervals, highest-LOD dataset per QTL.
  - Langlands-Perry 2022: 21 per-gene records merged into the 3 reported QTL.
  - Chartrain 2005 "Stb12" corrected to Stb7 (probable).
  - QStb.wai-3D moved from 3DS (Xgwm341) to 3DL (590.7-603.5 Mb).
  - Stb9 placed at Xwmc317 (2BL, about 784 Mb); Stb3 placed on 7AS.
- New in-silico PCR hits: gdm14 (6D 283.60 Mb), barc142 (5B 605.52 Mb).
- Skipped: duplicate PDFs, theses that repeat the papers (Karlstedt, Cobertera, Harilal, Seed, Ghaffary), Saini 2022 (meta-QTL), Brown 2015 (review).
- Not changed (outside this batch): A_STB_033 lists Stb17 on 3B; Stb17 is on 5AL.

## Part 5. Tan spot batch (Archive_1.zip), Nov 2026

Script: scripts/curation/ts52.py (records in ts52_recs.py). New in-silico PCR hits for 36 SSRs and three Tsc2 EST markers were added to sb19/ispcr_res3.json.

- 116 new records S52-001 to S52-116 (grade A 15, B 31, C 44, D 26) replace 61 older records from 12 sources (Chu 2008, Chu 2010, Dinglasan 2021, Faris and Friesen 2005, Faris 2012, Kalia 2018, Kariyawasam 2016 and 2018, Li 2011, Liu 2017, Shankar 2017, Zou 2017).
- New sources: Abeysekara 2010 (Tsc2), Effertz 2002 (Tsc1), Corsi 2020 (Tsc2 MAGIC QTL, printed v1.0 positions), Faris 2020 (Tsr7), Cheong 2004, Liu 2006 (SNB at Tsn1), Liu 2020 (race 2, durum; printed Mbp intervals), Singh 2008 (tsn2 and tsn5), Sun 2010. Zwart 2010 yellow spot QTL added to the existing Zwart source. Chu 2010 and Li 2011 now carry the correct year.
- Main corrections:
  - Dinglasan 2021: the old records used DArTseq marker IDs as Mb (for example 3033029 read as 3.03 Mb). Now placed at the gene models named in the paper (1A 4.66, 5B 339.92, 5D 152.75, 6A 609.25 Mb) and the printed 2B span (137.81-215.93 Mb).
  - Faris and Friesen 2005 QTs.fcu-3BL placed at Tsr7 (3B 466.3-470.0 Mb), following the re-analysis by Faris et al. (2020); the old position used Xbarc248 (162 Mb).
  - Tsc1-related QTL (Kariyawasam 2018, Liu 2017, Effertz 2002) placed at the fine-mapped Tsc1 interval; Tsn1-related QTL at the Tsn1 gene.
  - Shankar 2017: 23 QTL from both populations (old atlas had 2 unplaced); 2A QTL set by wPt-9320 and IWA5893 because the other interval markers place inconsistently.
  - Faris 2012 QTs.fcu-5D set at Xwmc215 because Xcfd183 amplifies only on 5DS.
- Two wide-gap liftovers fixed with 90K marker pairs (as in Part 1).
- Skipped: Chu 2007 (screening only), See 2023 (race survey), Liu 2020 meta-QTL and Saini 2021/2022 (meta-QTL), Stadlmeier 2019 and thesis, Hu 2019, Singh 2019, Liu 2004 and Virdi 2016 (already done), Bemister 2019 (re-read in batch 1), duplicate PDFs.
- Open: tsn2/tsn5 rest on Xgwm285 (3B 273 Mb), which lies on proximal 3BS although the genes are reported on 3BL; Zwart QYls.lrc-5A rests on Xgwm443 (5AS) although the paper relates it to 5AL.

## Part 6. Meta-QTL studies removed (user request, Nov 2026)

Script: scripts/curation/mq53.py. 497 records removed, all from meta-QTL studies:
- Liu et al. 2020 (tan spot meta-QTL, 215), Sharma et al. 2024 (powdery mildew meta-QTL, 278), Mao et al. 2010 (FHB QTL meta-analysis, 2);
- from Pundir et al. 2023, its 2 meta-QTL only. Its own 111 QTL (HUW 468 x C 306) and the 7 earlier CCN loci it anchored physically (QCcn-1B, QCre.srd-1B, QCre-ma2A, QCre.pau-2A, Cre5, QCcn-3A, QCcn-2A) are kept, as the original studies are not yet in the atlas.
A first pass also removed 5 records from studies that are not meta-QTL studies (Lhamo et al. 2024, 3; Arruda et al. 2016, 1; Goudemand et al. 2013, 1). These were restored the same day.

## Part 7. Spot blotch, black point and leaf spot GWAS batch (GWAS_SB_1, GWAS_BPR, GWAS_LS), Nov 2026

Script: scripts/curation/gw60.py. 19 unique papers checked against the atlas (Gurung 2014 and Bainsla 2020 were in two folders).
- New papers (118 records S60-001 to S60-118, all grade A, positions printed in the papers): Bainsla et al. 2020 (25 MTAs, FarmCPU and MLMM), Aggarwal et al. 2026 (7 stable MTAs), Lv et al. 2020 (27 black point SNPs and 3 linkage QTL), Tang et al. 2022 (11 black point QTL), Iqbal et al. 2022 (37 QTL for common bunt, leaf rust, leaf spot and stripe rust, RefSeq v2.0).
- Chen et al. 2025 (black point): the 5 auto-extracted records (one at 3A 0.87 Mb) replaced by the 6 QTL of Table 2.
- Juliana et al. 2022: 5B_68806 was a line-break truncation of 5B_688063050 (688.06 Mb, not 0.07 Mb); 5B_687326022 added.
- Grades corrected without moving positions: Juliana 2022 and Tomar 2021 GBS markers (position in marker name), Kamble 2023, Lozano-Ramirez 2022, Thirunarayanan 2024 and two Liu 2017 records raised to A (position printed); Gurung 2014 (25) and two Ayana 2018 records lowered to B (paper gives cM only).
- Thirunarayanan et al. 2024: citation corrected (it read "P et al.").
- Checked and unchanged: Li 2020, Liu 2017, Jamil 2018, Singh 2023, Kumar 2022, Tomar 2021 QTL, Ahirwar 2018, Singh 2020.
- Black point records stay under "Spot blotch" with "(black point)" in the name, as for Li 2020 and Liu 2017.
- Papers reviewed: 890.

## Part 8. GWAS batch 7 (STB, SNB, tan spot folders), Nov 2026 (gw70.py)
- 52 PDFs checked (37 distinct papers after duplicates). 21 were already in the atlas and were left unchanged.
- 15 papers added (202 records, S70-001 to S70-202; grades A 136, B 40, C 4, D 22): Alemu 2021, Aouini 2025, Chaudhary 2025 (RefSeq v2.1), Kassie 2025, Kidane 2017, Mahboubi 2022, Mikaberidze 2025 (resistance MTA only), Odilbekov 2019, Yates 2019, Miedaner 2013 (unplaced), Seed 2023 thesis, Liu 2015, Phan 2018, Berisha 2024 (2 MTAs from text; RefSeq v2.0), Mandal 2024.
- Mandal 2024 positions match the RefSeq v1.0 array SNP positions and are read as v1.0.
- Reviewed, no loci added: Thauvin 2024 (positions only in Table S5), Juliana 2017 (genomic prediction only), Gao 2016 (pathogen GWAS), Altameemi 2021 thesis (same data as AlTameemi 2021), Wang 2025 preprint (published version already in atlas).
- Papers reviewed: 890 -> 909.

## Part 9. Position fixes from re-uploaded papers, Nov 2026
- Manan et al. 2025 TAG (bacterial leaf streak; mn25.py, S80-001 to S80-009): all positions in the paper are on the durum Svevo v1 assembly, not Chinese Spring. The 5A MTAs moved from 657 to 698 Mb (RefSeq v1.0); 6A MTA to 1.50 Mb; QBls.zhl-6A.1 now the CIM interval shifted to CS; 6A.2 and 3A use Svevo coordinates as an approximation (grade C); QBls.zhl-5A start corrected (was a v2.1 value). Duplicate per-trait MTAs merged; V3250-V3260 removed.
- Pundir et al. 2023 BMC Plant Biol (cereal cyst nematode; pu23.py, S81-001 to S81-009): main-effect QTL now carry the closest-marker position from Table 2 as peak (grade A), with the Table S3 flanking interval; CIM and ICIM positions for 2A.1/2A.2 noted. Epistatic QCcnr.ccsu-2A.5 placed from the Table S3 positions of WSNP6862 and WSNP4934. The other ~100 epistatic loci have only cM positions and study-internal WSNP codes and stay unplaced. V3162-V3169 and V3332 removed.
- Pundir et al. 2022 Sci Rep (cereal cyst nematode, ITMI DH; pu22.py, S82-001 to S82-004): QCcn.ha-3B, -5D, -2B.1 and -1D placed at the peak-marker positions printed in Table S8 (grade A); the other 14 QTL regraded B -> A (same positions); QCcn.ha-4A stays unplaced (no BLAST hit in the paper).

## Part 10. Additions from four uploaded papers (ad4.py, Nov 2026)

- Pundir et al. (2022): 85 digenic epistatic loci for cereal cyst nematode from ESM Table S3 (36 grade B, 42 grade C, 7 grade D). Loci whose flanking markers lie more than 100 Mb apart are grade C.
- Manan et al. (2025): 14 bacterial leaf streak MTAs listed only in Table S5 (12 B, 2 D), plus QBls.zhl-6A in the BP025 population (B). Positions are Chinese Spring 90K positions, not the Svevo ones.
- Muellner et al. (2021): 10 single-population common and dwarf bunt intervals (MP-BLI, MP-BON) from Tables 2 and 3 (grade A).
- Pundir et al. (2023): nothing further to add.
- IDs S83-001 to S83-110.

- ad5.py: V1600 (auto-extracted "QTL_1D", Arora et al. 2026 preprint) corrected to Pch4, cloned NLR for eyespot (O. yallundae) on 1DS. Kept grade D: the paper gives only Ae. tauschii TA1675 coordinates (6.23-6.93 Mb).

- ad6.py: Martin et al. (2025) Crop Science, Fusarium crown rot. 9 auto records (V3423-V3431; 5 had no position, none lifted to v2.1) rebuilt from Table 4 and Table S3 printed CS positions as S84-001 to S84-018: 7 FCR QTL plus 4 white-head QTL (Table 4) and 5 FCR plus 2 white-head single-environment QTL (Table S3). 8 grade A, 10 grade C (interval > 100 Mb). Yield QTL not added.

- ad7.py: Hou et al. (2023) TAG 136:101, FCR GWAS. Positions are on RefSeq v2.0 (Methods), but the atlas stored them as v1.0/v1.x and unlifted (1B-4 shifted about 7 Mb on v2.1); 2 loci had no position and several auto records took another locus position. 8 records (R049, V1356-V1362) replaced by all 21 Table 4 loci, S85-001 to S85-021, grade A; KASP1869 primers added to 1B-4.

- ad8.py: Lin et al. (2023) Plant Disease 107:1151, FCR GWAS. Positions on RefSeq v2.0 (Methods); atlas had 4 of the 10 stable loci as v1.x (0.5-8 Mb off on v2.1) and two auto records (V1349 unplaced, V1350 7D at the 7B position). R117-R119, R154, V1349, V1350 replaced by the 10 Table 2 loci, S86-001 to S86-010, grade A, peak SNP from Table S4, KASP primers (Table S6) on 2A-1 and 7D.

- ad9.py: R158 Ra_c9427_300 (7A, adult FCR) moved from source key Rahman2020FCR (which carries the Lin 2023 citation) to new source Rahman2020MB (Rahman et al. 2020 Mol Breed 40:28, reference as given in Lin et al. 2023; DOI checked). Position: marker RefSeq v1.0 78.434 Mb, v2.1 82.085 Mb, matching the 82.09 Mb quoted by Lin et al. New ID S87-001, grade C (single marker).

- ad10.py: Li et al. (2010) TAG 121:941, Batavia/Ernie DH, crown rot and FHB. Auto records V1407 (3B, unplaced) and V1408 (2D; the 2-49 locus cited in the Discussion, not a result of this paper) replaced by the 4 Table 4 loci, S88-001 to S88-004: 3BL crown rot QTL (B, flanks wPt-9368 and wPt-1834 from the Fig 1-2 map order), FHB 3A and 6A (C, single DArT marker), FHB 5B (D, wPt-3661 unplaced). Plant-height QTL not added.

- ad11.py: Erginbas-Orakci et al. (2018) IJMS 19:2666, crown rot (F. culmorum), CIMMYT spring panel. Table 3 has 3 MTAs; wPt-2193 and wPt-2766 are coincident (3B 39.1 cM), so V0454 (wPt-2193, unplaced) merged into V0455. V0455 and V0456 regraded B to C (single DArT marker), stats and notes from the paper; positions unchanged.

- ad12.py: Li et al. (2024) TAG 137:43, FCR GWAS, 223 accessions. 35 per-SNP auto records (V0959-V0991, V1309, V1310), many labelled RefSeq v2.1 with wrong values and missing 3D-1/3D-2 positions, replaced by the 10 Table 4 QTL (S89-001 to S89-010, grade A, RefSeq v1.1 = v1.0 coordinates). KASP3D61373 and RIL validation on Qfcr.cau.3D-3. Candidate gene Z0256 (TraesCS3D02G524700) kept.

- ad13.py: Liu et al. (2011) Czech J Genet Plant Breed 47:S109 (summary paper). V1585 (auto QTL_3B, unplaced) replaced by S90-001 Qcrs.cpi-3B (CSCR6), grade C at gwm0181 (2.2 cM from the QTL; flanking DArT markers unplaced). Ernie 3BL QTL not duplicated (primary record S88-001, Li et al. 2010); barley 3HL QTL outside the wheat atlas.

- ad14.py: Alahmad et al. (2020) IJMS 21:5260, durum crown rot GWAS. V1627 (auto QTL_6B, unplaced) replaced by S91-001 qCR-6B (with co-located qWH1-6B, qWH2-6B), grade B, interval from DArTseq tags placed by BLAST (6B 667.20-674.33 Mb, v1.0); paper gives only Svevo cM. Plant height, stay-green and yield QTL not added.

- ad15.py: Neupane et al. (2024) BMC Plant Biol 24:852, FHB/FDK/DON, DH D8006W/Superior. 30 records (Z0191-Z0194, V0285-V0310; 25 unplaced, 2 on the wrong chromosome) replaced by 29 disease QTL from Table 3 and Table S3 (S92-001 to S92-029). CI of the most significant row converted to RefSeq v1.0 with the Table S5 map-marker positions; grade B (26) or C (3, span > 100 Mb). PHT, DA and AE QTL not added.

- ad16.py: Thambugala et al. (2020) TAG 133:2775, native FHB in Cutler/AC Barrie and AC Barrie/Reeder. 21 unplaced auto records (V0332-V0352) replaced by the 21 FHB QTL of Tables 5-6 (S93-001 to S93-021), grade B. Paper gives cM only; placed from the population linkage maps (Tables S5, S11) via the RefSeq v1.0 positions of 90K map markers around the strongest peaks. Anthesis and height QTL not added.

- ad17.py: in-silico sweep of the 762 records with no v2.1 position. 261 name a marker; 10 resolve on their own chromosome; 7 of those had been left unplaced on purpose (arm conflict or paralogous flanks) and stay so; the 3 Navathe et al. (2023) DArTseq markers placed (S94-001 to S94-003, grade C). Resolver normalisation kept in scripts/curation/lib/insilico_res.py.

- ad18.py: Lu et al. (2013) TAG 126:317 (cited as 2012), FHB in SHA3/CBRD x Naxos. 6 auto records (V1421-V1426) replaced by the 24 QTL of Tables 4-5 (S95-001 to S95-024): flanks placed in silico; 6 B, 13 C, 5 D (no placeable flank).

- ad19.py: He et al. (2016) PLoS ONE 11:e0158052, FHB in Soru#1 x Naxos. 7 unplaced auto records (V1552-V1558) replaced by the 9 disease QTL of Tables 4-6 (S96-001 to S96-009): flanks placed in silico after expanding the paper marker abbreviations; Rht-D1, Vrn-A1 at gene positions. 4 B, 4 C, 1 D (2DL flanks map to 6D in CS). KASP primers for 2DLc added.

- ad20.py: Nannuru et al. (2022) TAG 135:2247, FHB, NMBU Nordic spring wheat GWAS. Table 3 prints 13 QTL regions on RefSeq v1.0: 8 existing records regraded A with notes (R116, R150, V1341-V1344, V1346, V1347); V1345 (6A.1 with the 6B.1 interval) replaced by S97-004 (450-531 Mb); 1A.1, 3A.1, 3B.1, 6A.2 added (S97-001 to S97-005). KASP primers (S14) on 7A.2.

- ad21.py: Steiner et al. (2004) TAG 109:215, FHB in Frontana x Remus. 7 unplaced auto records (V1383-V1389) replaced by the 12 FHB QTL of Tables 4 and 6 (S98-001 to S98-012): SSRs placed in silico (3A at Xdupw227, flagged: paper says pericentromeric; 5A Xgwm129-Xbarc197; 6A Xgwm169), grade C; 9 AFLP-only QTL grade D.

## ad22 – Jia et al. (2005) Euphytica 146:183, FHB, Wangshuibai x Alondra's DH
- Deleted 7 unplaced auto records V1457–V1463 (1D and 7D were single-marker hits not significant in CIM; not kept).
- Added S99-001 to S99-014: 13 CIM QTL from Table 2 (same interval across years merged) and the epistasis-only 6A QTL (Table 3). SSRs placed in silico; grades B (4), C (9), D (1: 5B Xgwm443-Xbarc32, neither marker placed). 3BS placed at Xgwm493 only (Xgwm533 not resolved).

## ad23 – six uploaded papers (Nov 2026)
- Ballini et al. (2020) Euphytica 216:92, durum EPO GWAS: S100-001 to S100-016 (SNB, wheat blast, STB, leaf rust). Zavitan v2.0 positions; SNB_Tox1B placed at Snn1, WB_BR43_2A at Rmg7, WB_BR43_2B at Rmg8, STB_Pyc3A at TraesCS3A01G041800, STB_Nec3A 3.3 Mb distal of it; the rest use the Zavitan position as an approximate RefSeq v2.1 position (grade C).
- Arraiano et al. (2007) TAG 116:113: S101-001 Hobbit sib 5BS susceptibility region; S101-002 Stb6 in Bezostaya 1.
- Agostinelli et al. (2012) Euphytica 186:115: S102-001 Fhb1 and S102-002 QFhs.nau-2DL effects in the VA01W-476 cross.
- Arruda et al. (2016) Mol Breed 36:84 reprints the GWAS MTAs of Arruda et al. 2016 Plant Genome: S103-001 to S103-009 added to that study (A_FHB_024); Fhb1-LD SNPs at Fhb1, others grade D (CSS scaffold names, cM only).
- Not added: An et al. (2019) TAG (rye 4R addition) and Anderson et al. (2010) Euphytica (Lophopyrum substitutions): alien chromosomes, no wheat-genome position.

## ad24 – thirteen uploaded papers (Nov 2026)
- Added S104-001 (Ban & Suenaga 2000, B1-linked Sumai 3 gene, 5AL), S105-001 to S105-011 (Buerstmayr et al. 1997 monosomic analysis, chromosome-level, grade D), S106-001/002 (Bertucci et al. 2014, Tsn1 and Snn3-B1 survey), S107-001/002 (Bokore et al. 2017, Fhb1 and Fhb5 validation), S108-001/002 (Anh et al. 2015, Rmg8 and Rmg7), S109-001/002 (Burt & Nicholson 2011, Pch1 and Pch2), S110-001/002 (Ceoloni et al. 2017, Fhb-7EL alien, grade D), S111-001 (Bernardo et al. 2012, Fhb1 SNP map), S112-001 (Cao et al. 2001, snbTM 3A durum, grade D).
- Not added: Ban 2000 Breed Sci (RAPD QTL on unassigned linkage group), Ban & Suenaga 1997 (same B1-linked gene; RAPDs unplaced), Ansari et al. 2007 (transcript study), Belhadj Benyahia et al. 2020 (biocontrol).

## ad25 – eight uploaded papers (Nov 2026)
- Added S113-001 (Chaubey et al. 2023 preprint, Tsn1 spot blotch), S114-001 to 003 (Choupannejad et al. 2023 preprint, Tsn1/Snn1/Snn3-B1), S115-001 (Chartrain et al. 2005, Stb6 survey), S116-001 (Chapman et al. 2009, Pch2 candidates), S117-001 (Chen et al. 2007, Qfhs.ndsu-3AS at Xgwm2).
- Not added: Chu et al. 2006 (hybrid necrosis Ne1/Ne2), Mol Breed 2017 37:141 (agronomic GWAS), Trop Plant Pathol 2025 (cultivar screening, no locus).

## ad26 – seven uploads (Nov 2026)
- Added S118-001 to 006 (Cruppe et al. 2023, wheat blast QTL incl. 2NvS), S119-001 (Esvelt Klos et al. 2014, Pch2), S120-001 to 011 (Dinglasan et al. 2019, yellow spot GWAS, DArTseq placed by BLAST), S121-001 (Dutta et al. 2021, Stb6 SNP at 26.46 Mb, grade A).
- Not added: Duguid & Brule-Babel 2001 (inheritance only), Du et al. 2015 (alien 2Ns(2D) substitution).

## ad27 – eight uploads (Nov 2026)
- Added S122-001 to 004 (Friesen & Faris 2004: tsc2 at Tsc2 position, three RFLP-only minor QTL grade D), S123-001 to 005 (Faris et al. 1997: 1AS QTL at Tsc1 position, four RFLP-only QTL grade D), S124-001 to 007 (Faris et al. 1999 candidate genes: Lr34 region at Lr34, others grade D), S125-001 (Forte et al. 2014, Fhb-7el2 alien, grade D), S126-001 (Gholizadeh Vazvani et al. 2025, take-all, Xbarc197; new disease Take-all).
- Not added: Ferjaoui et al. 2021 preprint (already in atlas from BMC Genomics 2022), Gamba & Lamari 1998 (inheritance only), Faris et al. 2014 MGG (spike compactness).

## ad28 – ten uploads (Nov 2026)
- Added S127-001 to 006 (Grausgruber et al. 1997, substitution-line chromosome effects, grade D), S128-001/002 (Gyawali et al. 2023, Fhb1 and Fhb7), S129-001 (Guo et al. 2003, Fhb1 AFLP/STS), S130-001 (Hayden et al. 2004, Sr2 STMs), S131-001 to 003 (Gilsinger et al. 2005, FHB incidence/flower opening), S132-001 (Guo et al. 2015, Fhb7, alien), S133-001 to 003 (Herter et al. 2018, coded markers, grade D), S134-001 to 004 (Herter et al. 2019: Rht-D1 for FHB and STB, Fhb5, 4B), S135-001 to 003 (He et al. 2014 haplotyping: Wuhan 1 and CJ 9306 2D QTL, Fhb1).
- Not added: Hiebert et al. 2010 (telocentric gene stacking, no mapping).

## ad29 – ten uploads (Nov 2026)
- Added S136-001/002 (Jiang et al. 2015: Rht-D1 and Ppd-D1 associations with FHB, at the genes), S137-001 (Islam et al. 2024: Rmg11 on 7AS, Svevo position used as approximate v2.1), S138-001 to 006 (Zhang et al. 2025 preprint: leaf rust APR QTL; 1BL.2 and 3A unplaced), S139-001 (Jin et al. 2020: crown rot 5DL region 218.6-232.4 Mb, grade A), S140-001 (Kage et al. 2018: QTL-2DL candidate genes, wide region).
- S134-001 note corrected: Rht-D1 position confirmed from the gene table (TraesCS4D02G040400, 18.781 Mb on v1.0).
- Not added: Joshi & Chand 2002, Joshi et al. 2007 (no mapping), J Bangladesh Agric Univ 2019 (in silico survey), TAG 2017 130:471 (no loci), Hu et al. 2019 correction (copyright only; original already in atlas).

## ad30 – Kwiatek et al. 2015 (Nov 2026)
- Added S141-001 Pch1 survey of 372 breeding lines (EpD1b, XustSSR2001-7DL), from the abstract; placed at the Pch1 interval.

## ad31 – five uploads (Nov 2026)
- Added S142-001 to 004 (Kharina et al. 2025: Stb2, Stb9, Stb11, Stb13 survey at diagnostic SSRs), S143-001 (Kwiatek et al. 2012: Pch1 survey), S144-001 (Leonard et al. 2008: Pch1/Ep-D1 STS markers, placed across Xwmc14/Xbarc97/Xcfd175, grade B).
- Not added: Lendenmann et al. 2016 (pathogen QTL), Khodarahmi et al. 2020 (inheritance only).

## ad32 – Liu et al. 2008 (Nov 2026)
- Added S145-001 Fhb1 261 kb BAC interval and UMN10, from the abstract; placed at the cloned Fhb1 position.

## ad33 – eight uploads (Nov 2026)
- Added S146-001 SbWE35 (Li et al. 2026; 4BL 581.16-582.11 Mb from interval genes, grade A) and S146-002 Sb3 (WE35); S147-001/002 Qcrs.cpi-3B and Qcrs.cpi-4B (Ma et al. 2010); S148-001 Qcrs.cpi-3B NILs (Ma et al. 2012).
- Not added: Lorrain et al. preprint and Nature Plants 2026 (pathogen), Mackintosh et al. 2006 (transgenics), Ma et al. 2015 MGG (MDA QTL), Loffler et al. 2009 (meta-QTL, excluded by rule).

## ad34 (Nov 2026)
- Added S149-001 Stb7 (McCartney et al. 2003); S150-001 to S150-014 crown rot QTL (Martin et al. 2015; arm check discards markers placed on the wrong arm); S151-001 to S151-005 FHB alleles in Canadian backgrounds (McCartney et al. 2007); S152-001 to S152-003 stacked FHB QTL (Miedaner et al. 2006); S153-001 Pch1 markers (Meyer et al. 2011); S154-001 Fhb1 multiplex PCR (Matus-Cadiz et al. 2006).
- Not added: Mahboubi et al. preprint (published version already in atlas); McCartney et al. 2004 (haplotype survey, no new loci).

## ad35 (Nov 2026)
- Added S155-001 to S155-003 (Puskas et al. 2006; 3BS, 2BS, 6A FHB from Ning8331; from page images, Table 1 pages 629-630) and S156-001 to S156-003 (Kalih et al. 2014 triticale; 5A grade D, rye 5RL Ddw1 and 6R as alien grade D).

## ad36 (Nov 2026)
- Added S157 (Rahman et al. 2021, 3 crown rot resistance QTL), S158 (Nishio et al. 2008), S159 (Rawat et al. 2016, Fhb1 PFT), S160 (Mohler et al. 2002, 1DS FHB), S161 (Rosyara et al. 2009), S162 (Salameh et al. 2011), S163 (Miedaner et al. 2011, 9 associations incl. Rht-D1).
- Not added: Shen et al. 2004 (substitution lines), Mirdita et al. 2015 (genomic prediction, no loci), Running et al. 2025 (already in atlas), Muhovski et al. 2014 (expression, arm assignment only); Rahman 2021 tolerance QTL.

## ad37 / ad37b (Nov 2026)
- Added S164 (Soresi et al. 2017, Qfhs.ndsu-3AS), S165 (Simon et al. 2004, 5 STB QTL, RFLP, grade D), S166 (Simon et al. 2007, Stb5), S167 (Somers et al. 2005, 6 FHB QTL), S168 (Srinivasachary et al. 2008 TAG 116, 9 FHB QTL incl. Rht-D1), S169/S170 (Sheng et al. 2012, 2014, Ae. longissima eyespot QTL, alien D), S171 (Su et al. 2018), S172 (Slominska-Durdasiak et al. 2020, D), S173 (Tadesse et al. 2008, D), S174 (Shen and Ohm 2007, alien D), S175 (Singh et al. 2010, Tsr1, Tsr6).
- ad37b: erratum Euphytica 200:27; V1478-V1483 repointed to SzaboHever2014Euph (author names corrected).
- Not added: Song et al. TaNDC1 preprint and TaU11/U12-35K (functional gene studies), Stewart et al. 2016 (method), Steed et al. 2002 (conference version of Steed et al. 2005), Shi et al. 2008 (haplotype survey).

## ad38 (Nov 2026)
- Added S176 (Tessmann and Van Sanford 2018, 12 MTAs, D), S177 (Tessmann et al. 2019, 7 MTAs, D), S178 (Tadesse et al. 2010, tsr4), S179 (Waldron et al. 1999, 4 QTL), S180 (Tidakbi et al. 2026, 20 loci on v2.1 incl. Jamestown 1BL haplotype), S181 (Yoshioka et al. 2024, Rmg10), S182 (Tagle et al. 2015, Rmg7), S183 (Wang et al. 2024, Fhb8), S184 (Zhang et al. 2011, FhbLoP, alien D), S185 (Zanke et al. 2017, Pch1, 5A, 2AS).
- Not added: Wang et al. preprint (Elymus 1StL telosome addition), Zheng et al. 2021 (meta-QTL), Cardozo Tellez et al. 2022 (screening), Zhong et al. 2021 (pathogen), Wilhelm et al. 2013, Yao et al. 2019, Zhang et al. 2022 (already in), Zhao et al. 2018 (functional).

## ad39 (Nov 2026)
- Added S186 (Yang et al. 2010, Qcrs.cpi-4B soil-less assay), S187 (Wilde et al. 2007, Fhb1, Fhb5, Frontana 3A), S188 (Timonova et al. 2013, QLr.icg-5B/2A/1A, T. timopheevii, D), S189 (Xiao et al. 2016, Qfhs.pur-7EL and Lr19, alien D), S184-002 (Zhang et al. 2011, Lr19, alien D), S190 (Zhang et al. 2016, Fhb1 and Sr2), S191 (Zheng et al. 2017, Qcrs.cpi-3B/5D/2D), S192 (Wang et al. 2025, PFT not required for Fhb1), S193 (Yu et al. 2008, SSCP markers at Fhb1), S194 (Zhang et al. 2021, Snn3-D1 cloned).
- Not added: Walter et al. 2008 (expression), Yang et al. 2006 (SSR survey), Tessmann and Van Sanford 2019 Euphytica (morphology), Soltanloo et al. 2011 (diallel), Sucur et al. 2026 (Stb survey); duplicate PDF of Zhang et al. 2016.

## ad40 (Nov 2026)
- Buerstmayr and Buerstmayr (2015) Capo x Arina: from Online Resource 3 (Fig. S1 linkage maps and LOD profiles) the six auto records V1429-V1434 (five unplaced) were rebuilt with flanking markers (1BS B, 5AL at B1 B, 3B, 6BL, 3D, 7D C) and two QTL added (S195-001 2A, S195-002 4A, C). Main paper not uploaded; LOD and PVE not entered. Removed from the supplementary-needed list.

## ad41 (Nov 2026)
- Ferjaoui et al. (2022) Additional Table 4: Qstb7A, Qstb2B_1 and Qstb2B_2 (S44-023 to S44-025) placed from the Svevo positions of their flanking tags (approximate v2.1, grade C; Qstb2B_1 spans two inconsistent hits 230 Mb apart); Qstb1A note updated; Qstb5B is not in the table and stays D.
- Dinglasan et al. (2019) Online Resource 3 (89 MTAs on the DArTseq consensus map, cM only): no additional DArTseq clone could be placed near the 11 QTL, so records are unchanged; kept for the planned consensus-map conversion.
- Both removed from the supplementary-needed list.

## ad42 (Nov 2026)
- Herter et al. (2018): S133-001 to S133-003 placed by decoding the marker codes with Supplementary Table 2 (ESM2): 5B 587.3-589.1 Mb (B), 2D 636.0 Mb (C), 7D 501.0-517.3 Mb (B).
- Capo x Arina (Online Resource 3, first map): V1433 1BS gains Xgwm11 (216.0 Mb) as proximal anchor; S195-001 2A re-placed from the 2AS end to wPt-1657 (0-47.2 Mb), the map orientation fixed by Xgwm558, Xgwm312 and Xgwm294.
- Dinglasan et al. (2019) Online Resources 2-3: DArTseq IDs and consensus cM only; no new placements.

## ad43 (Nov 2026)
- Jin et al. (2020) Table S6 (ESM9): 104 crown rot SNPs (RefSeq v1.0 positions as printed) grouped into 13 regions; 12 new records S196-001 to S196-012 (grade A), 5DL summarised in S139-001 (35 SNPs, 189.7-246.8 Mb).

## ad44 (Nov 2026)
- Kokhmetova et al. (2021) Supplementary Table S5: 11 race-specific tan spot MTA SNPs (race 1: 8, race 5: 3) added as S197-001 to S197-010 (7B pair merged); 6 placed by DArT tag BLAST (1 B, 5 C), 4 unplaced (D). Main paper not uploaded.

## ad45 (Nov 2026)
- Martin et al. (2015) Supplementary Figure 1 (ESM1 EMF linkage maps): 9 crown rot QTL re-placed from the markers inside or flanking each QTL bar (S150-001, -002, -003, -004, -006, -007, -008, -013, -014; D to C for 2BS, 3BS, 4BS; C to B for 1AS, 1DL, 3BL x2, 6DL); two Janz-derived QTL added (S150-015 2BS, S150-016 2DS).

## ad46 (Nov 2026)
- Pankaj et al. (2022) ESM1 Table S3 primer sequences run by exact-match in-silico PCR on RefSeq v1.0: S19-059 (5B) placed at wmc173 119.165 Mb and S19-068 (1B) at cfd59 205.421 Mb (tentative, 1,464 bp product), both D to C. Primers for wmc109, gwm312, gwm213, gwm608, gwm614, gwm497 gave no product on the mapped chromosome; S19-064, S19-072, S19-073 stay D.
- Rahman et al. (2021) Supplementary Table S3: QTL_8 SNP Excalibur_c65830_82; S157-002 placed at 246.133 Mb (v1.0), D to C.

## ad47 (Nov 2026)
- Singh et al. (2018) IJMS Supplementary Table S1 (QTL with DH and PH as covariates): five covariate-only spot blotch QTL added, S198-001 BC 3AL (D), S198-002 BC 3B 21.86-23.95 Mb (B), S198-003 BC 5B (D), S198-004 BC 7A 448.59 Mb (C), S198-005 WC 2A 71.08 Mb (C). Other Table S1 QTL match existing S19-087 to S19-099 records; unchanged.

## ad48 (Nov 2026)
- Slominska-Durdasiak et al. (2020) Online Resource 2 primers (T3/T7 tails removed) by in-silico PCR: S172-001 Dicer1 placed in TraesCS5A02G516000 (5A 679.66 Mb) and S172-002 Ara6 in TraesCS4D02G267900 (4D 439.12 Mb), both D to C.

## ad49 (Nov 2026)
- Tessmann et al. (2019) Supplementary Table S3 decodes the panel SNP IDs to 90K names. Placed: S177-001 7B 673.96-678.64 Mb (B); S177-003 3A, S177-004 2B, S177-006 4A, S177-007 6B (C). Same IDs in Tessmann and Van Sanford (2018): S176-007 4A, S176-008 6B, S176-010 6B (C). S177-002 (SNP on 5D only) and S177-005 (no position) stay D; the other nine 2018 IDs are not in Table S3 and stay D.

## ad50 (Nov 2026)
- Virdi et al. (2016) Table S1 AL map checked against the AL tan spot QTL. S31-004 Tsc2 confirmed between Xiwa2304 (16.90 Mb) and Xwmc154 (36.45 Mb). S31-005 4B: map gap 9.7-44.2 cM, cannot be narrowed; grade B to C (242 Mb span). S31-006, S31-007, S31-008 unchanged.
