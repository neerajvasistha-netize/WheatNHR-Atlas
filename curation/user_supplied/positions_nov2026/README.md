# Physical positions from ten uploaded papers (Nov 2026)

Papers supplied by the atlas owner on 19 Nov 2026; records written to `web/seed_ad66.json` by `ad66_build.py`
(full record list with status in `records.tsv`). Positions were converted to RefSeq v1.0/v2.1 and the 13 pangenome
cultivars with the genome-wide gene anchors (`scripts/pangenome/convert.py`).

| Paper | Atlas source | What was done |
|---|---|---|
| Xu et al. 2024 TAG 137:49 (crown rot, C549/3642) | A_crown_rot_010 | 9 QTL placed from Table 2 flanking SNPs (named chr_position, v1.0). The 6B/7B SNP names in Table 2 lack one digit; the Discussion intervals (6B 281.6-327.0, 7B.1 315.6-432.5 Mb) were used and the six 6B/7B marker records corrected. |
| ElDoliefy et al. 2020 Euphytica 216:71 (FHB, Glenn) | V_ElDoliefy2020 | 4 QTL placed from the Fig. 2 DArT map markers (DArT tag positions, v1.0). 6BS interval is wide (151.8-466.8 Mb). |
| Farhad et al. 2026 PLOS One 21:e0349201 (wheat blast) | V_Farhad2026 | All 45 QTL of Table 3 with regions and leading SNPs (v2.1): 4 existing records placed, 41 added. |
| Zhao et al. 2023 Mol Breed 43:82 (spike length, FHB) | V_Zhao2023 | QSl.yaas-5A and QSl.yaas-5B added (Table 4, v2.1; positive allele lowers FHB severity). QFhb.yas-1AL is only cited; no position. |
| Semagn et al. 2023 Front Plant Sci 14:1190358 (FHB) | Z_1368 | 12 QTL checked against Table 2 (v2.0); 3 missing positions added, 3 wrong positions corrected, reference set to v2.0. |
| Yang et al. 2022 Front Plant Sci 13:990915 (STB) | Wai7A | 6 QTL placed from text regions or DArT-seq tag positions; 6 have only cM positions and stay unplaced. |
| Genievskaya et al. 2020 PeerJ 8:e9820 (rusts) | OA_Genievskaya2020genome | 29 of 44 MTAs placed from 90K SNP positions; 15 not placed (SNP maps to the homoeologous chromosome or has no position). |
| Running et al. 2022 Front Plant Sci 13:793925 (tan spot, Tsc1) | Running2022 | QTs.zhl-1A placed at the Tsc1 interval (the paper equates them). |
| Ren et al. 2019 TAG 132:241 (FHB, 7DL) | Z_1349 | QFhb.cau-1AS not placed: linked marker IWB63682 has no known physical position. |
| Shi et al. 2023 Agronomy 13:1712 (FHB GWAS) | Agron2023FHB / A_FHB_038 | 7 missing gFHB loci of Table 1 added (v1.1); cited Qfhb.nmbu.1A.1 placed at 3.1-13.3 Mb. |
