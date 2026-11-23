# Rust QTL and MTAs from Pal et al. 2022 (TAG) and the Plant Genome leaf rust meta-QTL study (tpg2.20185)

Files: 122_2022_4119_MOESM1_ESM.xlsx (Pal et al. 2022, rust multi-disease meta-QTL: Table S2 = 1,364 QTL with marker
intervals and references, Table S7 = 1,895 GWAS MTAs with physical positions); tpg220185 Tables S2 (leaf rust QTL,
genetic positions only), S4 (62 GWAS MTAs with physical positions), S5 (leaf rust genes).
Only original QTL and MTAs are used, credited to their own papers; meta-QTLs are not added.

- References resolved on Crossref (cr_resolve2.py: first author, year, rust title); earlier hand-checked matches reused.
  Kept only exact-year matches whose title fits the record type; papers already in the atlas (DOI) skipped,
  including those added the same day from Kumar et al. 2023. 27 papers kept (keep.json).
- Table S2 QTL placed in silico from their interval markers (existing marker table, e-PCR, Crop J Table S3);
  dropped when no marker had a position (90) or placed markers were >60 Mb apart (7).
- MTAs: the marker's published position where found, else the bp value in the table (taken as RefSeq v1.0).
- tpg Table S2 has no markers or physical positions, so it was not used; all Table S5 genes are already in the atlas.
Result: seed_ad81.json, 473 records (41 QTL, 432 MTAs) from 20 papers, grade C.
