# GrainGenes sequences supplied by the owner (WheatTraitMap_1.xlsx, 26 Nov 2026)

The owner filled primer pairs (88 markers) and clone/amplicon sequences (18) into the unpositioned-marker workbook from GrainGenes.
- `run.py`: in-silico PCR on all 21 RefSeq v1.0 chromosomes (exact 3' 15 bases of each primer, <= 4 mismatches over both primers,
  product 50-1500 bp) and minimap2 alignment of the sequences (>= 50% coverage, >= 90% identity). `post.py` keeps products on the chromosome the
  atlas record needs.
- `relax.py`: for primer pairs without a product on the target chromosome, seeds on the 3' 10 bases, <= 3 mismatches per primer.
- Result `final_pos.json`: 23 markers positioned (17 by PCR, 6 by sequence alignment). Most other primers give no product on the target
  chromosome of Chinese Spring (many amplify a homoeologue instead), so those loci stay unplaced.
- 17 atlas records placed (seed_ad96.json; S122-002, S44-027, S44-079, S138-004 also in the live database and seed_live.json). Rejected by the
  arm check: Xwmc120.2 (V3574), QPm.caas-2BS (K0379), Fhb7 S132-001 (text says 7DS); not used: Pm12 and Qkb.ksu-6BS.1 (Xwmc105 maps to 6BL),
  3AS QTL with Xwmc695b (3AL), Pm17 (rye 1RS).
- The filled workbook: `curation/WheatTraitMap_1_with_positions.xlsx`.
