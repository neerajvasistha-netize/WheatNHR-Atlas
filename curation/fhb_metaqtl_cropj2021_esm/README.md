# FHB QTL from the Crop Journal FHB meta-QTL study supplementary (S221451412030177X)

Files used: mmc2 (Table S1, 624 original FHB QTL with physical positions of interval and peak markers, RefSeq v1.0),
mmc5 (Table S4, reference list), mmc4 (Table S3, 32,074 consensus-map markers with RefSeq v1.0 positions).
Only the original QTL are used, credited to their own papers; meta-QTLs and the candidate-gene, marker-development,
expression and proteome tables (mmc6-mmc16) are not added.

- References [1]-[113] resolved on Crossref from the full citation (crbib.py; first author, year within 1, title word
  overlap >= 60%): 93 papers. 15 unresolved references were skipped.
- Skipped as already in the atlas: same DOI, or an FHB source in the atlas with the same first author and year
  (17 papers, skipped_possible_duplicates.txt).
- Result: seed_ad79.json, 417 QTL from 75 papers (grade C), named "QFhb <first author> <chr> (<marker interval>)".
- consensus_marker_positions_refseq_v1.tsv (from mmc4) is kept as a marker position source; it placed QPm.caas-7A
  (Xbarc174, 7A 116.00 Mb) from the Desktop/Papers batch (PY0018 in seed_ad75).
