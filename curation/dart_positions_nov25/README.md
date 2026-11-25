# DArT and DArTseq marker positions on RefSeq v1.0 (25 Nov 2026)

Input (supplied by the owner, not committed): `Wheat_Concensus_Map_4_0_with_sequences_3.csv` (DArT consensus map 4.0, 107,652 markers with
chromosome, cM and tag sequence) and `DArT_Wheat.fasta.txt` (7,818 wPt clone sequences).

- `aln.py`: every sequence (110,657: map 4.0 tags plus wPt clones not on the map) aligned with minimap2 (mappy 2.31, `sr` preset, best 3 hits)
  to each RefSeq v1.0 chromosome from Ensembl Plants, in 300 Mb chunks to stay within memory.
- `post.py`: a hit counts if it covers >= 80% of a tag (60% of a long clone) at >= 90% identity, or >= 60% (36%) at >= 95% identity. Best hit
  gives the position; equal hits resolved by the map 4.0 chromosome or treated as one region when within 2 Mb; otherwise ambiguous.
- Result `dart_marker_positions_refseq_v1.tsv.gz`: 71,962 markers placed (65,298 of 102,686 DArTseq tags; 6,664 of 7,971 wPt/tPt/rPt);
  physical chromosome agrees with the map 4.0 chromosome for 83.7%. Agreement with the earlier atlas DArT tag BLAST: 39,907 same, 354 differ.
- Added to `data/markerpos_all.json.gz` as "DArT consensus map 4.0 sequence alignment (this atlas)" (25,023 markers new to the table; existing
  positions keep priority). Re-running `curation/consensus_maps_nov25/place_all.py` placed 16 more atlas records (seed_ad96.json; 10 of them
  also updated in the live database and seed_live.json).
