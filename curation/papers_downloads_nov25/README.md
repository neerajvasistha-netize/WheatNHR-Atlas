# Loci from the owner's downloaded papers (Downloads/Paper_2026), 25 Nov 2026

Batch 1: the 160 older papers (1986-2015) in `Downloads/Paper_2026/Papers` that matched rows of the
"not downloaded loci papers" list (`tracker94.tsv`; matching from the list's "Removed - downloaded" sheet).

- 27 papers were skipped because their DOI is already a source in the atlas, 4 were duplicate files.
- Text was taken with pdftotext; five scanned PDFs were read by OCR (tesseract).
- `records94.py`: one entry per locus (gene, chromosome or arm, markers named by the paper, evidence quote).
- `build94.py`: places loci in silico on RefSeq v1.0 from the paper's markers (`look2.py`: atlas marker table
  `data/markerpos_all.json.gz`, then the consensus table of `curation/fhb_metaqtl_cropj2021_esm`); two placed markers
  within 30 Mb give an interval, otherwise the first placed marker; grade C. Loci whose markers have no known
  position (RFLP, RAPD, SCAR, AFLP, many SSRs) keep the chromosome arm (grade D). Four placements were rejected because
  the marker maps to the other arm than the paper states (listed in `build94.py`).
- Citations from OpenAlex (`oa.py`). Pangenome conversion with `conv94.py` (scratchpad copy of conv66.py).
- Result: `web/seed_ad94.json`, 165 records (54 grade C, 111 grade D) from 116 papers, ids PD6-001 to PD6-165.

Batch 2: the 2013-2026 papers in `Downloads/Paper_2026/Paper_2026_1` and the alien-introgression papers in
`Downloads/Paper_2026/Alien` (523 text files; `tracker95.tsv` lists every file and its status).

- 69 papers skipped because their DOI is already an atlas source, 41 more because their title matches an atlas
  citation, 6 had no extractable text; exact duplicate copies were not staged.
- 407 new papers were read (abstract and the sentences giving positions); loci in `records95.py`.
- `build95.py`: positions written by the paper are used as published (grade B interval, A peak; assembly as stated,
  otherwise RefSeq v1.0 for studies up to 2022 as noted in each record); otherwise in silico from named markers (grade C);
  otherwise the arm (grade D). Alien-introgression loci are kept only with the wheat chromosome that carries the segment
  (`ref` says "alien (...)"). Loci with no chromosome in the text read (19) and one student essay were left out.
  Citations from OpenAlex, or Crossref bibliographic search for PDFs without a DOI (six wrong matches corrected by hand).
- Result: `web/seed_ad95.json`, 307 records (15 A, 54 B, 25 C, 213 D) from 255 papers, ids PE6-001 to PE6-307.
  New disease added to the colour list: Flag smut.
- Multi-locus GWAS and QTL tables were not transcribed in full: only loci named in the abstract or with a position in
  the text are included. Supplementary files staged with this batch (pdis-02-25-0370, pdis-11-25-2203, pdis-12-25-2526,
  phyto-04-26-0105) are not yet used.
