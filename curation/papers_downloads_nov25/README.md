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
