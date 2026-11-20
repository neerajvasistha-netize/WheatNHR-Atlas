# Papers from the Zotero "Exported" folder (Nov 2026)

The owner's Zotero export (`~/Downloads/Exported /`, 660 PDFs of biotrophic disease papers) was read automatically
on 20 Nov 2026. PDFs are not kept here.

## What was done
1. Every PDF was converted to text with `pdftotext` (layout mode for tables, reading order for sentences).
2. Papers already cited in the atlas were set aside by title (94); 566 papers were new.
3. `rdf.py` read authors, journal, volume, pages and DOI from the Zotero RDF files.
4. `s1.py` and `s2.py` found the loci each paper reports: formally named QTL of the paper's own lab code
   (names found in the title or abstract) for rust, powdery mildew and bunt traits, and designated genes named in the title.
   For each locus the chromosome comes from the name or the text. The position was taken in this order:
   a Mb interval written right after the locus name; a base-pair position in the locus's table row after its chromosome;
   or the published physical positions of at least two linked markers named with the locus
   (`data/markerpos_all.json.gz`, at most 25 Mb apart). Positions on the wrong chromosome arm were rejected.
5. `build73.py` kept loci not already in the atlas under the same name and chromosome and wrote `web/seed_ad73.json`
   (199 records from 54 papers, ids PX0001 onward, grade C, curation "auto (paper text)").
   Each record carries the sentence or table row the position came from (`evidence`).
   Pangenome positions were added with `scripts/pangenome/convert.py` (as for ad66 and ad72).
6. All 566 papers were added to the paper tracker (ids `PX_<hash>`): "Added to atlas" when loci were placed,
   "To search" otherwise, with notes on loci named but not placed.

## Files
- `tracker73.tsv`: one row per new paper (status, records added, loci placed and unplaced).
- `unplaced_loci.tsv`: loci the reading found but could not place (no position in the text, no chromosome, or wrong arm).

## Limits
Records are machine-read and should be checked against each paper ("verify against paper").
Loci shown only in figures or supplementary files, alien-chromosome lines without wheat positions,
and GWAS marker-trait associations without locus names were not added.

## Second pass: markers read from figures (20 Nov 2026)
For the 404 loci (in 100 papers) that were named but had no position in the text and are not in the atlas,
the figure pages (pages with images or a figure caption, 385 pages) were rendered at 300 dpi and read with
Tesseract OCR (`pages.py`, `ocr.py`). `fig2.py` kept, for each locus, the marker names in the figure that have a published
physical position on the locus's chromosome (and arm, where the name or text gives the arm), choosing the markers drawn
nearest the locus label (or all markers of a figure whose caption names the locus), at least two markers within 15 Mb.
`addfig.py` added the 14 loci placed this way to `web/seed_ad73.json` (PX0200 onward, curation "auto (paper figure)");
see `figure_placements.tsv`. The other loci stayed unplaced: most figures are Manhattan plots or maps without
marker names on the locus chromosome, or their markers have no known physical position.

## Third pass: supplementary tables (20 Nov 2026)
`dl.py` downloaded the supplementary files of the 96 papers in `supplementary_needed.bib` from Europe PMC
(open-access papers) or the publisher page (Springer ESM, Frontiers); files were found for 66 papers
(`supplementary_download_log.json`; files not kept here). `rows.py` read every Excel, Word, PDF and text file into table rows and
`supp2.py` looked for the row of each unplaced locus: a physical position written next to the locus chromosome
(e.g. "3D: 606366088" or chromosome and bp columns), or, failing that, markers in the row with a published position on that
chromosome and arm. `addsupp.py` added the 106 loci placed this way to `web/seed_ad73.json`
(curation "auto (supplementary table)"; `supplementary_placements.tsv`).
`supplementary_manual.bib` / `.tsv`: the 30 papers whose supplements could not be fetched automatically (MDPI and some
other publishers block downloads; several older papers may have no supplement).

## Fourth pass: SNP names that carry positions (20 Nov 2026)
The owner supplied the supplement of Sharma et al. 2026 (Green Hammer x Lonerider). It holds genotype data and the linkage
map but no QTL positions; the QTL table of the paper itself gives flanking GBS markers named by their physical position
(e.g. S2B_45213215). `gbs.py` read such marker names from the QTL table rows (rejoining QTL names split over two lines)
for all papers with unplaced loci, and `addgbs.py` added 17 loci (15 from Sharma et al. 2026, 1 from Sharma et al. 2026
Baker's Ann, 1 from the QYraq.cau paper) to `web/seed_ad73.json` (`snpname_placements.tsv`).

## Fifth pass: markers without a known position (20 Nov 2026)
`unk.py` listed 133 markers named next to still-unplaced loci that have no published physical position in
`data/markerpos_all.json.gz` (`unknown_markers.json`): about 70 DArT (wPt), 48 SSR, 12 SNP-array and a few KASP names.
Web searches gave no position tables; GrainGenes blocks automated access. Sequences could be obtained for 16 SSRs
(primer pairs from the atlas primer table, STS records from NCBI); `epcr.py` searched them on the locus chromosomes of
RefSeq v2.1 (Ensembl Plants 63) by exact primer and 32-mer matching. Six gave a single position (`epcr_res.json`):
barc59 5B 670.70 Mb, barc261 5B 13.35 Mb, barc175 6D 411.88 Mb, barc202 6D 291.70 Mb, barc96 6D 418.28 Mb, barc121 7A 611.83 Mb.
None of them, with the other markers of its locus, gave a consistent interval (`place2.py`), so no locus was added.
