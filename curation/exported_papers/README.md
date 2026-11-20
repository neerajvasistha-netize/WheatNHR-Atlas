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
