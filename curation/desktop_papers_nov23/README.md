# Papers in ~/Desktop/Papers (read 23 Nov 2026)

184 PDFs (PDFs not kept here). 80 had already been read in an earlier batch (same file), 11 were duplicates, 7 had no
text (scanned), 16 were already cited in the atlas; 70 were new. Metadata from Crossref (DOI in the file name or the
PDF text, otherwise a bibliographic search on the first page; 27 papers have no metadata and use the title from the text).
Reading as in `curation/exported_papers/` (`s2.py`). `build75.py` wrote `web/seed_ad75.json`: 5 new loci (ids PY0001
onward, grade C, "auto (paper text)") and tracker entries for the 70 papers (`tracker75.tsv`).
Removed by hand: Lr34 placed on 2B (Lr34 is on 7DS and already in the atlas).
Most of these papers are 1986-2008 gene mapping studies with RFLP/SSR markers and cM maps; 54 named loci had no
position in the text and no linked marker with a published physical position.

## Second pass: in-silico placement from table markers (23 Nov 2026)
For loci named in the papers but not placed in the first pass, the marker interval or closest marker was read
from each paper's QTL table and placed at the published physical position of the marker (RefSeq v1.0), then
converted to RefSeq v2.1 and the 13 pangenome cultivars. 10 more loci were added (PY0007-PY0016, grade C), listed in
`insilico_placements.tsv`. `place75_candidates.py` lists candidate markers near each locus name; the final marker
for each locus was chosen by reading the table rows. Loci already in the atlas under the same name were skipped.
Loci whose markers have no known physical position (for example Xgwm768, Xgwm191, gwm707, cfd71b, gwm790b, gpw1106,
Xgwm304, Xbarc353b, Xbarc127, Xbarc174) remain unplaced.
