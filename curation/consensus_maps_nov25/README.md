# Placing arm-only loci with consensus genetic maps (25 Nov 2026)

Aim: give a RefSeq v1.0 position to atlas records that had a chromosome but no position, using the markers they name.

## Sources
Supplementary files downloaded by the owner (Downloads, 25 Nov 2026), plus one public map fetched here:

| Map | Markers | Anchored on RefSeq v1.0 |
|---|---|---|
| Integrated wheat consensus map, Supplementary Tables 2-8 of a grain yield/quality meta-QTL study (only its marker map is used; no meta-QTL enters the atlas) | 143,940 | 30,869 |
| Wen et al. 2017, Front Plant Sci 8:1389, Table S6 (four populations, 90K array) | 52,607 | 31,954 |
| Maccaferri et al. 2014, BMC Genomics, durum consensus framework map (Table S1) | 2,524 | 1,034 |
| Marone et al. 2012, Theor Appl Genet 125:1619, durum consensus map (ESM_2) | 1,898 | 640 |
| KOMUGI (NBRP Japan) Composite Wheat Map, ChromosomeMap.xls (public download) | 4,428 | 582 |
| Somers et al. 2004, Theor Appl Genet 109:1105, primer list of 1,235 SSRs (ESM) | primers | used for in-silico PCR |

The GrainGenes site (ITMI and Wheat Composite 2004 maps) refused automated access (bot check), so the KOMUGI composite map was used for the
old RFLP probes instead. Not used: Maccaferri et al. 2016 root QTL supplement (erw039) and Marone et al. 2018 review.

## Method
1. `parse.py` reads the maps into one table (map, chromosome, marker, cM); `anchor.py` adds the RefSeq v1.0 position of every marker that already
   has one in the atlas marker tables (`look2.py`): 90K/35K/660K SNPs, DArT tag BLAST, Zhao et al. 2019 SSRs, Crop J 2021 consensus table, e-PCR.
   Result: `consensus_maps_anchored.tsv.gz`.
2. `interp.py`: per map and chromosome, anchors out of collinear order are removed (longest non-decreasing run of Mb along cM). A marker without a
   position takes the Mb range of at least two anchors within 1 cM of it, otherwise the two flanking anchors (each within 2 cM); ranges wider than
   10 Mb are not used; placements from different maps must agree within 20 Mb.
   Leave-one-out test (`loo.py`): for collinear markers the median error is 0.7 Mb (1 of 648 > 10 Mb). Across all anchors about 1 in 10 markers
   sits away from its genetic neighbours (paralogues, mis-mapped SNPs), so interpolated positions are marked for checking.
3. Somers et al. 2004 primers of 66 unplaced SSRs were run as in-silico PCR on their target chromosomes (relaxed 3' match, `epcr/`); unique products
   placed wmc11, wmc500, barc10, wmc471, cfd36 and gwm295.
4. `place_records.py`: for each record, markers with a known position first, then e-PCR, then interpolation; several markers within 30 Mb give an
   interval. Placements on the arm opposite to the one the source states (approximate centromere positions) were rejected (8 records). Records whose
   notes record an earlier rejection, and S19-081 (paralogous), were not changed.
5. `apply.py` + `conv96.py`: positions written to the records (grade C, curation "auto (in-silico from markers, consensus maps Nov 2026)", note
   lists each marker, its position and its source); RefSeq v2.1 and 13 pangenome assemblies converted as for earlier batches.

## Result
127 records placed (`placements.tsv`): seed_cgsw 37, seed_ad62 28, seed_ad66 10, seed_ad94 4, seed_ad95 2, seed_ad60 1, and 45 live-database
records (updated in the live database and in `web/seed_live.json`). Of the 1,228 records without a position, 471 name markers at all; the others
(no markers, AFLP/RAPD only, alien chromosomes) cannot be placed this way.

## Abeysekara et al. 2012 (Phytopathology 102:94-104), PDF supplied by the owner, 25 Nov 2026
QSnb.fcu-7AS (S26-101) and its Catalogue duplicate K0507: interval read from Fig. 3B (peak at XBF293121, flanked by Xgwm60.1/Xbarc154
and Xgwm130.3/Xwmc83): Xgwm60 54.96 Mb to Xwmc83 89.97 Mb, peak XBF293121 61.43 Mb (in-silico PCR of the Table 2 primers), RefSeq v1.0.
K0507 also given marker XBF293121 and LOD 4.64; K0506 (QSnb.fcu-1A, Snn4) given LOD 7.03 (Table 5). S26-101 updated in the live database.
