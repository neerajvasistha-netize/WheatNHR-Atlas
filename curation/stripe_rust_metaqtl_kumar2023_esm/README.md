# Stripe rust QTL and MTAs from the Kumar et al. 2023 supplementary file

Source: 12864_2023_9336_MOESM1_ESM.xlsx, supplementary of the 2023 BMC Genomics stripe rust meta-QTL study (Kumar et al.).
Only the original QTL (Table S1) and GWAS MTAs (Table S5) are used, credited to their own papers; meta-QTLs, singletons
and hotspots are not added.

1. Each study ID in Table S1/S5 (e.g. Agenbag2012) was resolved on Crossref (cr_resolve.py: first author, year, title
   about stripe/yellow rust). Ambiguous IDs (a/b suffixes resolving to one paper, wrong-topic hits) were dropped.
   Agenbag 2012, Paillard 2012, Powell 2013 and Miedaner 2020 were confirmed by hand. Result: keep.json.
2. Studies already in the atlas (same DOI) were skipped.
3. Table S1 QTL: the interval markers were looked up for published physical positions (or e-PCR); a QTL is kept when at
   least one marker is placed and placed markers lie within 60 Mb (s1place.json). Table S1 gives no original QTL names,
   so records are named "QYr <study> <chr> (<marker interval>)".
4. Table S5 MTAs: published marker position if found, else the bp value of Table S5 (Aoun 2021 GBS tags and Miedaner
   2020 DArTseq markers use the Table S5 positions, taken as RefSeq v1.0). Arm letters in chromosome names were dropped.
Result: seed_ad78.json, 216 records (82 QTL, 134 MTAs) from 32 papers, grade C. Listed in added_records.tsv.

## Update with Jan et al. 2021 Table S1 (41598_2021_2049_MOESM1_ESM.xlsx, 353 stripe rust QTL with names)
- 77 QTL in seed_ad78 that were labelled "QYr <study> <chr> (<interval>)" now carry their published names, matched by
  paper DOI, chromosome and marker interval (generic names such as "QTL on 2AL" get the study added).
- With the real names, 21 of them turned out to be already in the atlas from the same paper (other seeds) and one sat on
  the wrong chromosome arm for its name; these 22 were removed (dropped_after_jan2021_names.json). seed_ad78 now holds 194 records.
- No further QTL were added from Jan et al. 2021: its QTL are either already in the atlas, from papers not identified
  reliably, or have no marker with a known physical position.
