# Powdery mildew QTL and MTAs from the Data Sheet 1 of a wheat powdery mildew meta-QTL study

File: Data_Sheet_1.xlsx (Table S8 = 222 original PM QTL with marker intervals and references; Table S6 = 283 GWAS MTAs
with physical positions). Only original QTL and MTAs are used, credited to their own papers; meta-QTLs, singletons,
hotspots and candidate genes are not added.
- References resolved on Crossref (cr_resolve3.py) plus DOIs checked by hand for papers the automatic search missed
  or mismatched (Hinterberger 2022, Kang 2020, Lu 2012, Keller 1999, Mingeot 2002, Lan 2010, Liu 2001, Desiderio 2021,
  Mohler 2011/2019, Mourad 2023, Nordestgaard 2021, Jakobson 2005, Hao 2014). Papers already in the atlas (DOI) skipped.
- Table S8 QTL placed in silico from their interval markers (existing marker table, e-PCR, Crop J Table S3); dropped when
  no marker had a position (76) or placed markers were >60 Mb apart (10).
- MTAs: the marker's published position where found, else the bp value of Table S6 (taken as RefSeq v1.0).
Result: seed_ad82.json, 207 records (48 QTL, 159 MTAs) from 22 papers, grade C.
