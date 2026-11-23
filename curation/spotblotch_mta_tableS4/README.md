# Spot blotch / black point MTAs from Table S4 of the 2024 Sci Rep meta-QTL study

Source file: Supplementary Table S4 (41598_2024_63924_MOESM4_ESM.xlsx) of "Meta-QTL analysis and identification of
candidate genes for multiple-traits associated with spot blotch resistance in bread wheat", Scientific Reports 2024
(doi 10.1038/s41598-024-63924-w). The table lists significant markers (MTAs) from earlier GWAS papers.

Only the original GWAS MTAs are used, credited to their own papers; the meta-QTLs are not added.

- Sections used: SBR (spot blotch) and BPR (black point). SG (stay green) is not a disease trait and is skipped.
- Papers added (were "To search" in the tracker): Adhikari et al. 2012, Ahirwar et al. 2018, Gurung et al. 2014,
  Jamil et al. 2018, Kumar et al. 2022, Lozano-Ramirez et al. 2022, Liu et al. 2017 (black point). 93 MTAs in seed_ad76.json, grade C.
- Papers skipped because they are already in the atlas from the original paper: Ayana 2018, Bainsla 2020, Juliana 2022,
  Tomar 2021, Tang 2022, Lv 2020.
- Not added: Singh et al. 2020 (paper not identified), Singh et al. 2023 (marker IDs are panel-specific and the Table S4
  values look like map positions), Li et al. 2019 (no marker names).
- Positions: the marker's published physical position from our marker table where available; otherwise the bp value in
  Table S4, used only for papers whose looked-up markers agree with it (Gurung, Jamil, Lozano-Ramirez, Liu, and Kumar
  whose KASP markers could not be looked up). For Adhikari and Ahirwar, Table S4 gives genetic-map values, so only
  looked-up markers were kept. One Lozano-Ramirez row whose position equals the marker ID was dropped.

Files: build76.py (python3 build76.py tableS4_rows_with_lookup.json crossref.jsonl out.json), added_mtas.tsv.
