# Correction log: 19 spot blotch / leaf spot papers re-read (11 Nov 2026)

Every locus below was re-read from the paper's own results table and placed on RefSeq v1.0 (or kept on v2.0 where the paper prints v2.0), then lifted to RefSeq v2.1. Records carry IDs `S19-001` to `S19-160` (S19-100 to S19-115, the Vasistha 2024 meta-QTL, were dropped); the 67 records they replace are retired (marked `deleted` in the live database, removed from `data/merged_v13f.json`). Scripts: `scripts/curation/sb19_fix.py` (the records) and `sb19_apply.py` (the merge).

SSR markers were placed by exact-match in-silico PCR of their primers on the IWGSC RefSeq v1.0 chromosomes from Ensembl Plants (`scripts/curation/sb19/ispcr.py`, results in `ispcr_res.json`). A marker is used only if it amplifies on the chromosome the paper states.

## Records per paper

| Paper | Old records | New records | Grade A / B / C / D |
|---|---|---|---|
| Gahtyari et al. 2021 | 5 | 19 | 19 / 0 / 0 / 0 |
| He et al. 2020 | 3 | 13 | 13 / 0 / 0 / 0 |
| Kaur et al. 2021 | 11 | 5 | 5 / 0 / 0 / 0 |
| Kumar et al. 2009 | 4 | 4 | 0 / 2 / 2 / 0 |
| Kumar et al. 2010 | 7 | 9 | 0 / 2 / 5 / 2 |
| Kumar et al. 2015 | 2 | 1 | 0 / 0 / 1 / 0 |
| Kumar et al. 2016 | 0 | 1 | 0 / 0 / 1 / 0 |
| Lillemo et al. 2013 | 2 | 4 | 0 / 2 / 2 / 0 |
| Lu et al. 2016 | 1 | 1 | 0 / 1 / 0 / 0 |
| Pankaj et al. 2022 | 2 | 16 | 0 / 4 / 7 / 5 |
| Roy et al. 2021 | 9 | 11 | 0 / 5 / 4 / 2 |
| Singh et al. 2016 | 2 | 2 | 0 / 1 / 1 / 0 |
| Singh et al. 2018 | 2 | 13 | 0 / 1 / 10 / 2 |
| Vasistha et al. 2024 | 7 | 0 (meta-QTL not kept) | – |
| Zhang et al. 2020 | 1 | 1 | 1 / 0 / 0 / 0 |
| Zhu et al. 2014 | 4 | 5 | 0 / 1 / 4 / 0 |
| Bemister et al. 2019 | 5 | 7 | 0 / 2 / 5 / 0 |
| Iqbal et al. 2023 | 0 | 32 | 32 / 0 / 0 / 0 |
| Gao et al. 2023 | 5 | 5 (kept, checked against Table 2) | unchanged |

## Errors found and corrected

- **Gahtyari et al. 2021**: the atlas held 5 of the 19 QTL in Table 4. Two were wrong: the KC 4D QTL (Rht-D1–BS00036421) had been put on **2B at 1.07 Mb** because the atlas marker dictionary had a false 2B entry for BS00036421 (true position 4D 32.35 Mb; dictionary fixed), and two separate 4B QTL had been merged into one 535–572 Mb record. Sb4 had been credited to this paper; it belongs to Zhang et al. 2020. All 19 QTL are now in, with the physical ranges printed in Table 4.
- **He et al. 2020**: 3 of 12 QTL were present; one was named "Chirya x CIMMYT" although the population is CASCABEL × CIANO T79, and one record labelled Sb2 carried the 315–447 Mb QSb.bhu-5B interval. All 12 Table 3 QTL added; the 5BL cluster (404.2–418.8 Mb) kept.
- **Sb2 (Kumar et al. 2015)**: the atlas showed Sb2 at 315.5–447.2 Mb, which is the coarse Kumar 2009 QSb.bhu-5B interval. Sb2 lies at Xgwm639, **504.30 Mb** (in-silico PCR), consistent with its 5BL1-0.55-0.75 bin and with Tsn1 marker Xfcp1 (546.15 Mb) mapping 1.1 cM away (Kumar et al. 2016).
- **Kaur et al. 2021**: 5 QTL were stored 11 times (duplicates from three extraction passes), one without a position, and a Jamil et al. 2018 SNP (S5B_679369233) that Kaur only cites was attributed to Kaur. Now 5 records. Table 2 and Table 3 give different SNPs for QSb.pau-3B (S3B_104700872 vs S3B_104700839); noted.
- **Lillemo et al. 2013**: the source key carried the citation of a different paper (Saxesena et al. 2022). Replaced by key `Lillemo2013TAG` with the correct citation. Sb1 placed on the Lr34 gene (47.41 Mb); the two minor 7DS QTL and the 1BL QTL at Lr46 added.
- **Pankaj et al. 2022**: the only record was Kaur 2021's S5B_703858864, which Pankaj only cites; the Research Square preprint duplicated it. Replaced by all 16 distinct QTL of Table 3. Many flanking SSRs do not amplify on the chromosome the paper assigns (e.g. Cfd71 → 4D, wmc173 → 4A for "5B"; gwm131 → 7B for "1B"); those QTL are graded C or D with the reason in the note.
- **Kumar et al. 2010**: QSb.bhu-2D, -3B and -7B were stored without their population (they are from Chirya 3 × Sonalika) and two had no position. The Chirya 3 QTL on 2B and 7D were missing. Xswm008 amplifies on 7D at 51.45 Mb (the old dictionary listed 4A/2B/3D), giving the Chirya 3 7D interval 13.46–51.45 Mb. Xbarc159 primary locus is 2D, so the 2A QTL is now anchored on Xgwm425 only.
- **Roy et al. 2021**: the two largest-effect QTL were missing: BC 4D (PVE up to 12.2%) and CC 5A 331–332 cM (PVE up to 12.6%, now 564.5–565.4 Mb). The CC 5B QTL flanks place 193 Mb apart for a 3 cM interval and are left unplaced.
- **Singh et al. 2018**: only two unplaced 4B/4D records existed; all 13 QTL added.
- **Singh et al. 2016**: QSb.iiwbr-7B is called 7BL in the paper but both flanks (Xwmc758 164.6 Mb, Xwmc335 233.2 Mb) are on 7BS. Xwmc653 amplifies only on 3B, not 7D.
- **Vasistha et al. 2024**: meta-QTL are not kept in the atlas (decision 11 Nov 2026). The 7 preprint records (Sharma & Vasistha 2022) were removed and the 16 MQTL were not added. Noted while reading: MQTL3-1B candidate genes (16.9–17.9 Mb) and MQTL8-5A genes (about 658 Mb) do not match the positions of their Table 1 flanking markers (604–623 Mb and 575.6 Mb).
- **Bemister et al. 2019**: the atlas had five of its markers as separate *tan spot* records taken second-hand from Liu et al. 2020. The paper scored the field leaf spot complex, not tan spot. Replaced by its 3 leaf spot, 2 leaf rust and 2 stripe rust QTL. New disease category "Leaf spot complex".
- **Iqbal et al. 2023**: not in the atlas; all 32 QTL added (RefSeq v2.0 positions as printed): leaf spot, leaf rust, stripe rust and common bunt.
- **Zhu et al. 2014**: QSb.cim-3B flank corrected to 990937–1123330 (the old record used 1002801); PVE added; unanchored QSb.cim-14 added, chromosome inferred (6B) and flagged.
- **Lu et al. 2016 (Sb3)** and **Zhang et al. 2020 (Sb4)** were already correct; confirmed by in-silico PCR (Sb3 6.09–6.45 Mb) and Ensembl gene coordinates (Sb4 from 581.07 Mb). Gao et al. 2023 matched Table 2.

## Still unplaced or needing data

- No primers available for Xgwm533, Xgwm1037, Xgwm263, Xgwm255, Xgwm815, Xgwm732, Xgwm1043, Xgwm1168, Xgwm1220, wmc728, cfd59 (GrainGenes is behind a browser check). QSb.bhu-3B and -7B (Chirya 3) therefore stay unplaced.
- Xbarc175 (6D 411.88 Mb) and Xbarc121 (7D 531.59 Mb) positions come from an earlier session; exact-match PCR now gives no product, so treat them as provisional.
- DArTseq tags not in the atlas dictionary leave several Singh 2018 and Roy 2021 flanks unplaced; their supplementary tag sequences would resolve them.

## Wider problems noticed

- The live database holds 180 documents (R001–R180) written on 31 Oct from release v1.2. The page lets these override the newer seed records with the same ID, so some hand corrections made since (e.g. Sb1 had been shown as "arm only") were hidden. Those 13 that belonged to these papers are now retired; the other ~167 still override the seed and should be reviewed.
- One existing atlas record lifts v1.0 3A 549.6 Mb to v2.1 112.5 Mb; this anchor was excluded from the new liftover and should be checked (records near 3A 549.6 Mb).
- The GitHub `web/index.html` is older than the published page (the published one has the curation-notice panel). The seed edits here were made on the GitHub copy; the live page was corrected through its database.
