# 2026 papers added on 24 Nov 2026 (seed_ad84)

The owner noted that recent 2026 studies were missing. Release v1.0 already held 33 studies from 2026; this pass looked for more.

## How
1. `cr26.py` searched Crossref for journal articles issued from 1 Jan 2026 with 27 wheat disease queries (3,207 hits).
   Titles were filtered for wheat, a disease and a mapping term, and DOIs already cited in the atlas were removed.
   78 mapping papers were kept by hand (reviews, cultivar surveys, pyramiding, pathogen studies and corrections dropped);
   3 of them were already in the atlas under another DOI, leaving 75.
2. `fetch.py` and `fetch2.py` collected each paper's abstract (OpenAlex, Europe PMC) and the open full text where the
   publisher allowed it (about 20 papers).
3. Loci were read by hand from the abstracts and open texts into `records84.py`:
   positions as published (grade A or B, with the assembly), the RefSeq v1.0 position of a candidate gene named by the
   authors when no interval was printed (grade C), or the chromosome or arm only (grade D).
   Where a paper did not name the assembly, the assumed one is stated in the record note.
4. `build84.py` wrote the records, sources and paper-tracker rows; pangenome positions were added with the conversion
   used for earlier seeds (other-chromosome results dropped). Output: `web/seed_ad84.json` (75 records, ids Y6-001 onward,
   36 papers).

## Not added
Insect pests (wheat stem sawfly, Russian wheat aphid, Hessian fly), loci on alien chromosomes without a wheat position
(Thinopyrum, Elymus, rye, Dasypyrum, Haynaldia and Aegilops segments), alleles of genes already in the atlas
(Pm4, Pm26, Pm36, Yr6, Yr52, Yr83, Fhb7), and GWAS papers whose marker tables are behind a paywall. These are in
`tracker84.tsv` as "To search".

## Files
- `tracker84.tsv`: the 75 papers with status and notes.
- `added_records.tsv`: the 75 records added.
