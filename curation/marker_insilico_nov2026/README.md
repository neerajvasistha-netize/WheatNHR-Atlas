# In-silico placement of loci from their molecular markers (19 Nov 2026)

Scope: every atlas record without a physical position that names markers (630 records).

1. **Markers with a known physical position** (atlas marker table `data/markerpos_all.json.gz`: 90K/35K/55K/660K SNP maps,
   DArT tag BLAST, Zhao et al. 2019 SSR positions) or with the position written in the marker name (chr_position SNPs).
   Earlier passes missed these because of name variants (Xiwa/Xiwb, numeric DArTseq clone IDs, wmc120.2-type suffixes).
2. **Marker sequences found and mapped on RefSeq v1.0** (`scanq.py`, every chromosome streamed from Ensembl Plants):
   - BARC SSR primer pairs (USDA-ARS BARC list, scabusa.org, `barc_primers.json`): in-silico PCR, product 60-800 bp,
     up to 2 mismatches per primer outside the 3' 12 bases. 7 of 18 needed pairs amplify.
   - RFLP probe sequences from GenBank (BCD/CDO/ABG/MWG/WG/PSR clones, `rflp_seqs.json`): 16-mer seed clusters in 5 kb
     windows; the copy on the record's chromosome is used (the homoeologous copies are listed in `scan_results.json`).
3. Rules: several markers within 30 Mb give an interval, otherwise the first-listed marker gives the position; grade B;
   placements that contradict the chromosome arm stated by the source were rejected (Xwmc120.2 1AS, QLr.fcu-3AL, QKb.cim-3DL);
   names that are marker IDs rather than positions (Navathe et al. chr:ID; Kazakh DArTseq IDs) and a locus earlier rejected as
   paralogous (S19-081) were not used.

Result: 75 records placed (`records.tsv`); 51 in `web/seed_ad68.json`, 24 live-database records updated in the live database
and in `web/seed_live.json`. Not placed: wPt DArT markers with no public sequence (about 100), SSRs without an accessible primer
sequence (GWM 8xx/9xx, many WMC, CFA, CFD), AFLP and KSU/FCP RFLP probes without sequence, and SNPs whose physical position
lies on a homoeologous chromosome.
