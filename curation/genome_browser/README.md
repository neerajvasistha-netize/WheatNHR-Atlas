# Genome browser and BLAST (Genomic tools)

Data in `web/gb/`, built by `build_gb.py`:
- `genes_<chr>.json`: gene models in `genes_all.json` order: v1.1 gene number, RefSeq v2.1 start (bp), length, strand,
  exons and CDS of the canonical .1 transcript (relative to the gene start). Gene structure and strand from the
  Ensembl Plants release 63 GFF3 (IWGSC RefSeq v1.1 annotation); the v2.1 start is the v2.1 gene midpoint
  (liftover cache) minus half the gene length.
- `cds_<chr>.txt`: coding sequences of the 106,7xx high-confidence genes, 2 bits per base (A0 C1 G2 T3) in base64,
  with a JSON header of gene numbers, offsets, lengths and any N runs.

Browser tracks: atlas loci (RefSeq v2.1), genes (candidate genes highlighted: inside a locus, or within 250 kb of a
single-position locus, with a defence-related gene class or annotation), markers (markers_lite), and the bases when
the view is 300 bp or shorter. Region FASTA (up to 5 Mb) comes from the atlas genome copy when hosted with it,
otherwise from NCBI (RefSeq v1.0, converted through anchor genes).

BLAST: blastn-style search run in the browser against the CDS store: exact seeds (word 11, 12 or 15), banded affine
gapped extension (reward 2, penalty -3, gap open 5, extend 2), Karlin-Altschul E-values (lambda 0.625, K 0.41) over
the whole database. Test: a 700-bp fragment of TraesCS6A02G051300 with two substitutions returns that gene first
(99.7%, CDS 101-800) followed by its homoeologs and paralogs.
