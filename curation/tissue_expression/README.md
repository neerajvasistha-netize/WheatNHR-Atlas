# Organ expression data for the enrichment plant figure (web/tissue_expr.json)

Source: Expression Atlas experiment E-MTAB-4484 (Pingault et al. 2015, Genome Biology 16:29,
doi 10.1186/s13059-015-0601-9): RNA-seq of Chinese Spring root, leaf, stem, spike and grain, three developmental stages each.
File used: https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-4484/E-MTAB-4484-tpms.tsv
(per assay group: minimum, lower quartile, median, upper quartile, maximum TPM of the replicates).

`build_tissue.py`: for each IWGSC v1.1 gene that is in web/genes_all.json, the median TPM of each stage is averaged over the
three stages of an organ (Root g11, g4, g9; Leaf g10, g1, g18; Stem g5, g7, g3; Spike g6, g8, g2; Grain g17, g19, g20) and stored
as log2(TPM+1) with one decimal; genes below 0.5 TPM in every organ are left out (read as 0).

Used by the GO and KEGG enrichment tool, Figure 5: each enriched term's list genes are averaged per organ and drawn on a
wheat plant (original drawing) and in a term-by-organ heatmap.
