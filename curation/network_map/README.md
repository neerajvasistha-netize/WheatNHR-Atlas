# Gene network and functional response map (Genomic tools)

Input: any gene list (all ID types of the enrichment tool), optional log2 fold change column.

Functional response map: genes are assigned to functional categories by keyword rules on the IWGSC gene
descriptions (the rules are the NMCAT list in web/index.html), grouped by the cell compartment where each
category mainly acts (apoplast and cell wall, plasma membrane, cytoplasm, nucleus, chloroplast). One cell per
gene, coloured by log2 fold change. The layout is the atlas's own design.

Gene network: genes are linked by
- shared function: 2 or more specific GO terms (each annotated to at most 300 genes) or a KEGG pathway;
- co-expression: Pearson r >= 0.8 of log2 fold changes across the atlas DEG contrasts (expr_degs.json),
  at least 3 shared contrasts, each gene keeping its 5 strongest partners;
- sequence similarity: CDS identity >= 90% (web/sim/similar_<chr>.json).
Layout: Fruchterman-Reingold per connected component, components shelf-packed by size. Edge list exports as CSV
for Cytoscape.
