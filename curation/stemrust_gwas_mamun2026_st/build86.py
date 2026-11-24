"""Al Mamun et al. (2026) Plant Disease stem rust GWAS (361 accessions, exome-capture SNPs): 34 MTAs from
Supplementary Table S10 (file supplied by the owner on 24 Nov 2026). Positions are IWGSC RefSeq v1.0
(Table S11 lists RefSeq v1.1 gene models, TraesCS..02G, at the same coordinates)."""
import openpyxl,json,sys
wb=openpyxl.load_workbook(sys.argv[1],read_only=True,data_only=True)
rows=list(wb['Table S10'].iter_rows(values_only=True))[2:]
genes={}
for r in list(wb['Table S11'].iter_rows(values_only=True))[2:]:
    if r[1]: genes.setdefault(r[1],[]).append(r[5])
recs=[]
for i,(race,mk,ch,bp,mb,known,amk,pmb,ref) in enumerate(rows,1):
    known=(known or '').strip(); rel='' if known in ('Likely novel','not provided by the literature') else f' Near previously reported {known} ({ref}).'
    recs.append({'id':f'MM6-{i:03d}','name':f'MTA Sr {race} {ch} {mb:.2f}','cls':'MTA','disease':'Stem rust','chr':ch,'start':None,'end':None,'peak':round(float(mb),6),
      'ref':'RefSeq v1.0','marker':mk,'pve':'','stat':'significant (GWAS)','source':'Y26_pdis1225','grade':'A','curation':'auto (supplementary table)',
      'evidence':f'Supplementary Table S10: race {race}, {mk}, chromosome {ch}, {bp:,} bp; status: {known}',
      'note':f'Seedling resistance to Pgt race {race} in 361 diverse accessions (302,524 exome-capture SNPs).{rel} Candidate genes in Table S11: {len(genes.get(mk,[]))}. Position as published (RefSeq v1.0, matching the v1.1 gene coordinates in Table S11).'})
json.dump(recs,open('mta86.json','w'),indent=0);print(len(recs))
