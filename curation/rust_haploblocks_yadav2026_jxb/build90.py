"""Yadav S et al. (2026) J Exp Bot 10.1093/jxb/erag073: consistent haploblocks for stripe rust (Supplementary Table S4)
and leaf rust (Table S5) in the OzWheat panel; positions on IWGSC RefSeq v1.0 (stated for Table S6)."""
import openpyxl,json,sys
wb=openpyxl.load_workbook(sys.argv[1],read_only=True,data_only=True)
recs=[]
def add(dis,blk,ch,s,e,n,var,known,lit,rank=None):
    tag='Yr' if dis=='Stripe rust' else 'Lr'
    kn='' if not known or str(known).strip() in ('None reported','None') else f' Overlap reported by the authors: {str(known).strip()}'+(f' ({str(lit).strip()})' if lit and 'No overlap' not in str(lit) else '')+'.'
    if 'novel' in str(known).lower(): kn=' The authors call it a putative novel locus.'
    recs.append({'name':f'Haploblock {blk} ({tag})','cls':'QTL','disease':dis,'chr':ch,'start':round(s,6),'end':round(e,6),'peak':None,'ref':'RefSeq v1.0','marker':f'{n} SNPs','pve':'','stat':f'variance of haplotype effects {var:.6f}',
      'grade':'B','evidence':f'Supplementary Table {"S4" if tag=="Yr" else "S5"}: {blk}, {ch}, {s}-{e} Mb, {n} markers'+(f', rank {rank}' if rank else ''),
      'note':f'Haplotype block with consistent effects on {dis.lower()} severity across field trials (OzWheat panel, local GEBV haplotype mapping).{kn} Positions on RefSeq v1.0.'})
for r in list(wb['Supplementary Table S4'].iter_rows(values_only=True))[2:]:
    if r[1]: add('Stripe rust',r[1],r[2],r[3],r[4],r[6],r[7],r[8],r[9])
for r in list(wb['Supplementary Table S5'].iter_rows(values_only=True))[2:]:
    if r[2]: add('Leaf rust',r[2],r[3],r[4],r[5],r[8],r[9],r[10],r[11],r[0])
json.dump(recs,open('hb90.json','w'),indent=0);print(len(recs))
