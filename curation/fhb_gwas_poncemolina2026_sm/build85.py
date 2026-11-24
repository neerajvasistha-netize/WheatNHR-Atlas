"""Ponce-Molina et al. (2026) Plant Disease FHB GWAS of PGRC accessions: QTN from Supplementary Tables S7, S9, S10
(supplement supplied by the owner on 24 Nov 2026). Positions are IWGSC RefSeq v1.0 (90K positions match published v1.0
90K coordinates). Only QTN for FHB traits are kept (incidence, severity, FHB index, FDK, DON); heading date and plant height dropped."""
import re,json,sys
T=open(sys.argv[1]).read().split('\n')
a=next(i for i,l in enumerate(T) if 'Supplementary Table S10' in l); b=next(i for i,l in enumerate(T) if 'Supplementary Table S11' in l)
FHB={'INC':'incidence','DS':'severity','IND':'FHB index','FDK':'Fusarium-damaged kernels','DON':'DON content'}
SITE={'19S':'Saskatchewan 2019','20S':'Saskatchewan 2020','21G':'Guelph 2021','21M':'2021 (M)','23G':'Guelph 2023','24S':'Saskatchewan 2024'}
ch=None;Q=[]
for l in T[a:b]:
    m=re.search(r'Chromosome\s+(\d[ABD])',l)
    if m: ch=m.group(1);continue
    if 'Unknow' in l: ch=None
    if not ch: continue
    m=re.match(r'\s*([A-Z]{2,3}\d\d[A-Z](?:,[A-Z]{2,3}\d\d[A-Z])*)\s+(AVRIG\s?\d+)\s+(\d{5,})\s*$',l)
    if m: Q.append(dict(chr=ch,marker=m.group(2).replace(' ',''),bp=int(m.group(3)),traits=m.group(1),array='40K'));continue
    m=re.match(r'\s+(\d{5,})\s+(\S+)\s+([A-Z]{2,3}\d\d[A-Z](?:,\s?[A-Z]{2,3}\d\d[A-Z])*)\s*$',l)
    if m: Q.append(dict(chr=ch,marker=m.group(2),bp=int(m.group(1)),traits=m.group(3).replace(' ',''),array='90K'))
# PVE from S7 and S9 (multi-trait QTN)
pve={}
for l in T:
    m=re.match(r'\s*(\S+(?: \d+)?)?\s+([SGM]\d\d)\s+(INC|DS|IND|FDK|DON|HD|PH)\s+(\d[ABD]|UN)\s+(\d+|-)\s+<0\.001\s+(-?[\d.]+)\s+([\d.]+)',l)
    if m and m.group(4)!='UN' and m.group(3) in FHB: pve.setdefault(int(m.group(5)),[]).append(float(m.group(7)))
recs=[];n=0
for q in Q:
    tr=[t for t in q['traits'].split(',') if re.sub(r'\d\d[A-Z]$','',t) in FHB]
    if not tr: continue
    n+=1
    desc='; '.join(f"{FHB[re.sub(r'\d\d[A-Z]$','',t)]} ({SITE.get(t[-3:],t[-3:])})" for t in tr)
    pv=pve.get(q['bp'])
    recs.append({'id':f'PM6-{n:03d}','name':f"QTN FHB {q['marker']}",'cls':'QTN','disease':'Fusarium head blight','chr':q['chr'],'start':None,'end':None,
      'peak':round(q['bp']/1e6,6),'ref':'RefSeq v1.0','marker':q['marker'],'pve':(f"{min(pv):.2f}-{max(pv):.2f}" if pv and len(set(pv))>1 else (f"{pv[0]:.2f}" if pv else '')),
      'stat':'P < 0.001','source':'Y26_pdis0325','grade':'A','curation':'auto (supplementary table)',
      'evidence':f"Supplementary Table S10 ({q['array']} SNP array): {q['marker']}, chromosome {q['chr']}, {q['bp']:,} bp, traits {q['traits']}",
      'note':f"QTN for {desc} in 300 PGRC hexaploid wheat accessions (multi-locus GWAS, {q['array']} array). Position as published; assembly not named in the supplement, but the 90K positions match the published RefSeq v1.0 coordinates of those markers. Not yet checked against the main paper's 13 QTL clusters."})
json.dump(recs,open('qtn85.json','w'),indent=0)
print(len(Q),len(recs))
