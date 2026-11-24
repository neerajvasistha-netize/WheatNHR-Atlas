"""Wang J et al. (2026) TAG 10.1007/s00122-026-05312-8: integrated yellow rust QTL from Supplemental Table 4 (ESM2).
One record per integrated QTL (MEQ = multi-environment, SEQ = single-environment). QTL shaded purple in the table
(omitted by the authors for linkage disequilibrium) and ChrUn QTL are left out, as are the 5 MEQ already added in
seed_ad88 from the paper text. Positions are Chinese Spring RefSeq v2.1 (stated in Supplemental Table 6)."""
import openpyxl,json,sys,re,collections
wb=openpyxl.load_workbook(sys.argv[1]); ws=wb['Supplemental Table 4']
have={'Q1A_3.92_10.73','Q2A_0.47_39.79','Q2B_112.61_207.39','Q5D_504.52_508.21','Q6A_607.21_617.66'}
Q=collections.OrderedDict(); purple=collections.Counter(); tot=collections.Counter()
for row in ws.iter_rows(min_row=4):
    v=[c.value for c in row]
    if not v[1]: continue
    qid=v[1]; tot[qid]+=1
    if row[1].fill and row[1].fill.fill_type and row[1].fill.fgColor.rgb=='FF7030A0': purple[qid]+=1; continue
    q=Q.setdefault(qid,{'chr':v[0],'cat':v[16],'rows':[]}); q['rows'].append(v)
recs=[];skipped=[]
for qid,q in Q.items():
    ch=q['chr'].replace('Chr','')
    if not re.match(r'^[1-7][ABD]$',ch) or qid in have: skipped.append(qid); continue
    m=re.match(r'Q\w+?_([\d.]+)(?:_([\d.]+))?$',qid); s=float(m.group(1)); e=float(m.group(2)) if m.group(2) else s
    best=max(q['rows'],key=lambda r:r[9] or 0)
    envs=sorted({str(r[20]) for r in q['rows'] if r[20]}); exps=sorted({r[19] for r in q['rows'] if r[19]})
    pves=[r[10] for r in q['rows'] if isinstance(r[10],(int,float))]
    lead=int(best[8].split('_')[1])/1e6
    meq='multi-environment QTL (MEQ)' if q['cat']=='MEQ' else 'single-environment QTL (SEQ)'
    recs.append({'name':qid,'cls':'QTL','disease':'Stripe rust','chr':ch,'start':s if e>s else None,'end':e if e>s else None,'peak':lead,
      'ref':'RefSeq v2.1','marker':best[8],'pve':f'{max(pves)*100:.1f}' if pves else '','stat':f'-log10(P) {best[9]:.2f}',
      'grade':'B' if e>s else 'A',
      'evidence':f'Supplemental Table 4: {qid}, {q["chr"]}, {q["cat"]}, lead SNP {best[8]} (-log10 P {best[9]:.2f}), {best[14]} model; {", ".join(exps)}; environments {", ".join(envs)}',
      'note':f'Yellow rust {meq} from GWAS of 5,840 hybrids and elite lines in central Europe; {len(q["rows"])} raw detections ({", ".join(exps)}). Interval from the integrated QTL ID (Mb, RefSeq v2.1); peak is the lead SNP. Dominance class: {best[15]}.'})
json.dump(recs,open('qtl89.json','w'),indent=0)
print(len(Q),len(recs),'purple-only',len([k for k in tot if purple[k]==tot[k]]),'skipped',skipped)
print(collections.Counter(r['note'].split(' QTL')[0][-16:] for r in recs))
