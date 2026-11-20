"""Correct Zhang et al. (2026) TAG 139 (doi 10.1007/s00122-025-05134-0) records against Tables 3-5 of the paper.
The auto-extraction (ad62) had merged the 'Interval' column markers into the marker field and built wide
ranges from them; positions now follow the paper's Position (Mb) column (RefSeq v1.0)."""
import json,glob,copy,sys
W='/home/claude/wheatnhr-atlas/web/'
R={}
for f in ['seed_ad62.json','seed_ad68.json']:
    for r in json.load(open(W+f))['reports']:
        if r.get('source')=='OA_Zhang2026unveil': R[r['id']]=r
SRC='Zhang et al. (2026) Theor Appl Genet 139'
T3='Table 3 (seedling stage, CYR32/CYR34 infection type)'
T4='Table 4 (adult-plant stage, field environments 22YL, 23YL, 24SP, 24TS and BLUP)'
# id: (name, marker, Mb, trait/env, P-value, R2, interval/reported, table)
FIX={
 'OA0356':('QYr.xjau-1A','RAC875_c19509_174',212.02,'CYR34-IT','5.95E-04','3.63','',T3),
 'OA0357':('QYr.xjau-2B.1','BS00066180_51',19.44,'CYR34-IT','6.77E-04','3.90','QYrsnb.nwafu-2BL-1 (AX-95658192 to AX-109400835)',T3),
 'OA0358':('QYr.xjau-2B.3','wsnp_Ra_c54256_57628288',706.73,'CYR32-IT','1.76E-04','3.34','QYrwdbj.nwafu-2B.4 (Xgwm388 to Xwmc317)',T3),
 'OA0359':('QYr.xjau-2B.4','RFL_Contig2324_729',794.89,'CYR34-IT','1.89E-04','4.50','',T3),
 'OA0360':('QYr.xjau-2D.2','Excalibur_c42512_584',577.11,'CYR34-IT','1.28E-04','4.50','',T3),
 'OA0361':('QYr.xjau-4B','IAAV558',588.44,'CYR32-IT','2.55E-04','0.89','',T3),
 'OA0362':('QYr.xjau-5A','wsnp_Ex_c6209_10838852',1.23,'CYR32-IT','1.66E-04','0.89','QYr.cim-5AL (wmc145.1 to wPt-4262)',T3),
 'OA0363':('QYr.xjau-5B.1','wsnp_Ex_c607_1204733',10.44,'CYR32-IT','7.54E-04','3.34','QYr.uga-5B (wPt-665267 to Xgdm152)',T3),
 'OA0364':('QYr.xjau-7A.1','BobWhite_c10229_415',244.47,'CYR34-IT','6.94E-04','3.30','',T3),
 'OA0365':('QYr.xjau-7A.2','wsnp_Ex_c23102_32328851',659.37,'CYR34-IT','2.21E-04','4.27','',T3),
 'OA0366':('QYr.xjau-7A.3','wsnp_Ra_c12773_20367106',731.27,'CYR32-IT','8.11E-04','3.34','',T3),
 'OA0367':('QYr.xjau-2A.1','Tdurum_contig14542_858',69.21,'IT; 22YL/23YL/BLUP','1.14E-04 to 5.01E-04','16.49–20.85','',T4),
 'OA0368':('QYr.xjau-2D.1','RAC875_c8069_1709',5.16,'IT; 22YL/23YL/BLUP','3.71E-06 to 8.89E-04','8.21–14','',T4),
 'OA0369':('QYr.xjau-3A','IACX5824',727.07,'DS&IT (first row); 22YL/23YL/24SP/24TS/BLUP','5.16E-06 to 5.96E-04','7.36–18.25','',T4),
 'OA0370':('QYr.xjau-3D','Kukri_c2289_635',596.72,'DS&IT (first row); 22YL/23YL/24SP/24TS/BLUP','2.22E-05 to 4.94E-04','6.64–14.24','',T4),
 'OA0371':('QYr.xjau-3D','Kukri_c2289_635',596.72,'DS&IT (second row); 22YL/24SP/24TS/BLUP','3.31E-05 to 9.05E-04','8.78–18.56','',T4),
 'OA0372':('QYr.xjau-5B.2','BS00110064_51',397.84,'DS; 24SP/24TS/BLUP','8.69E-05 to 9.77E-04','0.09–2.08','QYr.inra-5BL.1 (Xgwm499 to Xgwm639)',T4),
 'OA0373':('QYr.xjau-6A.1','RAC875_c4196_1689',95.80,'DS; 22YL/23YL/BLUP','5.86E-05 to 8.11E-04','4.17–6.21','QYrhm.nwafu-4B/6A (Xwmc201 to Xgwm570)',T4),
 'OA0374':('QYr.xjau-6A.2','RAC875_c64593_117',26.63,'IT; 22YL/23YL/BLUP','5.47E-05 to 8.48E-04','11.54–15.34','QYrhm.nwafu-4B/6A (Xwmc201 to Xgwm570)',T4),
 'OA0375':('QYr.xjau-6B.1','Kukri_rep_c103613_253',41.99,'DS&IT (first row); 22YL/23YL/24SP/BLUP','8.39E-05 to 9.23E-04','6.90–12.71','',T4),
 'OA0376':('QYr.xjau-6B.1','Kukri_rep_c103613_253',41.99,'DS&IT (second row); 22YL/23YL/BLUP','5.11E-05 to 9.99E-04','10.43–14','',T4),
 'OA0377':('QYr.xjau-6B.2','RAC875_c13920_836',42.78,'DS&IT (first row); 22YL/23YL/24SP/BLUP','4.69E-05 to 8.74E-04','6.94–12.61','',T4),
 'OA0378':('QYr.xjau-6B.2','RAC875_c13920_836',42.78,'DS&IT (second row); 22YL/23YL/24SP/BLUP','4.80E-06 to 8.04E-04','6.15–15.07','',T4),
 'OA0379':('QYr.xjau-6B.3','Tdurum_contig92921_193',44.84,'IT; 22YL/23YL/BLUP','1.01E-05 to 4.21E-04','0.23–3.13','QYrwdbj.nwafu-6B.2 (AX-110145647 to AX-109342799)',T4),
 'OA0380':('QYr.xjau-6B.4','BS00109036_51',663.53,'IT; 22YL/23YL/24SP/BLUP','1.19E-04 to 8.44E-04','0.01–0.1','',T4+'; marker from Table 5'),
}
NEW={'OA1093':('QYr.xjau-2B.1 (adult plant)','BS00067975_51',107.08,'IT; 22YL/23YL/BLUP','1.14E-04 to 5.01E-04','14.97–18.88','QYr.inra-2BL (Xwmc245 to Xwmc175)',T4)}
out=[]
for i,(nm,mk,mb,tr,p,r2,rep,tab) in list(FIX.items())+list(NEW.items()):
    base=copy.deepcopy(R.get(i) or R['OA0357'])
    base.update({'id':i,'name':nm,'marker':mk,'start':None,'end':None,'peak':mb,'ref':'RefSeq v1.0','pve':r2,'stat':'P-value '+p,
      'grade':'B','curation':'hand','evidence':tab,
      'note':f"{tab}: trait {tr}; position {mb} Mb as given in the paper (RefSeq v1.0)."+(f" Previously reported QTL/gene in the region: {rep}." if rep else "")+" Corrected against the paper, 20 Nov 2026 (ad72): the earlier entry had taken the 'Interval' markers as QTL markers."})
    if i=='OA1093': base['chr']='2B'
    base['asm']={}
    out.append(base)
json.dump([{'rec':r} for r in out],open(sys.argv[1],'w'))
print(len(out))

def finish(conv_path,out_path):
    """Keep the earlier pangenome positions where the paper's position matches the old in-silico one
    (within 0.1 Mb); use the fresh conversion for changed or new records."""
    conv={r['id']:r for r in json.load(open(conv_path))};res=[]
    for i,r in conv.items():
        old=R.get(i)
        if old and old.get('peak') is not None and abs(float(old['peak'])-r['peak'])<0.1 and old.get('asm') and 'RefSeq v2.1' in old['asm']:
            r['peak']=old['peak'];r['asm']=old['asm']
        else:
            a={}
            for k,v in r['asm'].items():
                if v.get('flag')=='other_chromosome': continue
                if v.get('flag')=='inverted': v=dict(v,flag='reoriented_checked')
                a[k]=v
            r['asm']=a
        res.append(r)
    json.dump({'reports':res,'sources':{}},open(out_path,'w'),ensure_ascii=False,separators=(',',':'))
    print('wrote',len(res))
if len(sys.argv)>3 and sys.argv[2]=='finish': finish(sys.argv[3],sys.argv[4])
