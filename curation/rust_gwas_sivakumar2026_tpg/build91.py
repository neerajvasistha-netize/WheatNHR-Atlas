"""Sivakumar T et al. (2026) Plant Genome 10.1002/tpg2.70188: single-trait rust MTAs from Supplemental Table S6
(supplied by the owner on 24 Nov 2026). Positions RefSeq v2.1 (35K Axiom array, as in the paper). AX-94681641 and
AX-94571854 were already added from the paper text (seed_ad84) and are skipped."""
import openpyxl,json,sys
wb=openpyxl.load_workbook(sys.argv[1],read_only=True,data_only=True)
D={'LR':'Leaf rust','SR':'Stem rust','YR':'Stripe rust'}
skip={'AX-94681641','AX-94571854'};recs=[]
for r in list(wb['Sheet1'].iter_rows(values_only=True))[2:]:
    if not r[0]: continue
    mk='AX-'+r[0][2:]; 
    if mk in skip: continue
    envs=[e.strip() for e in r[4].split(',')]; dis=D[envs[0]] if len(envs)==1 else 'Multiple rusts'
    recs.append({'name':f'MTA {"/".join(envs)} {mk}','cls':'MTA','disease':dis,'chr':r[1],'start':None,'end':None,'peak':r[2]/1e6,'ref':'RefSeq v2.1','marker':mk,'pve':'','stat':f'P {r[3]:.2e}',
      'grade':'A','evidence':f'Supplemental Table S6: {r[0]}, {r[1]}, {r[2]:,} bp, P {r[3]:.3g}, trait {r[4]}, models {r[5]}',
      'note':f'Seedling resistance to {" and ".join(D[e].lower() for e in envs)} (single-trait GWAS, {r[5]}) in an Indian wheat panel. Position as published on RefSeq v2.1. Read from Supplemental Table S6 supplied by the owner on 24 Nov 2026.'})
json.dump(recs,open('mta91.json','w'),indent=0);print(len(recs))
