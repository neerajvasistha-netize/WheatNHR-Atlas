import json,glob,sys,re,pandas as pd
sys.path.insert(0,'.');sys.path.insert(0,'../../p26'); from nm import *; from look2 import pos
S=pd.read_excel('../in/122_2004_1740_esm_data.xls',header=None).iloc[1:]
P={n2(r[1]):(str(r[3]).strip(),str(r[4]).strip(),str(r[2]).strip()) for _,r in S.iterrows() if isinstance(r[3],str)}
print('Somers primers',len(P))
need={}
for f in glob.glob('/home/claude/wheatnhr-atlas/web/seed*.json'):
    try: d=json.load(open(f))
    except: continue
    if not isinstance(d,dict): continue
    for r in d.get('reports',[]):
        if r.get('start') is None and r.get('peak') is None and r.get('chr'):
            for m in split_markers(r.get('marker')):
                k=n2(m)
                if k in P and not pos(m,r['chr'][:2]): need.setdefault(k,set()).add(r['chr'][:2])
print(len(need)); json.dump({k:[P[k][0],P[k][1],sorted(v),P[k][2]] for k,v in need.items()},open('somers_need.json','w'),indent=0)
print({k:(sorted(v),P[k][2]) for k,v in need.items()})
