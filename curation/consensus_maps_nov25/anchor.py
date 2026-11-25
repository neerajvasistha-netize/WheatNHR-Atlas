import sys,json,pandas as pd
sys.path.insert(0,'../../p26'); from look2 import pos
M=pd.read_csv('consensus_maps.tsv',sep='\t')
out=[]
for i,r in M.iterrows():
    p=pos(r.marker,r.chr)
    out.append(p[1] if p else None)
M['Mb']=out; M.to_csv('maps_anchored.tsv',sep='\t',index=False)
print(M.groupby('map').Mb.apply(lambda s:s.notna().sum()))
