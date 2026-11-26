import json,re,collections,pandas as pd
I=json.load(open('pdfidx.json'));F=json.load(open('first.json'))
d=json.load(open('../cmap/all_now.json'));S=d['sources'];A={r['id']:r for r in d['reports']}
X='/mnt/user-data/attachments/00583e43-a309-532b-9f1b-7957bcf4ea95/2ac25205-WheatTraitMap_unpositioned_records_and_missing_LOD.xlsx'
ids=set(pd.read_excel(X,'Unpositioned records')['Atlas ID'])|set(pd.read_excel(X,'QTL without LOD')['Atlas ID'])
srcs=collections.Counter(A[i]['source'] for i in ids if i in A)
def devpath(p):
    f=I['files'][p]
    return ('Downloads/Paper_2026/'+f.lstrip('./')) if p.startswith('p26txt') else f
M={};how=collections.Counter()
for s in srcs:
    v=S.get(s,{});u=(v.get('url') or '').lower();doi=re.sub(r'^https?://(dx\.)?doi\.org/','',u)
    cand=I['doi'].get(doi) or I['doi'].get(doi.rstrip('.')) if doi else None
    if cand:
        M[s]=devpath(sorted(set(cand),key=lambda p:-len(open(p,errors='ignore').read()))[0]);how['doi']+=1;continue
    cit=v.get('citation','');m=re.search(r'\(\d{4}[a-z]?\)\s*(.+?)\.\s',cit)
    t=re.sub(r'[^a-z0-9 ]',' ',(m.group(1) if m else '').lower());w=[x for x in t.split() if len(x)>3][:10]
    if len(w)>=4:
        best=max(F.items(),key=lambda kv:sum(1 for x in w if x in kv[1]))
        sc=sum(1 for x in w if x in best[1])
        if sc>=max(4,int(len(w)*0.8)): M[s]=devpath(best[0]);how['title']+=1;continue
    how['none']+=1
json.dump(M,open('src2dev.json','w'))
print(len(srcs),how)
