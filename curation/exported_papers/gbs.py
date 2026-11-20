"""Loci whose table rows give GBS/SNP markers named by physical position (e.g. S2B_45213215, chr2B_45213215)."""
import json,re,collections,sys
P={x['h']:x for x in json.load(open('s1.json'))}
need=json.load(open('need_fig.json'))
placed={(r['h'],r['name']) for r in json.load(open('../fig/fig_res.json'))}|{(r['h'],r['name']) for r in json.load(open('../supp/supp_res.json'))}
CEN={'1A':215,'1B':237,'1D':166,'2A':326,'2B':344,'2D':264,'3A':316,'3B':346,'3D':237,'4B':303,'4D':185,'5A':253,'5B':199,'5D':188,'6A':283,'6B':323,'6D':211,'7A':360,'7B':300,'7D':336}
GB=re.compile(r'(?<![A-Za-z0-9])(?:S|chr|Chr)?(\d[ABD])[_:](\d{4,9})(?![\d])')
QP=re.compile(r'Q[A-Z][a-z]{0,5}\.[A-Za-z0-9]{2,10}[-‐–](?:\d[ABD][SL]?(?:[.\-‐]\d{1,2})?)?')
def refof(t):
    tl=t.lower()
    if re.search(r'v\s?2\.1|version 2\.1',tl): return 'RefSeq v2.1'
    if re.search(r'v\s?2\.0|version 2\.0',tl): return 'RefSeq v2.0'
    return 'RefSeq v1.0'
want=collections.defaultdict(set)
for h,n,ch in need:
    n=n.replace('‐','-').replace('–','-')
    if (h,n) not in placed: want[h].add(n)
res=[]
for h,names in want.items():
    t=open('txt/'+h+'.txt',errors='ignore').read().replace('‐','-').replace('–','-');L=t.split('\n');ref=refof(open('txr/'+h+'.txt',errors='ignore').read())
    rows=collections.defaultdict(list)
    cur=None
    for i,l in enumerate(L):
        m=QP.search(l)
        if m:
            nm=m.group(0)
            if nm.endswith('-') and i+1<len(L):
                s=re.match(r'\s*(\d[ABD][SL]?(?:[.\-]\d{1,2})?)\b',L[i+1])
                if s: nm=nm+s.group(1)
            cur=nm if nm in names else None
        elif l.strip()=='' : cur=None
        if cur:
            for c,p in GB.findall(l): rows[cur].append((c,int(p)))
    for nm,ms in rows.items():
        ch=re.search(r'-(\d[ABD])',nm).group(1);a=re.search(r'-\d[ABD]([SL])',nm)
        ps=sorted(p/1e6 for c,p in ms if c==ch)
        if a and ch in CEN:
            c0=CEN[ch]*(1.02 if 'v2' in ref else 1);ps=[x for x in ps if (x<c0+40 if a.group(1)=='S' else x>c0-40)]
        if not ps: continue
        if ps[-1]-ps[0]>150: continue
        res.append(dict(h=h,name=nm,chr=ch,start=round(ps[0],4),end=round(ps[-1],4),ref=ref,markers=sorted({f'S{c}_{p}' for c,p in ms if c==ch})[:6]))
json.dump(res,open('gbs_res.json','w'),indent=0)
print(len(res));[print(r['name'],r['chr'],r['start'],r['end'],r['ref'],r['markers'][:3]) for r in res]
