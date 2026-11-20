"""Second marker placement for still-unplaced loci: published positions + positions found by e-PCR of SSR primers/STS on RefSeq v2.1."""
import json,re,collections,sys,statistics
sys.path.insert(0,'../pp')
from mp import look
src=open('unk.py').read();exec(src[src.index('MK=re.compile'):].split('\n')[0])
need=json.load(open('need_fig.json'))
done={(r['h'],r['name']) for f in ['../fig/fig_res.json','../supp/supp_res.json','gbs_res.json'] for r in json.load(open(f))}
EP=json.load(open('../genome/epcr_res.json'))
CEN={'1A':215,'1B':237,'1D':166,'2A':326,'2B':344,'2D':264,'3A':316,'3B':346,'3D':237,'4B':303,'4D':185,'5A':253,'5B':199,'5D':188,'6A':283,'6B':323,'6D':211,'7A':360,'7B':300,'7D':336}
def gene_chr(t,name):
    c=collections.Counter();a=collections.Counter()
    for m in re.finditer(re.escape(name),t):
        s=t[max(0,m.start()-150):m.end()+150]
        for x in re.findall(r'\b(?:chromosome|arm|chromosomes)\s+(\d[ABD])([SL])?',s): c[x[0]]+=2
        for x in re.findall(r'\b(\d[ABD])([SL])\b',s): c[x[0]]+=1;a[x]+=1
    if not c: return None,None
    ch=c.most_common(1)[0][0];ar=[k[1] for k,v in a.most_common() if k[0]==ch];return ch,(ar[0] if ar else None)
def pos(tok,ch):
    """position on RefSeq v2.1 or v1.0: (Mb, ref, how)"""
    for c in (tok,tok[1:] if tok[:1]=='X' else None,tok+'_51' if re.fullmatch(r'BS\d{8}',tok) else None):
        if not c: continue
        L=[x for x in look(c,ch) if x[2] in('RefSeq v1.0','RefSeq v2.1')]
        if L: v1=[x for x in L if x[2]=='RefSeq v1.0'];x=(v1 or L)[0];return x[1],x[2],'published position'
    e=[x for x in EP.get(tok,[]) if x[0]==ch]
    if e and max(x[1] for x in e)-min(x[1] for x in e)<1: return e[0][1],'RefSeq v2.1','e-PCR of its primers or STS sequence on RefSeq v2.1'
    return None
res=[];seen=set()
for h,n,ch in need:
    n=n.replace('‐','-').replace('–','-')
    if (h,n) in done or (h,n) in seen: continue
    seen.add((h,n))
    t=re.sub(r'\s+',' ',open('txr/'+h+'.txt',errors='ignore').read())
    g,arm=gene_chr(t,n);ch=ch or g
    m=re.search(r'-\d[ABD]([SL])',n);arm=m.group(1) if m else arm
    if not ch: continue
    hits={}
    for mm in re.finditer(re.escape(n),t):
        w=t[max(0,mm.start()-150):mm.end()+250]
        for tok in MK.findall(w):
            p=pos(tok,ch)
            if p and (not arm or ch not in CEN or ((p[0]<CEN[ch]*1.02+40) if arm=='S' else (p[0]>CEN[ch]-40))): hits[tok]=p
    if not any('e-PCR' in v[2] for v in hits.values()): continue   # only loci newly placeable thanks to e-PCR
    by=collections.defaultdict(list)
    for k,v in hits.items(): by[v[1]].append((k,v))
    ref,L=max(by.items(),key=lambda x:len(x[1]))
    if len(by)>1: # mix of refs: approximate v1.0 vs v2.1 on same scale only if within 15 Mb
        L=[(k,v) for vs in by.values() for k,v in vs]
    med=statistics.median_low([v[0] for k,v in L]);L=[(k,v) for k,v in L if abs(v[0]-med)<=15]
    if len({k for k,v in L})<2 and not re.match(r'(Yr|Lr|Sr|Pm)',n): continue
    ps=[v[0] for k,v in L];ref=collections.Counter(v[1] for k,v in L).most_common(1)[0][0]
    res.append(dict(h=h,name=n,chr=ch,start=round(min(ps),4) if len(ps)>1 else None,end=round(max(ps),4) if len(ps)>1 else None,peak=None if len(ps)>1 else round(ps[0],4),ref=ref,markers=sorted(k for k,v in L),how='; '.join(sorted({f"{k}: {v[2]} {v[0]:.2f} Mb" for k,v in L}))))
json.dump(res,open('epcr_place.json','w'),indent=0);print(len(res))
for r in res: print(r['name'],r['chr'],r['start'],r['end'],r['peak'],r['how'])
