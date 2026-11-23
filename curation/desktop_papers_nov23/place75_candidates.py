import json,re,collections,sys,statistics
sys.path.insert(0,'../pp')
from mp import look
src=open('/home/claude/wheatnhr-atlas/curation/exported_papers/unk.py').read();exec(src[src.index('MK=re.compile'):].split('\n')[0])
EP=json.load(open('../genome/epcr_res.json'))
CEN={'1A':215,'1B':237,'1D':166,'2A':326,'2B':344,'2D':264,'3A':316,'3B':346,'3D':237,'4A':320,'4B':303,'4D':185,'5A':253,'5B':199,'5D':188,'6A':283,'6B':323,'6D':211,'7A':360,'7B':300,'7D':336}
fix={'Lr34':'7D','Yr18':'7D','LrW2':'1D','PmI6':'4A','Sr15':'7A','Pm21':'6A'}
CH={nm:fix.get(nm,c) for nm,c,s in json.load(open('unpl_chr.json'))}
def pos(tok,ch):
    for c in (tok,tok[1:] if tok[:1]=='X' else None):
        if not c: continue
        L=[x for x in look(c,ch) if x[2] in('RefSeq v1.0','RefSeq v2.1')]
        if L: v1=[x for x in L if x[2]=='RefSeq v1.0'];x=(v1 or L)[0];return x[1],x[2],'published position'
    e=[x for x in EP.get(tok,[])+EP.get(tok[1:] if tok[:1]=='X' else '',[]) if x[0]==ch]
    if e and max(x[1] for x in e)-min(x[1] for x in e)<1: return e[0][1],'RefSeq v2.1','e-PCR on RefSeq v2.1'
    return None
res=[];nomk=[];unknown=collections.defaultdict(set);seen=set()
for h,n,why in json.load(open('s2log.json')):
    n=n.replace('‐','-')
    if (h,n) in seen: continue
    seen.add((h,n));ch=CH.get(n)
    t=re.sub(r'\s+',' ',open('txt/'+h+'.txt',errors='ignore').read())
    m=re.search(r'-\d[ABD]([SL])',n);arm=m.group(1) if m else None
    hits={};toks=set()
    for mm in re.finditer(re.escape(n)+r'(?![\w.])',t):
        w=t[max(0,mm.start()-150):mm.end()+250]
        for tok in MK.findall(w):
            toks.add(tok);p=pos(tok,ch)
            if p and (not arm or ((p[0]<CEN[ch]+40) if arm=='S' else (p[0]>CEN[ch]-40))): hits[tok]=p
            elif not p: unknown[tok].add(ch)
    if not hits: nomk.append((n,ch,sorted(toks)));continue
    med=statistics.median_low([v[0] for v in hits.values()]);L=[(k,v) for k,v in hits.items() if abs(v[0]-med)<=15]
    ps=[v[0] for k,v in L];ref=collections.Counter(v[1] for k,v in L).most_common(1)[0][0]
    res.append(dict(h=h,name=n,chr=ch,start=round(min(ps),4) if len(ps)>1 else None,end=round(max(ps),4) if len(ps)>1 else None,peak=None if len(ps)>1 else round(ps[0],4),ref=ref,markers=sorted(k for k,v in L),how='; '.join(sorted(f"{k}: {v[2]} {v[0]:.2f} Mb" for k,v in L))))
json.dump(dict(res=res,nomk=nomk,unknown={k:sorted(v) for k,v in unknown.items()}),open('place75.json','w'),indent=0)
for r in res: print('PLACED',r['name'],r['chr'],r['start'],r['end'],r['peak'],'|',r['how'])
for x in nomk: print('NO',x)
print(len(unknown),'unknown markers')
