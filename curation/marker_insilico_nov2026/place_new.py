"""Third marker placement for still-unplaced loci using new marker position tables from attached supplements
(Crop Journal FHB meta-QTL Table S3, 32k markers on RefSeq v1.0), together with the existing marker table and e-PCR."""
import json,re,sys,csv,glob,collections,statistics
sys.path.insert(0,'pp')
from mp import look
src=open('/home/claude/wheatnhr-atlas/curation/exported_papers/unk.py').read();exec(src[src.index('MK=re.compile'):].split('\n')[0])
EP=json.load(open('genome/epcr_res.json'))
NEW={}
for r in csv.DictReader(open('/home/claude/wheatnhr-atlas/curation/fhb_metaqtl_cropj2021_esm/consensus_marker_positions_refseq_v1.tsv'),delimiter='\t'):
    try: NEW.setdefault(r['marker'].upper(),[]).append((r['chr'],float(r['pos_bp'])/1e6))
    except: pass
CEN={'1A':215,'1B':237,'1D':166,'2A':326,'2B':344,'2D':264,'3A':316,'3B':346,'3D':237,'4A':320,'4B':303,'4D':185,'5A':253,'5B':199,'5D':188,'6A':283,'6B':323,'6D':211,'7A':360,'7B':300,'7D':336}
def keyv(t):
    t=t.upper();t=t[1:] if t.startswith('X') and re.match(r'X(GWM|WMC|BARC|CFD|CFA|GDM|WPT|GPW|CFB|PSP|WMS)',t) else t
    return t.replace('_','-') if t.startswith('WPT') else t
def pos(tok,ch):
    for c in (tok,tok[1:] if tok[:1]=='X' else None):
        if not c: continue
        L=[x for x in look(c,ch) if x[2] in('RefSeq v1.0','RefSeq v2.1')]
        if L: v1=[x for x in L if x[2]=='RefSeq v1.0'];x=(v1 or L)[0];return x[1],x[2],'published position','old'
    e=[x for x in EP.get(tok,[]) if x[0]==ch]
    if e and max(x[1] for x in e)-min(x[1] for x in e)<1: return e[0][1],'RefSeq v2.1','e-PCR','old'
    k=keyv(tok);n=[p for c,p in NEW.get(k,[])+NEW.get(k.replace('-',''),[]) if c==ch]
    if n: return n[0],'RefSeq v1.0','consensus-map marker table (Crop J FHB meta-QTL Table S3)','new'
    return None
def gene_chr(t,name):
    c=collections.Counter()
    for m in re.finditer(re.escape(name),t):
        s=t[max(0,m.start()-150):m.end()+150]
        for x in re.findall(r'\b(?:chromosome|arm|chromosomes)\s+(\d[ABD])',s): c[x]+=1
    return c.most_common(1)[0][0] if c else None
jobs=[]
for h,n,ch in json.load(open('exp/need_fig.json')): jobs.append(('exp',h,n,ch))
for b in ('exp4','exp5'):
    for h,n,w in json.load(open(b+'/s2log.json')): jobs.append((b,h,n,None))
fix={'Lr34':'7D','Yr18':'7D','LrW2':'1D','PmI6':'4A','Sr15':'7A','Pm21':'6A'}
res=[];seen=set();stat=collections.Counter()
for b,h,n,ch in jobs:
    n=n.replace('‐','-').replace('–','-')
    if (h,n) in seen: continue
    seen.add((h,n))
    try: t=re.sub(r'\s+',' ',open(f'{b}/txr/{h}.txt',errors='ignore').read())
    except FileNotFoundError: stat['notext']+=1;continue
    m=re.search(r'-(\d[ABD])([SL])?',n) or re.search(r'(\d[ABD])([SL])?$',n)
    ch=ch or (m.group(1) if m else None) or fix.get(n) or gene_chr(t,n)
    arm=m.group(2) if m else None
    if not ch or ch not in CEN: stat['nochr']+=1;continue
    hits={}
    for mm in re.finditer(re.escape(n)+r'(?![\w.])',t):
        w=t[max(0,mm.start()-150):mm.end()+250]
        for tok in MK.findall(w)+re.findall(r'\b(?:IW[AB]\d{2,6})\b',w):
            p=pos(tok,ch)
            if p and (not arm or ((p[0]<CEN[ch]+40) if arm=='S' else (p[0]>CEN[ch]-40))): hits[tok]=p
    if not any(v[3]=='new' for v in hits.values()): stat['nonew']+=1;continue
    med=statistics.median_low([v[0] for v in hits.values()]);L=[(k,v) for k,v in hits.items() if abs(v[0]-med)<=15]
    gene=re.match(r'(Yr|Lr|Sr|Pm|Fhb|Stb|Tsn|Snn|Sb|Bt|Ltn)\d',n)
    if len({k for k,v in L})<2 and not gene: stat['one']+=1;continue
    if not any(v[3]=='new' for k,v in L): stat['nonew2']+=1;continue
    ps=[v[0] for k,v in L];ref=collections.Counter(v[1] for k,v in L).most_common(1)[0][0]
    stat['placed']+=1
    res.append(dict(batch=b,h=h,name=n,chr=ch,start=round(min(ps),4) if len(ps)>1 else None,end=round(max(ps),4) if len(ps)>1 else None,peak=None if len(ps)>1 else round(ps[0],4),ref=ref,markers=sorted(k for k,v in L),how='; '.join(sorted(f"{k}: {v[2]} {v[0]:.2f} Mb" for k,v in L))))
json.dump(res,open('place_new.json','w'),indent=0);print(stat)
for r in res: print(r['batch'],r['name'],r['chr'],r['start'],r['end'],r['peak'],'|',r['how'][:150])
