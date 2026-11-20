import json,re,csv,sys,collections,statistics,glob,os
sys.path.insert(0,'../pp')
from mp import look
E='../exp/'
CHL={'1A':595,'1B':690,'1D':496,'2A':781,'2B':802,'2D':652,'3A':751,'3B':831,'3D':616,'4A':745,'4B':674,'4D':510,'5A':710,'5B':714,'5D':567,'6A':618,'6B':721,'6D':474,'7A':737,'7B':751,'7D':639}
CEN={'1A':215,'1B':237,'1D':166,'2A':326,'2B':344,'2D':264,'3A':316,'3B':346,'3D':237,'4B':303,'4D':185,'5A':253,'5B':199,'5D':188,'6A':283,'6B':323,'6D':211,'7A':360,'7B':300,'7D':336}
def refof(t):
    tl=t.lower()
    if re.search(r'v\s?2\.1|version 2\.1|refseq v2\b',tl): return 'RefSeq v2.1'
    if re.search(r'v\s?2\.0|version 2\.0',tl): return 'RefSeq v2.0'
    return 'RefSeq v1.0'
def nz(s): return re.sub(r'[^a-z0-9]','',str(s).lower())
P={x['h']:x for x in json.load(open(E+'s1.json'))};meta=json.load(open(E+'rdfmeta.json'))
need=json.load(open(E+'need_fig.json'));placed={(r['h'],r['name']) for r in json.load(open('../fig/fig_res.json'))}
# key -> h
keyh={}
for r in csv.DictReader(open(E+'supplementary_needed.tsv'),delimiter='\t'):
    doi=r['doi_link'].replace('https://doi.org/','')
    for h,pp in P.items():
        if meta[pp['path']]['doi']==doi: keyh[r['key']]=h
exec(open('../fig/fig2.py').read().split('def caption')[0].split("def gene_chr")[1].join(["def gene_chr",""]) if False else '')
def gene_chr(h,name):
    t=open(E+'txr/'+h+'.txt',errors='ignore').read().replace('\n',' ');c=collections.Counter();a=collections.Counter()
    for m in re.finditer(re.escape(name),t):
        s=t[max(0,m.start()-150):m.end()+150]
        for x in re.findall(r'\b(?:chromosome|arm|chromosomes)\s+(\d[ABD])([SL])?',s): c[x[0]]+=2
        for x in re.findall(r'\b(\d[ABD])([SL])\b',s): c[x[0]]+=1;a[x]+=1
    if not c: return None,None
    ch=c.most_common(1)[0][0];arms=[k[1] for k,v in a.most_common() if k[0]==ch]
    return ch,(arms[0] if arms else None)
def armok(a,ch,x,ref):
    if not a or ch not in CEN: return True
    c=CEN[ch]*(1.02 if 'v2' in ref else 1);return x<c+40 if a=='S' else x>c-40
def mk(tx,ch):
    tx=str(tx).strip().strip('.,;:()[]*')
    if len(tx)<4 or re.fullmatch(r'[\d.,]+',tx): return None
    for c in (tx,tx[1:] if tx[:1] in 'Xx' else None,re.sub(r'^K[-_]','',tx)):
        if c:
            L=[x for x in look(c,ch) if x[2] in('RefSeq v1.0','RefSeq v2.1')]
            if L: v1=[x for x in L if x[2]=='RefSeq v1.0'];return (v1 or L)[0]
    return None
res=[];why=collections.Counter()
seen=set()
for h0,name,ch in need:
    name=name.replace('‐','-').replace('–','-')
    if (h0,name) in placed or (h0,name) in seen: continue
    seen.add((h0,name))
    keys=[k for k,h in keyh.items() if h==h0]
    if not keys or not os.path.exists('rows/'+keys[0]+'.json'): why['no supplementary file']+=1;continue
    R=json.load(open('rows/'+keys[0]+'.json'))
    if not R: why['supplement unreadable or empty']+=1;continue
    g,arm=gene_chr(h0,name)
    m=re.search(r'-\d[ABD]([SL])',name)
    if m: arm=m.group(1)
    ch=ch or g
    if not ch: why['no chromosome']+=1;continue
    ref=refof(open(E+'txr/'+h0+'.txt',errors='ignore').read())
    hit=[r for r in R if any(nz(name)==nz(c) or (nz(name) in nz(c) and len(nz(c))<=len(nz(name))+3) for c in r[2:])]
    if not hit: why['locus not named in supplement']+=1;continue
    s=e=p=None;how='';mks=[]
    bps=[]
    for r in hit:
        cs=[str(c).strip() for c in r[2:]]
        for i,c in enumerate(cs):
            m2=re.fullmatch(r'(?:chr)?(\d[ABD])[SL]?\s*[:_]\s*(\d{1,3}(?:,\d{3}){1,2}|\d{5,9})(?:\s*[-–]\s*(\d{1,3}(?:,\d{3}){1,2}|\d{5,9}))?',c)
            vals=[]
            if m2:
                if m2.group(1)!=ch: continue
                vals=[m2.group(2)]+([m2.group(3)] if m2.group(3) else [])
            elif re.fullmatch(r'\d{1,3}(?:,\d{3}){1,2}|\d{6,9}(?:\.0)?',c) and any(re.fullmatch(r'(?:chr)?'+ch+r'[SL]?',x,re.I) for x in cs[max(0,i-2):i]):
                vals=[c.split('.')[0]]
            for x in vals:
                v=int(x.replace(',',''))
                if 1e5<=v<=CHL[ch]*1e6+5e6 and armok(arm,ch,v/1e6,ref): bps.append(v/1e6)
        for c in r[2:]:
            for tok in re.split(r'[\s,;/]+',c):
                q=mk(tok,ch)
                if q and armok(arm,ch,q[1],q[2]): mks.append((tok.strip('.,;:()*'),q))
    if bps:
        if max(bps)-min(bps)>250: bps=[]
    if bps:
        s,e=min(bps),max(bps);how='physical position in the supplementary table row of the locus';mk_s=''
    elif len({nz(t) for t,q in mks})>=1 and mks:
        rf=collections.Counter(q[2] for t,q in mks).most_common(1)[0][0];mks=[(t,q) for t,q in mks if q[2]==rf]
        med=statistics.median_low([q[1] for t,q in mks]);mks=[(t,q) for t,q in mks if abs(q[1]-med)<=15]
        ps=[q[1] for t,q in mks];s,e=min(ps),max(ps);ref=rf
        how='in-silico placement of the markers in the supplementary table row of the locus'
    else: why['no position or placeable marker in its row']+=1;continue
    if e-s<0.01: p,s,e=s,None,None
    res.append(dict(h=h0,key=keys[0],name=name,chr=ch,start=s,end=e,peak=p,ref=ref,markers=sorted({t for t,q in mks})[:6],how=how,file=hit[0][0],row=' | '.join(hit[0][1:])[:240]))
    why['placed']+=1
json.dump(res,open('supp_res.json','w'),indent=0);print(why)
