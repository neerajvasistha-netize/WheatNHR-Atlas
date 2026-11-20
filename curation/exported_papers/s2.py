import json,re,collections,sys,statistics
sys.path.insert(0,'../pp')
from mp import look,norm,P
from load2 import load
R,S,base=load()
CHL={'1A':595,'1B':690,'1D':496,'2A':781,'2B':802,'2D':652,'3A':751,'3B':831,'3D':616,'4A':745,'4B':674,'4D':510,'5A':710,'5B':714,'5D':567,'6A':618,'6B':721,'6D':474,'7A':737,'7B':751,'7D':639}
def tdis(code):
    c=code.lower()
    if c.startswith('yr') or c=='pst': return 'Stripe rust'
    if c.startswith('lr'): return 'Leaf rust'
    if c.startswith('sr'): return 'Stem rust'
    if c.startswith('pm'): return 'Powdery mildew'
    if c in('cbt','bt'): return 'Common bunt'
    if c=='kb': return 'Karnal bunt'
    if c=='stb': return 'Septoria tritici blotch'
    return None
GD={'yr':'Stripe rust','lr':'Leaf rust','sr':'Stem rust','pm':'Powdery mildew','ml':'Powdery mildew'}
# existing names
def nk(s): return re.sub(r'[^a-z0-9]','',s.lower().split(' (')[0])
EX=collections.defaultdict(set)
for r in R.values():
    for part in re.split(r'[/,;]| - ',r['name']):
        EX[nk(part)].add(r.get('chr'))
    for g in re.findall(r'\b(?:Yr|Lr|Sr|Pm)\d+[a-z]?\b',r['name']): EX[nk(g)].add(r.get('chr'))
QRX=r'\bQ([A-Z][a-z]{0,5})\.([A-Za-z0-9]{2,10})[-‐–](\d[ABD])([SL])?([.\-‐–]\d{1,2})?\b'
GRX=r'\b(Yr|Lr|Sr|Pm|Ml)(\d{1,3}[a-z]?|[A-Z][A-Za-z0-9]{1,9})\b'
NOTG={'Yr','Lr','Sr','Pm','Ml'}
MBR=re.compile(r'(\d{1,3}(?:\.\d+)?)\s*(?:Mb)?\s*(?:–|-|‐|to|and)\s*(\d{1,3}(?:\.\d+)?)\s*Mb\b')
MB1=re.compile(r'(\d{1,3}\.\d+)\s*Mb\b')
BP=re.compile(r'(?<![\d.,\-_A-Za-z#])(\d{1,3}(?:,\d{3}){1,2}|\d{6,9})(?![\d.,_#A-Za-z])')
TOK=re.compile(r'[A-Za-z][A-Za-z0-9_.\-]{3,40}')
def refof(t):
    tl=t.lower()
    if re.search(r'v\s?2\.1|version 2\.1|refseq v2\b',tl): return 'RefSeq v2.1'
    if re.search(r'v\s?2\.0|version 2\.0',tl): return 'RefSeq v2.0'
    if re.search(r'refseq|iwgsc|chinese spring.{0,30}(genome|assembly)|v\s?1\.0|cs v1',tl): return 'RefSeq v1.0'
    return None
CEN={'1A':215,'1B':237,'1D':166,'2A':326,'2B':344,'2D':264,'3A':316,'3B':346,'3D':237,'4B':303,'4D':185,'5A':253,'5B':199,'5D':188,'6A':283,'6B':323,'6D':211,'7A':360,'7B':300,'7D':336}
def arm(name,ch,s,e,p,ref):
    m=re.search(r'\d[ABD]([SL])\b',name) or re.search(r'-\d[ABD]([SL])',name)
    if not m or ch not in CEN: return None
    x=p if p is not None else (s+e)/2
    c=CEN[ch]*(1.02 if ref and 'v2' in ref else 1)
    return (x<c+40) if m.group(1)=='S' else (x>c-40)
def markers_in(text,ch):
    hits=[]
    for tk in set(TOK.findall(text)):
        tk=tk.strip('.-')
        if len(tk)<4 or re.match(QRX,tk): continue
        L=[x for x in look(tk,ch) if x[2] in('RefSeq v1.0','RefSeq v2.1')]
        if not L and tk.lower().startswith('x'): L=[x for x in look(tk[1:],ch) if x[2] in('RefSeq v1.0','RefSeq v2.1')]
        if L:
            v1=[x for x in L if x[2]=='RefSeq v1.0']; hits.append((tk,(v1 or L)[0]))
    return hits
def cluster(hits):
    if not hits: return None
    ps=sorted(h[1][1] for h in hits); med=statistics.median(ps)
    keep=[h for h in hits if abs(h[1][1]-med)<=25]
    if not keep: return None
    return keep
def ctx(t,name):
    # lines and sentences mentioning the name
    lines=[l for l in t.split('\n') if name in l]
    flat=re.sub(r'\s+',' ',t)
    sents=[s for s in re.split(r'(?<=[.;])\s+(?=[A-Z])',flat) if name in s]
    return lines,sents
papers=json.load(open('s1.json'))
OUT=[];LOG=[]
for pp in papers:
    t=open('txt/'+pp['h']+'.txt',errors='ignore').read()
    t=t.replace(' ',' ').replace(' ',' ')
    try: tr=open('txr/'+pp['h']+'.txt',errors='ignore').read().replace('\u2009',' ').replace('\u00a0',' ')
    except FileNotFoundError: tr=t
    tr=re.sub(r'-\n(?=[a-z])','',tr)
    ref=refof(t)
    pdis=pp['dis'][0][0] if pp['dis'][0][1]>0 else None
    head=pp['ti']+' '+t[:5000]
    loci={}
    # QTL
    qs=collections.Counter(m.group(0) for m in re.finditer(QRX,t))
    codes=collections.Counter()
    for q,n in qs.items(): codes[re.match(QRX,q).group(2).lower()]+=n
    if codes:
        hq={re.match(QRX,m.group(0)).group(2).lower() for m in re.finditer(QRX,re.sub(r'\s+',' ',pp['ti']+' '+tr[:5000]))}
        own=hq or {codes.most_common(1)[0][0]}
        for q,n in qs.items():
            m=re.match(QRX,q); d=tdis(m.group(1))
            if not d or m.group(2).lower() not in own: continue
            loci[q]=dict(cls='QTL',chr=m.group(3),disease=d,n=n)
    # genes in title/abstract
    for m in re.finditer(GRX,pp['ti']):
        g=m.group(0); 
        if g in loci: continue
        d=GD[m.group(1).lower()]
        loci[g]=dict(cls='DG' if re.match(r'\d',m.group(2)) else 'DG',chr=None,disease=d,n=t.count(g))
    rows=[]
    for name,L in loci.items():
        lines,sents=ctx(t,name)
        ch=L['chr']
        if not ch:
            c=collections.Counter()
            for s in sents[:60]:
                for x in re.findall(r'\b(?:chromosome|arm|chromosomes)\s+(\d[ABD])([SL])?',s): c[x[0]]+=2
                for x in re.findall(r'\b(\d[ABD])([SL])\b',s): c[x[0]]+=1
            if c: ch=c.most_common(1)[0][0]
        if not ch: LOG.append((pp['h'],name,'no chromosome'));continue
        s=e=p=None;how='';pref=ref;mk='';ev=''
        # windows after each occurrence: rest of table line (layout) and flat text window
        LW=[];FW=[]
        for m in re.finditer(re.escape(name),t):
            le=t.find('\n',m.end()); LW.append(t[m.end():le if le>0 else m.end()+200][:220])
        flat=re.sub(r'\s+',' ',tr)
        STOP=re.compile(QRX[:-2]+r'|\b(?:Yr|Lr|Sr|Pm)\d|\. [A-Z]')
        for m in re.finditer(re.escape(name),flat):
            w=flat[m.end():m.end()+220]; st=STOP.search(w)
            FW.append(w[:st.start()] if st else w)
        for w in FW[:60]:
            for m in MBR.finditer(w):
                a,b=float(m.group(1)),float(m.group(2))
                if 0<=a<b<=CHL[ch]+5 and b-a<150: s,e=a,b;how='Mb interval in text';ev=name+w[:m.end()];break
            if how: break
        if not how:
            for w in LW[:40]:
                cm=re.search(r'(?<![A-Za-z0-9])(?:chr)?'+ch+r'[SL]?(?![A-Za-z0-9])',w)
                if not cm: continue
                w=w[cm.end():]
                bs=[int(x.replace(',','')) for x in BP.findall(w)]
                bs=[x for x in bs if 1e5<=x<=CHL[ch]*1e6+5e6]
                if bs:
                    if len(bs)>=2 and max(bs)-min(bs)<150e6: s,e=min(bs)/1e6,max(bs)/1e6
                    else: p=bs[0]/1e6
                    how='physical position (bp) in table';ev=name+' ... '+w[:120];break
        if not how:
            hits=cluster(markers_in('\n'.join([STOP.split(w)[0] for w in LW[:30]]+FW[:30]),ch))
            if hits and (len({h[0].lower().lstrip('x') for h in hits})>=2 or L['cls']=='DG'):
                pref=hits[0][1][2]
                hits=[h for h in hits if h[1][2]==pref]; ps=[h[1][1] for h in hits]
                mk=', '.join(sorted({h[0] for h in hits})[:4])
                if len(ps)>1 and max(ps)-min(ps)>0.05: s,e=min(ps),max(ps)
                else: p=ps[0]
                how='in-silico placement of linked markers ('+mk+')'
        if how and arm(name,ch,s,e,p,pref)==False: LOG.append((pp['h'],name,'position on wrong arm'));continue
        if not how: LOG.append((pp['h'],name,'no position'));continue
        if not pref: pref='RefSeq v1.0'
        # pve / lod from table lines
        pve='';stat=''
        for ln in lines[:10]:
            pass
        rows.append(dict(name=name,cls=L['cls'],chr=ch,disease=L['disease'],start=s,end=e,peak=p,ref=pref,marker=mk,how=how,ev=re.sub(r'\s+',' ',ev)[:260],n=L['n'],ctx=(sents[:1] or lines[:1] or [''])[0][:300]))
    OUT.append(dict(h=pp['h'],rows=rows,ref=ref,pdis=pdis))
json.dump(OUT,open('s2.json','w'));json.dump(LOG,open('s2log.json','w'))
n=sum(len(o['rows']) for o in OUT);print('loci placed',n,'papers',sum(1 for o in OUT if o['rows']),'log',collections.Counter(x[2] for x in LOG))
new=0
for o in OUT:
    for r in o['rows']:
        if r['chr'] not in EX.get(nk(r['name']),set()): new+=1
print('new vs atlas',new)
