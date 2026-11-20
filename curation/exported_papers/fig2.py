import json,re,glob,csv,collections,statistics,subprocess,sys,hashlib,os
sys.path.insert(0,'../pp');sys.path.insert(0,'../exp')
from mp import look
E='../exp/'
P={x['h']:x for x in json.load(open(E+'s1.json'))}
need=json.load(open(E+'need_fig.json'))
BASE='/mnt/user-data/uploads/Downloads/Exported /'
QRX=r'\bQ([A-Z][a-z]{0,5})\.([A-Za-z0-9]{2,10})[-‐–](\d[ABD])'
def nz(s): return re.sub(r'[^a-z0-9]','',s.lower())
def words(f):
    k=hashlib.md5(f.encode()).hexdigest()[:10];W=collections.defaultdict(list);dims={}
    for t in sorted(glob.glob('ocr/'+k+'_*.tsv')):
        pg=int(t[-7:-4])
        for r in csv.DictReader(open(t),delimiter='\t',quoting=csv.QUOTE_NONE):
            if r['level']=='1': dims[pg]=(int(r['width']),int(r['height']))
            tx=(r.get('text') or '').strip()
            if tx and float(r['conf'] or -1)>30: W[pg].append((tx,int(r['left'])+int(r['width'])/2,int(r['top'])+int(r['height'])/2))
    return W,dims
def mk(tx,ch):
    tx=tx.strip('.,;:()[]*†‡')
    if len(tx)<4 or tx.isdigit(): return None
    for c in (tx,tx[1:] if tx[:1] in 'Xx' else None,re.sub(r'^K[-_]','',tx)):
        if not c: continue
        L=[x for x in look(c,ch) if x[2] in('RefSeq v1.0','RefSeq v2.1')]
        if L:
            v1=[x for x in L if x[2]=='RefSeq v1.0'];return (v1 or L)[0]
    return None
def gene_chr(h,name):
    t=open(E+'txr/'+h+'.txt',errors='ignore').read().replace('\n',' ')
    c=collections.Counter()
    for m in re.finditer(re.escape(name),t):
        s=t[max(0,m.start()-150):m.end()+150]
        for x in re.findall(r'\b(?:chromosome|arm|chromosomes)\s+(\d[ABD])([SL])?',s): c[x[0]]+=2
        for x in re.findall(r'\b(\d[ABD])([SL])\b',s): c[x[0]]+=1
    a=collections.Counter()
    for m in re.finditer(re.escape(name),t):
        s=t[max(0,m.start()-150):m.end()+150]
        for x in re.findall(r'\b(\d[ABD])([SL])\b',s): a[x]+=1
    if not c: return None
    ch=c.most_common(1)[0][0];arms=[k[1] for k,v in a.most_common() if k[0]==ch]
    ARM[name]=arms[0] if arms else None
    return ch
def caption(f,pg):
    t=subprocess.run(['pdftotext','-f',str(pg),'-l',str(pg),f,'-'],capture_output=True,text=True).stdout
    m=re.search(r'(Fig\.?|Figure)\s*(\d+)',t);return (m.group(0) if m else 'figure on page %d'%pg),t
res=[];why=collections.Counter();cache={};ARM={};seen=set()
CEN={'1A':215,'1B':237,'1D':166,'2A':326,'2B':344,'2D':264,'3A':316,'3B':346,'3D':237,'4B':303,'4D':185,'5A':253,'5B':199,'5D':188,'6A':283,'6B':323,'6D':211,'7A':360,'7B':300,'7D':336}
def armok(name,ch,x,ref):
    m=re.search(r'-\d[ABD]([SL])',name);a=m.group(1) if m else ARM.get(name)
    if not a or ch not in CEN: return True
    c=CEN[ch]*(1.02 if 'v2' in ref else 1)
    return x<c+40 if a=='S' else x>c-40
for h,name,ch in need:
    name=name.replace('‐','-').replace('–','-')
    if (h,name) in seen: continue
    seen.add((h,name))
    ch0=ch
    g=gene_chr(h,name)
    if not ch: ch=g
    elif g!=ch: ARM.pop(name,None)
    if not ch: why['no chromosome']+=1;continue
    f=BASE+P[h]['path'][2:]
    if f not in cache: cache[f]=words(f)
    W,dims=cache[f]
    best=None
    for pg,ws in W.items():
        wd,ht=dims.get(pg,(2500,3300))
        lab=[w for w in ws if nz(name) in nz(w[0]) and len(nz(w[0]))<=len(nz(name))+4]
        if not lab:  # label split across OCR tokens, e.g. "QYr.wpg" "-1B.2"
            tail=nz(name.split('-')[-1]) if '-' in name else None
            if tail: lab=[w for w in ws if nz(w[0]) in (tail,tail.replace('.',''))]
            if lab and not any(nz(name.split('-')[0])==nz(w[0]) for w in ws): lab=[]
        cand=[(w,mk(w[0],ch)) for w in ws];cand=[(w,p) for w,p in cand if p and armok(name,ch,p[1],p[2])]
        if not cand: continue
        sel=[];how=''
        if lab:
            for lw in lab:
                near=[(w,p) for w,p in cand if abs(w[1]-lw[1])<0.22*wd and abs(w[2]-lw[2])<0.25*ht]
                near=sorted(near,key=lambda q:((q[0][1]-lw[1])/wd)**2+(2*(q[0][2]-lw[2])/ht)**2)[:6]
                if near:
                    md=statistics.median(p[1] for w,p in near[:3]);near=[(w,p) for w,p in near if abs(p[1]-md)<=15]
                if len({nz(w[0]) for w,p in near})>len({nz(w[0]) for w,p in sel}): sel=near
            how='markers next to the locus label in '
        else:
            cap,t=caption(f,pg)
            if name in t and len(cand)>=2: sel=cand;how='markers in '
        if len({nz(w[0]) for w,p in sel})<2: continue
        ref=collections.Counter(p[2] for w,p in sel).most_common(1)[0][0]
        sel=[(w,p) for w,p in sel if p[2]==ref]
        med=statistics.median(p[1] for w,p in sel);sel=[(w,p) for w,p in sel if abs(p[1]-med)<=15]
        if len({nz(w[0]) for w,p in sel})<2: continue
        if best is None or (bool(lab),len(sel))>(best[4],len(best[1])): best=(pg,sel,how,ref,bool(lab))
    if not best: why['no markers found in figures']+=1;continue
    pg,sel,how,ref,_l=best;cap,_=caption(f,pg);ps=[p[1] for w,p in sel]
    if not armok(name,ch,(min(ps)+max(ps))/2,ref): why['wrong arm']+=1;continue
    res.append(dict(h=h,name=name,chr=ch,start=round(min(ps),4),end=round(max(ps),4),ref=ref,markers=sorted({w[0].strip('.,;:()') for w,p in sel}),fig=cap,how=how+cap))
    why['placed']+=1
json.dump(res,open('fig_res.json','w'));print(why)
