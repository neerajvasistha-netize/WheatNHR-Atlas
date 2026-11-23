"""Strict pass: a locus name followed (within 80 characters) by a marker interval 'M1-M2' (or 'M1 and M2', 'flanked by M1 and M2')
names that locus' own interval; both markers are looked up (existing table, e-PCR, Crop J Table S3)."""
import json,re,sys,glob,collections
exec(open('place_new.py').read().split('def gene_chr')[0])
MKP=r'(X?(?:gwm|wmc|barc|cfd|cfa|gdm|wms|psp|cfb|gpw|swes|wPt|wpt|rPt|tPt|IW[AB]|BS|AX|Kukri|Excalibur|RAC875|BobWhite|Tdurum|wsnp|Ku|Ex|IAAV|TA|tplb|CAP|GENE|Ra|D_)[A-Za-z0-9_\-]*\d[A-Za-z0-9_]*)'
PAIR=re.compile(MKP+r'\s*(?:[–\-—/]|and|to)\s*'+MKP)
EX=collections.defaultdict(set)
for f in glob.glob('/home/claude/wheatnhr-atlas/web/seed_*.json'):
    for r in json.load(open(f)).get('reports',[]):
        for p in re.split(r'[/,;() ]',r.get('name','')):
            if p: EX[p.lower()].add(r.get('chr'))
jobs=[]
for h,n,ch in json.load(open('exp/need_fig.json')): jobs.append(('exp',h,n,ch))
for b in ('exp4','exp5'):
    for h,n,w in json.load(open(b+'/s2log.json')): jobs.append((b,h,n,None))
fix={'Lr34':'7D','Yr18':'7D','LrW2':'1D','PmI6':'4A','Sr15':'7A','Pm21':'6A','Lr38':'6D'}
res=[];seen=set();st=collections.Counter()
for b,h,n,ch in jobs:
    n=n.replace('‐','-').replace('–','-')
    if (h,n) in seen: continue
    seen.add((h,n))
    if n.lower() in EX: st['atlas']+=1;continue
    try: t=re.sub(r'\s+',' ',open(f'{b}/txr/{h}.txt',errors='ignore').read())
    except FileNotFoundError: continue
    m=re.search(r'-(\d[ABD])([SL])?',n)
    ch=ch or (m.group(1) if m else None) or fix.get(n)
    if not ch or ch not in CEN: st['nochr']+=1;continue
    arm=m.group(2) if m else None
    cand=collections.Counter()
    for mm in re.finditer(re.escape(n)+r'(?![\w.])',t):
        w=t[mm.end():mm.end()+80]
        if re.search(r'\bQ[A-Z][a-z]*\.|\b(Yr|Lr|Sr|Pm)\d',w[:w.find(PAIR.search(w).group(0))] if PAIR.search(w) else ''): continue  # another locus name in between
        p=PAIR.search(w)
        if p: cand[(p.group(1),p.group(2))]+=1
    if not cand: st['nopair']+=1;continue
    (a,c),_=cand.most_common(1)[0]
    P=[(k,pos(k,ch)) for k in (a,c)];P=[(k,v) for k,v in P if v]
    if arm: P=[(k,v) for k,v in P if ((v[0]<CEN[ch]+40) if arm=='S' else (v[0]>CEN[ch]-40))]
    if not P: st['nopos']+=1;continue
    if len(P)==2 and abs(P[0][1][0]-P[1][1][0])>40: st['far']+=1;continue
    ps=[v[0] for k,v in P];st['placed']+=1
    res.append(dict(batch=b,h=h,name=n,chr=ch,interval=f'{a}-{c}',start=round(min(ps),4) if len(ps)>1 else None,end=round(max(ps),4) if len(ps)>1 else None,peak=None if len(ps)>1 else round(ps[0],4),
        ref=collections.Counter(v[1] for k,v in P).most_common(1)[0][0],markers=[k for k,v in P],how='; '.join(f"{k}: {v[2]} {v[0]:.2f} Mb" for k,v in P)))
json.dump(res,open('place_strict.json','w'),indent=0);print(st)
for r in res: print(r['batch'],r['name'],r['chr'],r['interval'],r['start'],r['end'],r['peak'],'|',r['how'][:120])
