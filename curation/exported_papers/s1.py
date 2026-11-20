import json,re,collections,os
T=json.load(open('../exported_triage.json'))
new={x[0]:x for x in T['new']}
idx=dict(l.rstrip('\n').split('\t')[::-1] for l in open('txt/index.tsv'))  # path->hash
DIS=[('Stripe rust',r'stripe rust|yellow rust|Puccinia striiformis|\bYr\d'),('Leaf rust',r'leaf rust|brown rust|Puccinia triticina|\bLr\d'),
('Stem rust',r'stem rust|black rust|Puccinia graminis|\bSr\d|Ug99'),('Powdery mildew',r'powdery mildew|Blumeria|\bPm\d'),('Crown rust',r'crown rust'),
('Fusarium head blight',r'Fusarium head blight|head scab'),('Septoria tritici blotch',r'Septoria tritici|Zymoseptoria'),('Spot blotch',r'spot blotch|Bipolaris'),('Tan spot',r'tan spot|Pyrenophora tritici')]
out=[]
for p,(path,au,yr,ti) in new.items():
    h=idx.get(p)
    if not h: print('miss',p);continue
    t=open('txt/'+h+'.txt',errors='ignore').read()
    doi=re.search(r'(?:doi\.org/|doi:\s*|DOI\s*)(10\.\d{4,9}/[^\s,;]+)',t,re.I)
    doi=doi.group(1).rstrip('.)') if doi else ''
    head=ti+' '+t[:6000]
    sc={d:len(re.findall(rx,head,re.I))*3+len(re.findall(rx,t,re.I)) for d,rx in DIS}
    qtl=collections.Counter(re.findall(r'\bQ[A-Z][a-z]{0,5}\.[A-Za-z0-9]{2,10}[-‐–]\d[ABD](?:[SL])?(?:\.\d{1,2})?\b',t))
    gen=collections.Counter(re.findall(r'\b(?:Yr|Lr|Sr|Pm)(?:\d{1,3}[a-z]?|[A-Z][A-Za-z0-9]{1,8})\b',head))
    out.append(dict(h=h,path=path,au=au,yr=yr,ti=ti,doi=doi,n=len(t),dis=sorted(sc.items(),key=lambda x:-x[1])[:3],qtl=qtl.most_common(40),gen=gen.most_common(12)))
json.dump(out,open('s1.json','w'))
print(len(out),sum(1 for o in out if o['qtl']),sum(1 for o in out if o['doi']),sum(1 for o in out if o['n']<3000))
