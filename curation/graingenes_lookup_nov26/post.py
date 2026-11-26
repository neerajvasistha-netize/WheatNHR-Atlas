import json,glob,collections
P=json.load(open('primers.json'));SQ=json.load(open('seqs.json'))
R={}
for f in glob.glob('res_*.json'):
    ch=f[4:-5];d=json.load(open(f))
    for m,hs in d['epcr'].items():
        for h in hs: R.setdefault(m,{'e':[],'a':[]})['e'].append((ch,)+tuple(h))
    for m,hs in d['aln'].items():
        for h in hs: R.setdefault(m,{'e':[],'a':[]})['a'].append((ch,)+tuple(h))
OUT={}
for m,(f,r,chs) in P.items():
    hs=[h for h in R.get(m,{}).get('e',[]) if h[4]+h[5]<=4]
    on=[h for h in hs if h[0] in chs]
    allc=sorted({h[0] for h in hs})
    if not on: OUT[m]=dict(status='no product on '+'/'.join(chs)+(f' (products on {", ".join(allc)})' if allc else ''),hits=hs);continue
    best=min(h[4]+h[5] for h in on); top=[h for h in on if h[4]+h[5]==best]
    by=collections.defaultdict(list)
    for h in top: by[h[0]].append(h)
    res={}
    for c,L in by.items():
        if max(h[1] for h in L)-min(h[1] for h in L)<2e6: h=min(L,key=lambda x:x[1]); res[c]=(h[1]+1,h[2],h[3],best)
    OUT[m]=dict(status=('unique' if len(on)==1 else f'{len(on)} products on target chromosome(s); best match used'),pos=res,other=sorted(set(allc)-set(chs)),hits=hs)
for m,(s,chs) in SQ.items():
    L=len(s); hs=[h for h in R.get(m,{}).get('a',[]) if (h[8]-h[7])>=0.5*L and h[4]/h[5]>=0.9]
    on=[h for h in hs if h[0] in chs]
    if not on: OUT['SEQ:'+m]=dict(status='no alignment on '+'/'.join(chs)+(f' (hits on {", ".join(sorted({h[0] for h in hs}))})' if hs else ''));continue
    b=max(on,key=lambda h:h[4]-h[6]); res={b[0]:(b[1]+1,b[2],round(b[4]/b[5]*100,1),round((b[8]-b[7])/L*100))}
    OUT['SEQ:'+m]=dict(status='aligned (identity %, coverage %)',pos=res)
json.dump(OUT,open('post.json','w'),default=list)
c=collections.Counter(v['status'].split(' (')[0].split(' on ')[0] for v in OUT.values());print(c)
