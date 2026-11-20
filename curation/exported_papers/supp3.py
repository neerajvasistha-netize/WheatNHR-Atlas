"""Loci whose markers are named in the paper text and whose marker positions are listed in the supplement."""
import json,re,collections,statistics,sys,os
exec(open('supp2.py').read().split('res=[];why')[0])
keys=sys.argv[1:]
CH=re.compile(r'^(?:chr)?(\d[ABD])$',re.I)
res=[]
for k in keys:
    h0=keyh.get(k)
    if not h0 or not os.path.exists('rows/'+k+'.json'): continue
    R=json.load(open('rows/'+k+'.json'))
    idx={}
    for r in R:
        cs=[str(c).strip() for c in r[2:]]
        chs=[i for i,c in enumerate(cs) if CH.match(c)]
        if not chs: continue
        ch=CH.match(cs[chs[0]]).group(1)
        nums=[c for c in cs[chs[0]+1:chs[0]+5] if re.fullmatch(r'\d{6,9}(\.0)?',c)]
        if not nums: continue
        v=float(nums[0]);v=v/1e6 if v>1e4 else v
        if v>CHL[ch]+5: continue
        for c in cs[:chs[0]]:
            if len(c)>=4 and not re.fullmatch(r'[\d.]+',c): idx[nz(c)]=(ch,v,c)
    t=re.sub(r'\s+',' ',open(E+'txr/'+h0+'.txt',errors='ignore').read());ref=refof(t)
    for hh,name,ch in need:
        name=name.replace('‐','-').replace('–','-')
        if hh!=h0 or (hh,name) in placed: continue
        g,arm=gene_chr(hh,name);ch=ch or g
        ma=re.search(r'-\d[ABD]([SL])',name);arm=ma.group(1) if ma else arm
        if not ch: continue
        hits={}
        for m in re.finditer(re.escape(name),t):
            w=t[max(0,m.start()-200):m.end()+300]
            for tok in re.findall(r'[A-Za-z][A-Za-z0-9_.\-]{3,40}',w):
                if re.match(r'Q[A-Z][a-z]{0,5}\.',tok): continue
                q=idx.get(nz(tok.strip('.,;')))
                if q and q[0]==ch and armok(arm,ch,q[1],ref): hits[q[2]]=q[1]
        if not hits: continue
        ps=sorted(hits.values());med=statistics.median_low(ps);sel={a:b for a,b in hits.items() if abs(b-med)<=15}
        ps=sorted(sel.values())
        res.append(dict(h=hh,key=k,name=name,chr=ch,start=ps[0] if len(ps)>1 and ps[-1]-ps[0]>0.01 else None,end=ps[-1] if len(ps)>1 and ps[-1]-ps[0]>0.01 else None,peak=None if len(ps)>1 and ps[-1]-ps[0]>0.01 else ps[0],ref=ref,markers=sorted(sel)[:6],how='markers named with the locus in the text; positions from the supplementary marker table',file='supplement',row=', '.join(f'{a} {b:.2f} Mb' for a,b in sel.items())[:240]))
json.dump(res,open('supp3_res.json','w'),indent=0);print(len(res));[print(r['name'],r['chr'],r['start'],r['end'],r['peak'],r['row'][:120]) for r in res]
