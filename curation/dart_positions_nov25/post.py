"""Choose one RefSeq v1.0 position per DArT/DArTseq marker from the per-chromosome minimap2 (sr) alignments."""
import json,glob,re,mappy,collections,csv
L={h.split('|')[0]:(h.split('|')[1],h.split('|')[2],len(s)) for h,s,_ in mappy.fastx_read('q.fa')}
H=collections.defaultdict(list)
for f in glob.glob('hits_*.json'):
    ch=f[5:-5]
    for h,hs in json.load(open(f)).items():
        for x in hs: H[h.split('|')[0]].append((ch,)+tuple(x))
rows=[];st=collections.Counter()
for m,(mch,cm,qlen) in L.items():
    hs=H.get(m,[])
    # score: matched bases minus mismatches/gaps
    sc=[(x[4]-x[6],x) for x in hs]
    f=1 if qlen<150 else 0.6
    sc=[(s,x) for s,x in sc if ((x[9]-x[8])>=0.8*qlen*f and x[4]/x[5]>=0.9) or ((x[9]-x[8])>=0.6*qlen*f and x[4]/x[5]>=0.95)]   # coverage and identity
    sc.sort(key=lambda t:-t[0])
    if not sc: rows.append([m,mch,cm,qlen,'','','','','','',0,'no hit']);st['no hit']+=1;continue
    best=sc[0][0]; top=[t for t in sc if t[0]>=best-2]          # equally good copies
    on=[t for t in top if t[1][0]==mch]
    nchr=len({t[1][0] for t in top})
    if len(top)==1:
        x=top[0][1]; status='unique' + ('' if mch in ('NA',x[0]) else ' (other chromosome than the map)')
    elif on and len(on)==1:
        x=on[0][1]; status=f'map chromosome copy chosen ({len(top)} equal hits on {nchr} chromosomes)'
    elif len({t[1][0] for t in top})==1 and max(t[1][1] for t in top)-min(t[1][1] for t in top)<2e6:
        x=top[0][1]; status='unique region (tandem copies within 2 Mb)'
    else:
        rows.append([m,mch,cm,qlen,'','','','','','',len(top),f'ambiguous ({len(top)} equal hits on {nchr} chromosomes)']);st['ambiguous']+=1;continue
    ch,rs,re_,strand,mlen,blen,nm,mq,qs,qe=x
    ident=round(mlen/blen*100,1)
    cov=round((qe-qs)/qlen*100)
    if cov<80 and qlen<150: status+=' (partial tag alignment)'
    rows.append([m,mch,cm,qlen,ch,rs+1,re_,'+' if strand>0 else '-',ident,blen,len(top),status]);st[status.split(' (')[0]]+=1
json.dump(rows,open('dart_positions.json','w'))
print(len(rows),st.most_common())
