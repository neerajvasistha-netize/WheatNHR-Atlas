"""Relaxed in-silico PCR for markers without a product on their target chromosome: seeds on the 3' 10 bases of each primer,
<=3 mismatches over each whole primer, the two primers facing each other 50-1500 bp apart."""
import gzip,urllib.request,json,os,sys,time,collections
U='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.chromosome.%s.fa.gz'
P=json.load(open('primers.json'));O=json.load(open('post.json'))
need=collections.defaultdict(list)
for m,v in O.items():
    if m.startswith('SEQ:') or v.get('pos'): continue
    for c in P[m][2]:
        if len(c)==2: need[c].append(m)
comp=str.maketrans('ACGTN','TGCAN');rc=lambda s:s.translate(comp)[::-1]
for ch,ms in sorted(need.items()):
    out=f'rel_{ch}.json'
    if os.path.exists(out): continue
    t=time.time();fn=f'r{ch}.fa.gz'
    for a in range(4):
        try: urllib.request.urlretrieve(U%ch,fn);break
        except Exception as e: time.sleep(15)
    S=''.join(l.rstrip() for l in gzip.open(fn,'rt') if l[0]!='>').upper(); os.remove(fn)
    def sites(p,fwd):
        core=p[-10:] if fwd else rc(p)[:10]; q=p if fwd else rc(p); L=[];i=S.find(core)
        while i>=0:
            st=i-(len(p)-10) if fwd else i
            mm=sum(1 for a,b in zip(S[st:st+len(p)],q) if a!=b)
            if mm<=3: L.append((st,st+len(p),mm))
            i=S.find(core,i+1)
        return L
    R={}
    for m in ms:
        f,r,_=P[m]; hits=[]
        for a,b in ((f,r),(r,f)):
            A=sites(a,True); B=sites(b,False)
            for s1,e1,m1 in A:
                for s2,e2,m2 in B:
                    if 50<=e2-s1<=1500: hits.append((s1,e2,e2-s1,m1,m2))
        R[m]=hits
    json.dump(R,open(out,'w')); del S
    print(ch,{m:len(v) for m,v in R.items()},round(time.time()-t),flush=True)
print('ALLDONE',flush=True)
