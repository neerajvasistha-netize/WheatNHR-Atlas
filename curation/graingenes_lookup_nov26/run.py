import mappy,sys,json,os,urllib.request,time,gzip,gc
U='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.chromosome.%s.fa.gz'
P=json.load(open('primers.json'));SQ=json.load(open('seqs.json'))
comp=str.maketrans('ACGTN','TGCAN');rc=lambda s:s.translate(comp)[::-1]
CH=300_000_000;OV=20_000
for ch in sys.argv[1:]:
    out=f'res_{ch}.json'
    if os.path.exists(out): continue
    t=time.time();fn=f'chr{ch}.fa.gz'
    for a in range(4):
        try: urllib.request.urlretrieve(U%ch,fn);break
        except Exception as e: print('retry',ch,e,flush=True);time.sleep(15)
    parts=[]
    with gzip.open(fn,'rt') as g:
        for line in g:
            if line[0]!='>': parts.append(line.rstrip())
    os.remove(fn); S=''.join(parts).upper(); del parts; gc.collect()
    def find(core):
        o=[];i=S.find(core)
        while i>=0: o.append(i);i=S.find(core,i+1)
        return o
    E={}
    for m,(f,r,_) in P.items():
        hits=[]
        for a,b,lab in ((f,r,'F/R'),(r,f,'R/F')):
            for i in find(a[-15:]):
                st=i-(len(a)-15); win=S[i:i+1500]; j=win.find(rc(b)[:15])
                if j>=0:
                    end=i+j+len(b);prod=end-st
                    mmA=sum(1 for x,y in zip(S[st:st+len(a)],a) if x!=y);mmB=sum(1 for x,y in zip(S[end-len(b):end],rc(b)) if x!=y)
                    if 50<=prod<=1500: hits.append((st,end,prod,mmA,mmB,lab))
        if hits: E[m]=hits
    A={}
    for b in range(0,len(S),CH):
        cf=f'c_{ch}.fa'; open(cf,'w').write(f'>{ch}\n{S[b:b+CH+OV]}\n')
        al=mappy.Aligner(cf,preset='sr',best_n=3)
        for m,(s,_) in SQ.items():
            for x in al.map(s): A.setdefault(m,[]).append((x.r_st+b,x.r_en+b,x.strand,x.mlen,x.blen,x.NM,x.q_st,x.q_en))
        del al; os.remove(cf); gc.collect()
    json.dump({'epcr':E,'aln':A},open(out,'w')); del S; gc.collect()
    print(ch,len(E),len(A),round(time.time()-t),flush=True)
print('ALLDONE',flush=True)
