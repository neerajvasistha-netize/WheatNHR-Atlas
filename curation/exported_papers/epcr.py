import json,gzip,urllib.request,os,re,collections
D=json.load(open('ssr_seqs.json'));S=D['seqs'];C=D['chr']
def rc(s): return s.translate(str.maketrans('ACGTacgtN','TGCAtgcaN'))[::-1]
res=collections.defaultdict(list)
chrs=sorted({c for v in C.values() for c in v})
for ch in chrs:
    f=f'Triticum_aestivum.IWGSC.dna.chromosome.{ch}.fa.gz'
    if not os.path.exists(f):
        urllib.request.urlretrieve('https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/triticum_aestivum/dna/'+f,f)
    with gzip.open(f,'rt') as fh:
        fh.readline();g=''.join(l.strip() for l in fh).upper()
    for k,v in S.items():
        if ch not in C[k]: continue
        hits=[]
        if 'primers' in v:
            p1,p2=v['primers'][0].upper(),v['primers'][1].upper()
            for a,b in [(p1,p2),(p2,p1)]:
                # full primer a on + strand, 3' 14-mer of primer b (as reverse complement) downstream within 1500 bp
                for A,B,lab in [(a,rc(b[-14:]),'primer pair (one exact, one 3prime 14-mer)'),(a[-14:],rc(b),'primer pair (one exact, one 3prime 14-mer)')]:
                    i=g.find(A);n=0
                    while i!=-1 and n<2000:
                        j=g.find(B,i,i+1500)
                        if j!=-1: hits.append((i,lab,'+',j-i))
                        i=g.find(A,i+1);n+=1
            hits=list({h[0]//1000:h for h in hits}.values())
        if not hits and 'sts' in v:
            s=v['sts'].upper();ks=[s[x:x+32] for x in range(20,max(21,len(s)-40),max(1,(len(s)-60)//4))][:5]
            pos=collections.Counter()
            for kmer in ks:
                for q in (kmer,rc(kmer)):
                    i=g.find(q);n=0
                    while i!=-1 and n<5: pos[i//10000]+=1;i=g.find(q,i+1);n+=1
            for b,cnt in pos.items():
                if cnt>=3: hits.append((b*10000,'STS 32-mers',str(cnt),0))
        for h in hits: res[k].append((ch,round(h[0]/1e6,4),h[1],h[2],h[3]))
        print(ch,k,hits[:3],flush=True)
    del g
json.dump(res,open('epcr_res.json','w'),indent=1);print('DONE')
