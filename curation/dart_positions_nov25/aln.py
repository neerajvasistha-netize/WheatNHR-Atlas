import mappy,sys,json,os,urllib.request,time,gzip,gc
U='https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.chromosome.%s.fa.gz'
Q=[(h,s) for h,s,_ in mappy.fastx_read('q.fa')]
CH=300_000_000;OV=20_000
for ch in sys.argv[1:]:
    out=f'hits_{ch}.json'
    if os.path.exists(out): continue
    t=time.time(); fn=f'chr{ch}.fa.gz'
    for a in range(4):
        try: urllib.request.urlretrieve(U%ch,fn);break
        except Exception as e: print('retry',ch,e,flush=True);time.sleep(15)
    # split into chunks written as plain fasta
    seq=[];n=0;base=0;chunks=[]
    def flush(s,base):
        cf=f'c_{ch}_{base}.fa'; open(cf,'w').write(f'>{ch}\n{s}\n'); chunks.append((cf,base))
    parts=[];ln=0
    with gzip.open(fn,'rt') as g:
        for line in g:
            if line[0]=='>': continue
            l=line.strip();parts.append(l);ln+=len(l)
            if ln>=CH+OV:
                buf=''.join(parts); flush(buf[:CH+OV],base); base+=CH; buf=buf[CH:]; parts=[buf]; ln=len(buf); del buf
    buf=''.join(parts); parts=None
    if buf: flush(buf,base)
    del buf; os.remove(fn); gc.collect()
    H={}
    for cf,b in chunks:
        A=mappy.Aligner(cf,preset='sr',best_n=3)
        for h,s in Q:
            for x in A.map(s):
                H.setdefault(h,[]).append((x.r_st+b,x.r_en+b,x.strand,x.mlen,x.blen,x.NM,x.mapq,x.q_st,x.q_en))
        del A; gc.collect(); os.remove(cf)
    for h in H:   # drop duplicates from chunk overlaps, keep top 3 by score
        u=sorted(set(map(tuple,H[h])),key=lambda x:-(x[3]-x[5])); H[h]=u[:3]
    json.dump(H,open(out,'w'))
    print(ch,len(H),round(time.time()-t),flush=True)
print('ALLDONE',flush=True)
