"""Build the genome-browser data files in web/gb/.
gb/genes_<chr>.json : genes in genes_all.json order: [v1.1 number, v2.1 start (bp), length, strand, exons, cds]
   exons/cds = flat [start,end,...] relative to the gene start (forward strand), from the canonical .1 transcript
   (Ensembl Plants release 63 GFF3, IWGSC RefSeq v1.1 annotation); v2.1 start = v2.1 gene midpoint (liftover cache) - length/2.
gb/cds_<chr>.txt : line 1 JSON {"n":[gene numbers],"o":[offsets],"l":[lengths],"N":[[gene index, pos, run]...]},
   line 2 base64 of the CDS sequences packed 2 bits per base (A0 C1 G2 T3; N stored as A and listed in "N")."""
import gzip,json,re,base64,collections,sys
S='/tmp/claude-0/-home-claude/00583e43-a309-532b-9f1b-7957bcf4ea95/scratchpad/'
W='/home/claude/wheatnhr-atlas/web/'
GFF=S+'gb/dl/Triticum_aestivum.IWGSC.63.gff3.gz'
c=json.load(open('/mnt/user-data/uploads/WheatNHR_liftover/output/_ref_cache.json'))
v21=c['v21']
G=json.load(open(W+'genes_all.json'))
gene={};tx=collections.defaultdict(lambda:{'ex':[],'cds':[]})
for ln in gzip.open(GFF,'rt'):
    if ln[0]=='#':continue
    p=ln.rstrip('\n').split('\t')
    if len(p)<9:continue
    t=p[2]
    if t=='gene':
        m=re.search(r'ID=gene:([^;]+)',p[8])
        if m: gene[m.group(1)]=(p[0],int(p[3]),int(p[4]),p[6])
    elif t in('exon','CDS'):
        m=re.search(r'Parent=transcript:([^;]+)',p[8])
        if m and m.group(1).endswith('.1'):
            tx[m.group(1)[:-2]]['ex' if t=='exon' else 'cds'].append((int(p[3]),int(p[4])))
import os
os.makedirs(W+'gb',exist_ok=True)
CODE={'A':0,'C':1,'G':2,'T':3}
seqs={}
cur=None;buf=[]
for ln in open(S+'pan/cs_query_all.fa'):
    if ln[0]=='>':
        if cur:seqs[cur]=''.join(buf)
        cur=ln[1:].split()[0];buf=[]
    else: buf.append(ln.strip().upper())
if cur:seqs[cur]=''.join(buf)
tot=0;miss=0
for ch,rows in G['g'].items():
    out=[];n=[];o=[];l=[];N=[];packed=bytearray();bits=[]
    seqcat=[]
    for r in rows:
        num=r[1];gid=f'TraesCS{ch}02G{num:06d}'
        g=gene.get(gid)
        if not g: out.append([num,round(r[0]*1e6),0,0,[],[]]);miss+=1;continue
        gs,ge,st=g[1],g[2],g[3];L=ge-gs+1
        mid=v21.get(gid);start21=round(mid[1]*1e6-L/2) if mid and mid[0]==ch else round(r[0]*1e6)
        t=tx.get(gid,{'ex':[],'cds':[]})
        ex=[x for a,b in sorted(t['ex']) for x in (a-gs,b-gs)]
        cd=[x for a,b in sorted(t['cds']) for x in (a-gs,b-gs)]
        out.append([num,start21,L,1 if st=='+' else -1,ex,cd])
        sq=seqs.get(gid)
        if sq:
            n.append(num);o.append(sum(l));l.append(len(sq));seqcat.append(sq)
    s=''.join(seqcat);tot+=len(s)
    for m in re.finditer(r'[^ACGT]+',s):
        # gene index for position
        import bisect
        gi=bisect.bisect_right(o,m.start())-1
        N.append([gi,m.start()-o[gi],m.end()-m.start()])
    s2=re.sub(r'[^ACGT]','A',s)
    b=bytearray((len(s2)+3)//4)
    for i,chh in enumerate(s2):
        b[i>>2]|=CODE[chh]<<((i&3)*2)
    json.dump(out,open(W+f'gb/genes_{ch}.json','w'),separators=(',',':'))
    with open(W+f'gb/cds_{ch}.txt','w') as f:
        f.write(json.dumps({'n':n,'o':o,'l':l,'N':N},separators=(',',':'))+'\n'+base64.b64encode(bytes(b)).decode()+'\n')
    print(ch,len(rows),len(n),len(s))
print('total nt',tot,'genes without GFF',miss)
