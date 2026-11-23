"""Rust QTL (Table S2, marker intervals placed in silico) and GWAS MTAs (Table S7) from Pal et al. 2022 (TAG, rust MDR
meta-QTL; 122_2022_4119_MOESM1_ESM.xlsx) and leaf rust GWAS MTAs (Table S4) from the Plant Genome leaf rust meta-QTL
study (tpg2.20185), credited to the original papers. Meta-QTLs are not added. Only papers identified on Crossref and
not yet in the atlas (DOI) are used.
usage: python3 build81.py rows.json keep.json out.json  (run from the scratchpad, needs pp/mp.py and place_new.py)"""
import json,sys,re,collections
exec(open('place_new.py').read().split('def gene_chr')[0])
R=json.load(open(sys.argv[1]));keep=json.load(open(sys.argv[2]))
def key(ref):
    m=re.match(r"([A-Za-zÀ-ž'\-]+).*?(\d{4})([a-z]?)",str(ref).strip());return f"{m.group(1)}{m.group(2)}{m.group(3)}" if m else None
DIS={'stripe':'Stripe rust','yellow':'Stripe rust','leaf':'Leaf rust','stem':'Stem rust'}
def dis(s):
    s=str(s).lower()
    for k,v in DIS.items():
        if k in s: return v
    return None
def cite(c):
    return {'citation':f"{c[1]} ({c[2]}) {c[3].rstrip('.')}. {c[4]} {c[5]}:{c[6]}".replace(':None','').replace(' None',''),'url':'https://doi.org/'+c[0],
            'short':c[1].split(',')[0].split(' ')[0].capitalize()+(' et al.' if ',' in c[1] else '')+f" ({c[2]})"}
src={};out=[];seen=set();n=0;cnt=collections.Counter();st=collections.Counter()
for r in R['S2']:
    k=key(r[13]);c=keep.get(k)
    if not c: continue
    sid='P2_'+re.sub(r'\W','',c[0])[-12:];src[sid]=cite(c)
    ch=str(r[2]).strip()[:2];iv=str(r[8]).strip();d=dis(r[1])
    if ch not in CEN or not d: st['bad']+=1;continue
    mks=[x.strip() for x in re.split(r'\s*[–\-—/]\s*(?=[A-Za-z])|\s+and\s+|,\s*',iv) if x.strip()]
    P=[(m,pos(m,ch)) for m in mks];P=[(m,p) for m,p in P if p]
    if not P: st['nopos']+=1;continue
    if len(P)>1 and max(p[0] for m,p in P)-min(p[0] for m,p in P)>60: st['far']+=1;continue
    v1=[(m,p) for m,p in P if p[1]=='RefSeq v1.0'];L=v1 or P;rf=L[0][1][1];ps=[p[0] for m,p in L]
    kk=(sid,ch,iv)
    if kk in seen: continue
    seen.add(kk);n+=1;cnt[sid]+=1;two=len(ps)>1
    ab={'Stripe rust':'QYr','Leaf rust':'QLr','Stem rust':'QSr'}[d]
    out.append(dict(id='P2-%04d'%n,name=f"{ab} {src[sid]['short'].split(' ')[0]} {ch} ({iv})",cls='QTL',disease=d,chr=ch,start=round(min(ps),4) if two else None,end=round(max(ps),4) if two else None,peak=None if two else round(ps[0],4),
        ref=rf,marker=', '.join(m for m,p in L),pve=str(r[4] or ''),stat=f"LOD {r[3]}" if r[3] else '',source=sid,grade='C',curation='auto (supplementary compilation, in-silico)',
        evidence=f"Table S2 of Pal et al. 2022: QTL {r[0]}, {r[1]}, interval {iv}, cross {r[9]}",
        note=f"QTL from the original paper ({r[9]}, peak {r[5]} cM), as compiled in Table S2 of Pal et al. (2022, rust multi-disease meta-analysis). Placed in silico at the published physical position of the interval marker(s) {', '.join(f'{m} {p[0]:.2f} Mb ({p[2]})' for m,p in L)} ({rf}). Not yet checked by hand against the original paper."))
def mta(ref,mk,ch,bp,d,tab):
    global n
    k=key(ref);c=keep.get(k)
    if not c or not mk: return
    sid='P2_'+re.sub(r'\W','',c[0])[-12:];src[sid]=cite(c);ch=str(ch).strip()[:2];mk=str(mk).strip()
    if ch not in CEN: return
    kk=(sid,mk,ch)
    if kk in seen: return
    p=pos(mk,ch)
    if p and p[3]=='old': v,rf,how=p[0],p[1],'published physical position of the marker'
    else:
        try: v,rf,how=float(bp)/1e6,'RefSeq v1.0',f'bp position given in {tab}'
        except: return
    seen.add(kk);n+=1;cnt[sid]+=1
    out.append(dict(id='P2-%04d'%n,name=f"{mk} ({d.lower()})",cls='MTA',disease=d,chr=ch,start=None,end=None,peak=round(v,4),ref=rf,marker=mk,pve='',stat='',source=sid,grade='C',
        curation='auto (supplementary compilation)',evidence=f"{tab}: {mk} {ch}",
        note=f"Significant GWAS marker from the original paper, as compiled in {tab}. Position: {how} ({rf}). Not yet checked by hand against the original paper."))
for r in R['S7']:
    d=dis(r[4])
    if d: mta(r[0],r[1],r[2],r[3],d,'Table S7 of Pal et al. (2022)')
for r in R['T4']: mta(r[9],r[2],r[3],r[4],'Leaf rust','Table S4 of the Plant Genome leaf rust meta-QTL study (tpg2.20185)')
json.dump(dict(reports=out,sources=src,counts={src[k]['short']:v for k,v in cnt.items()}),open(sys.argv[3],'w'),indent=0)
print(n,st,collections.Counter(r['cls'] for r in out));print({src[k]['short']:v for k,v in cnt.items()})
