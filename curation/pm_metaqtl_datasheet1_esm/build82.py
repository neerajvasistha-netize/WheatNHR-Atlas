"""Powdery mildew QTL (Table S8, marker intervals placed in silico) and GWAS MTAs (Table S6) from the supplementary
Data Sheet 1 of a wheat powdery mildew meta-QTL study, credited to the original papers. Meta-QTLs are not added.
Only papers identified on Crossref (keep.json, automatic matches plus hand-checked DOIs) and not yet in the atlas.
usage (from the scratchpad): python3 build82.py rows.json keep.json out.json"""
import json,sys,re,collections
exec(open('place_new.py').read().split('def gene_chr')[0])
R=json.load(open(sys.argv[1]));keep=json.load(open(sys.argv[2]))
def key(ref):
    m=re.match(r"([A-Za-zÀ-ž'\-]+).*?(\d{4})([a-z]?)",str(ref).strip());return f"{m.group(1)}{m.group(2)}{m.group(3)}" if m else None
def cite(c):
    return {'citation':f"{c[1]} ({c[2]}) {c[3].rstrip('.')}. {c[4]} {c[5]}:{c[6]}".replace(':None','').replace(' None',''),'url':'https://doi.org/'+c[0],
            'short':c[1].split(',')[0].split(' ')[0].capitalize()+(' et al.' if ',' in c[1] else '')+f" ({c[2]})"}
src={};out=[];seen=set();n=0;cnt=collections.Counter();st=collections.Counter()
for r in R['S8']:
    c=keep.get(key(r[12]))
    if not c: continue
    sid='PM_'+re.sub(r'\W','',c[0])[-12:];src[sid]=cite(c)
    ch=str(r[1]).strip()[:2];iv=re.sub(r'\s+',' ',str(r[4])).strip()
    if ch not in CEN: st['bad']+=1;continue
    mks=[x.strip() for x in re.split(r'\s*[–\-—/]\s*(?=[A-Za-z])|\s+and\s+|,\s*',iv) if x.strip()]
    P=[(m,pos(m,ch)) for m in mks];P=[(m,p) for m,p in P if p]
    if not P: st['nopos']+=1;continue
    if len(P)>1 and max(p[0] for m,p in P)-min(p[0] for m,p in P)>60: st['far']+=1;continue
    v1=[(m,p) for m,p in P if p[1]=='RefSeq v1.0'];L=v1 or P;rf=L[0][1][1];ps=[p[0] for m,p in L]
    kk=(sid,ch,iv)
    if kk in seen: continue
    seen.add(kk);n+=1;cnt[sid]+=1;two=len(ps)>1
    out.append(dict(id='PM-%04d'%n,name=f"QPm {src[sid]['short'].split(' ')[0]} {ch} ({iv})",cls='QTL',disease='Powdery mildew',chr=ch,start=round(min(ps),4) if two else None,end=round(max(ps),4) if two else None,peak=None if two else round(ps[0],4),
        ref=rf,marker=', '.join(m for m,p in L),pve=str(r[9] or ''),stat=f"LOD {r[8]}" if r[8] else '',source=sid,grade='C',curation='auto (supplementary compilation, in-silico)',
        evidence=f"Table S8 of the powdery mildew meta-QTL Data Sheet 1: {r[0]}, {r[2]}, interval {iv}, cross {r[3]}",
        note=f"QTL ({r[2]}, {r[3]}, peak {r[5]} cM) from the original paper, as compiled in Table S8 of the powdery mildew meta-QTL Data Sheet 1. Placed in silico at the published physical position of the interval marker(s) {', '.join(f'{m} {p[0]:.2f} Mb ({p[2]})' for m,p in L)} ({rf}). Not yet checked by hand against the original paper."))
for r in R['S6']:
    c=keep.get(key(r[0]));mk=str(r[1]).strip();ch=str(r[2]).strip()[:2]
    if not c or ch not in CEN: continue
    sid='PM_'+re.sub(r'\W','',c[0])[-12:];src[sid]=cite(c);kk=(sid,mk,ch)
    if kk in seen: continue
    p=pos(mk,ch)
    if p and p[3]=='old': v,rf,how=p[0],p[1],'published physical position of the marker'
    else:
        m=re.match(r'\s*(\d+)',str(r[3]))
        if not m: continue
        v,rf,how=int(m.group(1))/1e6,'RefSeq v1.0','bp position given in Table S6'
    seen.add(kk);n+=1;cnt[sid]+=1
    out.append(dict(id='PM-%04d'%n,name=f"{mk} (powdery mildew)",cls='MTA',disease='Powdery mildew',chr=ch,start=None,end=None,peak=round(v,4),ref=rf,marker=mk,pve='',stat='',source=sid,grade='C',
        curation='auto (supplementary compilation)',evidence=f"Table S6 of the powdery mildew meta-QTL Data Sheet 1: {mk} {ch}",
        note=f"Significant GWAS marker from the original paper, as compiled in Table S6 of the powdery mildew meta-QTL Data Sheet 1. Position: {how} ({rf}). Not yet checked by hand against the original paper."))
json.dump(dict(reports=out,sources={k:v for k,v in src.items() if cnt[k]},counts={src[k]['short']:v for k,v in cnt.items()}),open(sys.argv[3],'w'),indent=0)
print(n,st,collections.Counter(r['cls'] for r in out));print({src[k]['short']:v for k,v in cnt.items()})
