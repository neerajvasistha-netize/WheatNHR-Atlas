"""FHB QTL from Table S1 (mmc2) of the Crop Journal FHB meta-QTL study (S221451412030177X), credited to the original papers
(reference list in Table S4, mmc5). The table already gives physical positions (RefSeq v1.0) of the QTL interval markers
and of the peak marker. Meta-QTLs are not added. Papers already in the atlas (same DOI) are skipped.
usage: python3 build79.py fhb.json fhbcr.json atlas_dois.txt out.json"""
import json,sys,re,collections
J=json.load(open(sys.argv[1]));CR=json.load(open(sys.argv[2]));AT=set(open(sys.argv[3]).read().split())
def words(s): return set(w for w in re.findall(r'[a-z]{4,}',s.lower()))
ok={}
for k,c in CR.items():
    if not c: continue
    t=J['refs'][k];ov=len(words(c[3])&words(t))/max(1,len(words(c[3])))
    if ov>=0.6 and c[0] not in AT: ok[int(k)]=c
src={};out=[];n=0;cnt=collections.Counter();seen=set()
for q in J['Q']:
    c=ok.get(q['rn'])
    if not c: continue
    sid='CJ_'+re.sub(r'\W','',c[0])[-12:]
    src[sid]={'citation':f"{c[1]} ({c[2]}) {c[3].rstrip('.')}. {c[4]} {c[5]}:{c[6]}".replace(':None','').replace(' None',''),'url':'https://doi.org/'+c[0],
              'short':c[1].split(',')[0].split(' ')[0]+(' et al.' if ',' in c[1] else '')+f" ({c[2]})"}
    try: l,r=float(q['l'])/1e6,float(q['r'])/1e6
    except: continue
    s,e=min(l,r),max(l,r)
    try: pk=float(q['pkp'])/1e6
    except: pk=None
    k=(sid,q['chr'],q['ci'])
    if k in seen: continue
    seen.add(k);n+=1;cnt[src[sid]['short']]+=1
    ci=re.sub(r'\s+',' ',str(q['ci'])).strip()
    out.append(dict(id='CJ-%03d'%n,name=f"QFhb {src[sid]['short'].split(' ')[0]} {q['chr']} ({ci})",cls='QTL',disease='Fusarium head blight',chr=q['chr'],start=round(s,4),end=round(e,4),peak=round(pk,4) if pk and s<=pk<=e else None,
        ref='RefSeq v1.0',marker=ci.replace(' - ',', '),pve=str(q['pve'] or ''),stat=f"LOD {q['lod']}" if q['lod'] not in('',None) else '',source=sid,grade='C',curation='auto (supplementary compilation)',
        evidence=f"Table S1 of the Crop Journal FHB meta-QTL study: {q['src']}, type {q['type']} resistance, interval {ci}, peak marker {q['pk']}",
        note=f"Type {str(q['type']).strip()} FHB resistance QTL from {str(q['src']).strip()} ({q['pop']}), as compiled in Table S1 of the Crop Journal FHB meta-QTL study (physical positions of the interval and peak markers on RefSeq v1.0). Not yet checked by hand against the original paper."))
json.dump(dict(reports=out,sources=src,counts=cnt),open(sys.argv[4],'w'),indent=0)
print(n,len(src),len(ok))
