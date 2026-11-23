"""Stripe rust QTL (Table S1, marker intervals placed in silico) and GWAS MTAs (Table S5) from the supplementary file of
Kumar et al. 2023 (BMC Genomics stripe rust meta-QTL; 12864_2023_9336_MOESM1_ESM.xlsx), credited to the original papers.
Meta-QTLs are not added. Only studies whose paper was identified on Crossref (first author, year, stripe/yellow rust title)
and that are not already in the atlas (DOI check) are used.
QTL: interval markers looked up for published physical positions (or e-PCR); kept if >=1 marker placed and placed markers lie
within 60 Mb. MTA: published marker position if found, else the bp value of Table S5.
usage: python3 build78.py s1place.json s5rows.json keep.json out.json"""
import json,sys,re,collections
s1=json.load(open(sys.argv[1]));s5=json.load(open(sys.argv[2]));keep=json.load(open(sys.argv[3]))
def cite(c):
    t=c[3].rstrip('.');return {'citation':f"{c[1]} ({c[2]}) {t}. {c[4]} {c[5]}:{c[6]}".replace(' None:None','').replace(':None',''),'url':'https://doi.org/'+c[0],
        'short':c[1].split(',')[0].split(' ')[0].capitalize()+(' et al.' if ',' in c[1] else '')+f" ({c[2]})"}
src={};out=[];seen=set();cnt=collections.Counter();n=0
for x in s1:
    c=keep[x['key']];sid='K8_'+re.sub(r'\W','',c[0])[-12:];src[sid]=cite(c)
    k=(sid,x['chr'],x['interval'])
    if k in seen: continue
    seen.add(k)
    ok=x['ok'];v1=[(m,p) for m,p in ok if p[1]=='RefSeq v1.0'];L=v1 or ok;rf=L[0][1][1]
    ps=[p[0] for m,p in L];n+=1;cnt[sid]+=1
    two=len(ps)>1
    out.append(dict(id='K8-%03d'%n,name=f"QYr {x['qid'].split('_')[0]} {x['chr']} ({x['interval']})",cls='QTL',disease='Stripe rust',chr=x['chr'],
        start=round(min(ps),4) if two else None,end=round(max(ps),4) if two else None,peak=None if two else round(ps[0],4),ref=rf,
        marker=', '.join(m for m,p in L),pve=str(x['r2'] or ''),stat=f"LOD {x['lod']}" if x['lod'] else '',source=sid,grade='C',curation='auto (supplementary compilation, in-silico)',
        evidence=f"Table S1 of Kumar et al. 2023: {x['qid']}, interval {x['interval']}",
        note=f"QTL {x['qid']} (trait {x['trait']}, peak {x['cM']} cM) from the original paper, as compiled in Table S1 of Kumar et al. (2023, stripe rust meta-QTL). Placed in silico at the published physical position of the interval marker(s) {', '.join(m for m,p in L)} ({rf}). Not yet checked by hand against the original paper."))
for x in s5:
    c=keep['S5:'+x['ref']];sid='K8_'+re.sub(r'\W','',c[0])[-12:];src[sid]=cite(c)
    k=(sid,x['mk'],x['chr'])
    if k in seen or x['mk'] in('NA','-',''): continue
    if x['v1'] is not None: p,rf,how=x['v1'],'RefSeq v1.0','published physical position of the marker'
    elif x['v2'] is not None: p,rf,how=x['v2'],'RefSeq v2.1','published physical position of the marker'
    elif x['tab']: p,rf,how=x['tab'],'RefSeq v1.0','bp position given in Table S5'
    else: continue
    seen.add(k);n+=1;cnt[sid]+=1
    out.append(dict(id='K8-%03d'%n,name=f"{x['mk']} (stripe rust)",cls='MTA',disease='Stripe rust',chr=x['chr'],start=None,end=None,peak=round(p,4),ref=rf,marker=x['mk'],pve='',stat='',source=sid,grade='C',
        curation='auto (supplementary compilation)',evidence=f"Table S5 of Kumar et al. 2023: {x['mk']} {x['chr']}",
        note=f"Significant GWAS marker from the original paper, as compiled in Table S5 of Kumar et al. (2023, stripe rust meta-QTL). Position: {how} ({rf}). Not yet checked by hand against the original paper."))
json.dump(dict(reports=out,sources=src,counts={src[k]['short']:v for k,v in cnt.items()}),open(sys.argv[4],'w'),indent=0)
print(n,{src[k]['short']:v for k,v in cnt.items()})
