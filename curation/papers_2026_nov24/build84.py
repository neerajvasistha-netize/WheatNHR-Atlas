import json,re,unicodedata,collections,sys
S='/tmp/claude-0/-home-claude/00583e43-a309-532b-9f1b-7957bcf4ea95/scratchpad'
exec(open('records84.py').read())
sel=json.load(open(S+'/y26/sel.json')); meta=json.load(open(S+'/y26/meta.json'))
allr=json.load(open(S+'/ms/all.json'))
used=collections.Counter(re.sub(r'[a-z]?\)$','',x.get('short') or '') for x in allr)
taken=set(x.get('short') for x in allr)
def ini(n):
    p=n.replace('.',' ').split(); 
    if not p: return ''
    sur=p[-1]; ii=''.join(x[0] for x in p[:-1] if x[0].isalpha()).upper()
    return f'{sur} {ii}'
YR={v['doi']:v['year'] for v in sel}
LAT=lambda L:[x for x in L if all(ord(c)<0x2E80 for c in x)]
def cit(doi):
    m=meta[doi]; au=LAT(m['authors']); a=', '.join(ini(x) for x in au[:6])+(' et al.' if len(au)>6 else '')
    b=m.get('biblio') or {}; vp=''
    if b.get('volume'): vp=f" {b['volume']}"+(f":{b['first_page']}"+(f"-{b['last_page']}" if b.get('last_page') and b.get('last_page')!=b.get('first_page') else '') if b.get('first_page') else '')
    t=re.sub(r'\s+',' ',re.sub('<[^>]+>','',m['title'] or '')).strip()
    return f"{a} ({YR[doi]}) {t}. {m['journal']}{vp}"
def short(doi):
    au=LAT(meta[doi]['authors']); sur=au[0].split()[-1]; y=YR[doi]
    base=f"{sur} et al. ({y})" if len(au)>2 else (f"{sur} and {au[1].split()[-1]} ({y})" if len(au)==2 else f"{sur} ({y})")
    if base not in taken: taken.add(base); return base
    for L in 'bcdefghij':
        s=base[:-1]+L+')'
        if s not in taken: taken.add(s); return s
src={};papers={};reports=[];n=0
for r in R:
    i,name,cls,dis,ch,s,e,p,ref,mk,pve,stat,g,ev,note=r
    doi=sel[i]['doi']; key='Y26_'+re.sub(r'[^a-z0-9]','',doi.split('/',1)[1])[:24]
    if key not in src: src[key]={'citation':cit(doi),'url':'https://doi.org/'+doi,'short':short(doi)}
    n+=1
    rec={'id':f'Y6-{n:03d}','name':name,'cls':cls,'disease':dis,'chr':ch,'start':s,'end':e,'peak':p,'ref':ref or ('RefSeq v1.0' if g!='D' else ''),'marker':mk,'pve':pve,'stat':stat,'source':key,'grade':g,
         'curation':'auto (2026 paper abstract/text, hand-read)','evidence':ev,'note':(note+' Read from the abstract or open full text on 24 Nov 2026; not yet checked against the full paper.').strip()}
    if g=='C': rec['curation']='auto (2026 paper text, candidate-gene position)'
    reports.append(rec)
# tracker for all screened papers
added=collections.Counter(r[0] for r in R)
for i,v in enumerate(sel):
    doi=v['doi'];m=meta.get(doi,{})
    k='Y26P_'+re.sub(r'[^a-z0-9]','',doi.split('/',1)[1])[:24]
    au=LAT(m.get('authors') or [v.get('first','')])
    papers[k]={'id':k,'authors':', '.join(ini(x) for x in au[:6])+(' et al.' if len(au)>6 else '')+f" ({v['year']})",'title':re.sub(r'\s+',' ',re.sub('<[^>]+>','',v['title'])).strip(),'journal':m.get('journal') or v['journal'],
      'check':'A','assignee':'','type':'','disease':'','in_atlas':bool(added[i]),'status':'Added to atlas' if added[i] else 'To search','link':'https://doi.org/'+doi,
      'where':'abstract / open full text','notes':(f"{added[i]} loci added (seed_ad84, Crossref 2026 search)" if added[i] else 'Found by Crossref 2026 search; no locus position in the abstract or open text (alien chromosome, insect pest, known gene, or positions only in tables/supplements)')}
json.dump({'reports':reports,'sources':src,'papers':papers},open('seed84_pre.json','w'),ensure_ascii=False,indent=0)
print(len(reports),len(src),len(papers),collections.Counter(r['grade'] for r in reports))
for k,v in src.items(): print(v['short'],'|',v['citation'][:130])
