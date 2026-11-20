"""Records for seed_ad73.json from the automatic reading of the Zotero 'Exported' papers (Nov 2026)."""
import json,re,sys,collections
from load2 import load
R,S,_=load()
def nk(s): return re.sub(r'[^a-z0-9]','',s.lower().split(' (')[0])
EX=collections.defaultdict(set)
for r in R.values():
    for part in re.split(r'[/,;]| - ',r['name']): EX[nk(part)].add(r.get('chr'))
    for g in re.findall(r'\b(?:Yr|Lr|Sr|Pm)\d+[a-z]?\b',r['name']): EX[nk(g)].add(r.get('chr'))
meta=json.load(open('rdfmeta.json')); P={x['h']:x for x in json.load(open('s1.json'))}
out=[];src={};seen=set();n=0
for x in json.load(open('s2.json')):
    pp=P[x['h']];m=meta[pp['path']]
    rows=[]
    for r in x['rows']:
        name=r['name'].replace('‐','-').replace('–','-')
        k=(nk(name),r['chr'])
        if r['chr'] in EX.get(nk(name),set()) or k in seen: continue
        seen.add(k);rows.append((name,r))
    if not rows: continue
    key='PX_'+x['h']
    au=m['au'];yr=m['yr'] or pp['yr']
    names=', '.join(f"{s} {''.join(w[0] for w in re.split(r'[ .-]+',g) if w)}" for s,g in au[:6])+(' et al.' if len(au)>6 else '')
    short=(au[0][0] if au else pp['au'].split(' et al')[0])+(' et al.' if len(au)>2 else (' and '+au[1][0] if len(au)==2 else ''))+f' ({yr})'
    cit=f"{names} ({yr}) {m['ti'] or pp['ti']}. {m['j']}{' '+m['vol'] if m['vol'] else ''}{':'+m['pg'] if m['pg'] else ''}".strip()
    src[key]={'citation':cit,'url':'https://doi.org/'+m['doi'] if m['doi'] else '','short':short}
    for name,r in rows:
        n+=1
        meth=r['how']
        note=(f"Read automatically from the paper's full text (Nov 2026): {meth}. "
              +("Position as given in the paper" if not meth.startswith('in-silico') else "Position from the published physical positions of the linked markers named with this locus in the paper")
              +f" ({r['ref']}). Not yet checked by hand against the paper's tables and figures.")
        out.append({'id':'PX%04d'%n,'name':name,'cls':('MTA' if r['cls']=='QTL' and re.search(r'genome-wide assoc|association|GWAS',m['ti'] or pp['ti'],re.I) else r['cls']),'disease':r['disease'],'chr':r['chr'],
            'start':None if r['start'] is None else round(r['start'],4),'end':None if r['end'] is None else round(r['end'],4),
            'peak':None if r['peak'] is None else round(r['peak'],4),'ref':r['ref'],'marker':r['marker'],'pve':'','stat':'',
            'source':key,'grade':'C','curation':'auto (paper text)','evidence':r['ev'] or meth,'note':note})
json.dump([{'rec':r} for r in out],open(sys.argv[1],'w'));json.dump(src,open('src73.json','w'),ensure_ascii=False)
print(len(out),'records',len(src),'sources',collections.Counter(r['disease'] for r in out))
# ---- paper tracker entries for every new paper in the Exported folder
LOG=collections.defaultdict(list)
for h,name,why in json.load(open('s2log.json')): LOG[h].append((name.replace('‐','-'),why))
S2={x['h']:x for x in json.load(open('s2.json'))}
used=collections.Counter(r['source'] for r in out)
papers={};trk=[]
for h,pp in P.items():
    m=meta[pp['path']];au=m['au'];yr=m['yr'] or pp['yr']
    ti=m['ti'] or pp['ti'];tl=ti.lower()
    typ='GWAS' if re.search(r'genome-wide assoc|association|gwas',tl) else 'Fine mapping' if 'fine map' in tl else 'QTL' if re.search(r'\bqtl|quantitative trait',tl) else 'Cytogenetics' if re.search(r'cytogenet|translocation|addition|substitution|introgression',tl) else 'Gene' if re.search(r'gene|\b(yr|lr|sr|pm)\d',tl) else 'Other'
    placed=[r['name'] for r in S2[h]['rows']]
    added=used.get('PX_'+h,0);unpl=[n for n,w in LOG[h]]
    if added: st='Added to atlas'
    elif placed: st='Added to atlas'
    else: st='To search'
    notes=[]
    if added: notes.append(f'{added} locus records added automatically (seed_ad73, grade C; check against the paper)')
    if placed and len(placed)>added: notes.append(f'{len(placed)-added} placed loci were already in the atlas under the same name')
    if unpl: notes.append('named but not placed (no position found in the text): '+', '.join(sorted(set(unpl))[:12])+(' ...' if len(set(unpl))>12 else ''))
    if not placed and not unpl: notes.append('no named wheat locus found by the automatic reading (alien chromosome lines, marker screening or figure-only maps); read by hand')
    dis=pp['dis'][0][0] if pp['dis'][0][1]>0 else ''
    papers['PX_'+h]={'id':'PX_'+h,'disease':dis,'type':typ,'authors':(', '.join(f"{s} {''.join(w[0] for w in re.split(r'[ .-]+',g) if w)}" for s,g in au[:5])+(' et al.' if len(au)>5 else ''))+f' ({yr})',
        'title':ti,'journal':(m['j']+(' '+m['vol'] if m['vol'] else '')+(':'+m['pg'] if m['pg'] else '')).strip(),'link':'https://doi.org/'+m['doi'] if m['doi'] else '',
        'check':'A','in_atlas':bool(placed),'status':st,'assignee':'','where':'Main text and tables (automatic reading)' if placed else '','notes':'; '.join(notes)}
    trk.append([pp['path'],yr,ti,m['doi'],typ,dis,st,added,len(placed),len(set(unpl)),'; '.join(notes)])
json.dump(papers,open('papers73.json','w'),ensure_ascii=False)
import csv
with open('tracker73.tsv','w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['pdf','year','title','doi','type','disease','status','records_added','loci_placed','loci_unplaced','notes']);w.writerows(trk)
print('papers',len(papers),collections.Counter(p['status'] for p in papers.values()))
