"""Records for seed_ad75.json: papers downloaded into Zotero (WheatRLA_biotroph_not_downloaded) and Downloads/Papers, read Nov 2026."""
import json,re,sys,collections
sys.path.insert(0,'../exp')
from load2 import load
R,S,_=load()
def nk(s): return re.sub(r'[^a-z0-9]','',s.lower().split(' (')[0])
EX=collections.defaultdict(set)
for r in R.values():
    for part in re.split(r'[/,;]| - ',r['name']): EX[nk(part)].add(r.get('chr'))
    for g in re.findall(r'\b(?:Yr|Lr|Sr|Pm)\d+[a-z]?\b',r['name']): EX[nk(g)].add(r.get('chr'))
for f in ['seed_ad73','seed_ad74']:
    for r in json.load(open('/home/claude/wheatnhr-atlas/web/'+f+'.json'))['reports']: EX[nk(r['name'])].add(r['chr'])
meta=json.load(open('meta5.json'));P={x['h']:x for x in json.load(open('s1.json'))}
def names(au): return ', '.join(f"{s} {''.join(w[0] for w in re.split(r'[ .-]+',g) if w)}" for s,g in au[:6])+(' et al.' if len(au)>6 else '')
def short(au,yr): return (au[0][0] if au else 'Anon.')+(' et al.' if len(au)>2 else (' and '+au[1][0] if len(au)==2 else ''))+f' ({yr})'
out=[];src={};seen=set();n=0
rows_by=json.load(open(sys.argv[2])) if len(sys.argv)>2 else None
for x in json.load(open('s2.json')):
    m=meta[x['h']];rows=[]
    for r in x['rows']:
        name=r['name'].replace('‐','-').replace('–','-');k=(nk(name),r['chr'])
        if r['chr'] in EX.get(nk(name),set()) or k in seen: continue
        seen.add(k);rows.append((name,r))
    if not rows: continue
    key='PY_'+x['h'];yr=m['yr'] or P[x['h']]['yr']
    src[key]={'citation':f"{names(m['au'])} ({yr}) {m['ti']}. {m['j']}{' '+m['vol'] if m['vol'] else ''}{':'+m['pg'] if m['pg'] else ''}".strip(),'url':'https://doi.org/'+m['doi'] if m['doi'] else '','short':short(m['au'],yr)}
    gw=re.search(r'genome-wide assoc|association|GWAS',m['ti'],re.I)
    for name,r in rows:
        n+=1;meth=r['how']
        out.append({'id':'PY%04d'%n,'name':name,'cls':('MTA' if gw and r['cls']=='QTL' else r['cls']),'disease':r['disease'],'chr':r['chr'],
            'start':None if r['start'] is None else round(r['start'],4),'end':None if r['end'] is None else round(r['end'],4),'peak':None if r['peak'] is None else round(r['peak'],4),
            'ref':r['ref'],'marker':r['marker'],'pve':'','stat':'','source':key,'grade':'C','curation':'auto (paper text)','evidence':r['ev'] or meth,
            'note':f"Read automatically from the paper's full text (Nov 2026): {meth}. "+("Position as given in the paper" if not meth.startswith('in-silico') else "Position from the published physical positions of the linked markers named with this locus in the paper")+f" ({r['ref']}). Not yet checked by hand against the paper."})
# paper tracker
LOG=collections.defaultdict(list)
for h,nm,why in json.load(open('s2log.json')): LOG[h].append(nm.replace('‐','-'))
S2={x['h']:x for x in json.load(open('s2.json'))};used=collections.Counter(r['source'] for r in out);papers={}
for h,pp in P.items():
    m=meta[h];yr=m['yr'] or pp['yr'];ti=m['ti'] or pp['ti'];tl=ti.lower()
    typ='GWAS' if re.search(r'genome-wide assoc|association|gwas',tl) else 'Fine mapping' if 'fine map' in tl else 'QTL' if re.search(r'\bqtl|quantitative trait',tl) else 'Cytogenetics' if re.search(r'cytogenet|translocation|addition|substitution|introgression|telocentric|monosomic',tl) else 'Gene' if re.search(r'gene|\b(yr|lr|sr|pm)\d',tl) else 'Other'
    placed=[r['name'] for r in S2[h]['rows']];added=used.get('PY_'+h,0);unpl=sorted(set(LOG[h]))
    notes=[]
    if added: notes.append(f'{added} locus records added automatically (seed_ad75, grade C; check against the paper)')
    if placed and len(placed)>added: notes.append(f'{len(placed)-added} placed loci were already in the atlas under the same name')
    if unpl: notes.append('named but not placed (no position found in the text): '+', '.join(unpl[:12])+(' ...' if len(unpl)>12 else ''))
    if not placed and not unpl: notes.append('no named wheat locus found by the automatic reading; read by hand')
    papers['PY_'+h]={'id':'PY_'+h,'disease':pp['dis'][0][0] if pp['dis'][0][1]>0 else '','type':typ,'authors':names(m['au'])+f' ({yr})','title':ti,
      'journal':(m['j']+(' '+m['vol'] if m['vol'] else '')+(':'+m['pg'] if m['pg'] else '')).strip(),'link':'https://doi.org/'+m['doi'] if m['doi'] else '',
      'check':'A','in_atlas':bool(placed),'status':'Added to atlas' if placed else 'To search','assignee':'','where':'Main text and tables (automatic reading)' if placed else '','notes':'; '.join(notes)}
json.dump([{'rec':r} for r in out],open(sys.argv[1],'w'));json.dump({'sources':src,'papers':papers},open('srcpap75.json','w'),ensure_ascii=False)
print(len(out),'records',len(src),'sources',len(papers),'papers',collections.Counter(r['disease'] for r in out))
