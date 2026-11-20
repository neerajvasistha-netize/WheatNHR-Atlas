import json,re
W='/home/claude/wheatnhr-atlas/web/seed_ad73.json'
d=json.load(open(W));meta=json.load(open('rdfmeta.json'));P={x['h']:x for x in json.load(open('s1.json'))}
sr=json.load(open('../supp/supp_res.json'))
have={(r['name'],r['chr']) for r in d['reports']}
n=max(int(r['id'][2:]) for r in d['reports'])
GD={'yr':'Stripe rust','lr':'Leaf rust','sr':'Stem rust','pm':'Powdery mildew','ml':'Powdery mildew','cb':'Common bunt','bt':'Common bunt'}
def dis(name):
    m=re.match(r'Q?([A-Za-z]{2})',name);return GD.get(m.group(1).lower())
out=[]
for r in sr:
    if (r['name'],r['chr']) in have or not dis(r['name']): continue
    have.add((r['name'],r['chr']))
    pp=P[r['h']];m=meta[pp['path']];key='PX_'+r['h']
    if key not in d['sources']:
        au=m['au'];yr=m['yr'] or pp['yr']
        names=', '.join(f"{s} {''.join(w[0] for w in re.split(r'[ .-]+',g) if w)}" for s,g in au[:6])+(' et al.' if len(au)>6 else '')
        short=(au[0][0])+(' et al.' if len(au)>2 else (' and '+au[1][0] if len(au)==2 else ''))+f' ({yr})'
        d['sources'][key]={'citation':f"{names} ({yr}) {m['ti']}. {m['j']}{' '+m['vol'] if m['vol'] else ''}{':'+m['pg'] if m['pg'] else ''}".strip(),'url':'https://doi.org/'+m['doi'] if m['doi'] else '','short':short}
    gw=re.search(r'genome-wide assoc|association|GWAS',m['ti'],re.I)
    n+=1;mk=', '.join(r['markers'][:6])
    insil=r['how'].startswith('in-silico')
    out.append({'id':'PX%04d'%n,'name':r['name'],'cls':('MTA' if gw else 'QTL') if r['name'].startswith('Q') else 'DG','disease':dis(r['name']),'chr':r['chr'],
      'start':None if r['start'] is None else round(r['start'],4),'end':None if r['end'] is None else round(r['end'],4),'peak':None if r['peak'] is None else round(r['peak'],4),
      'ref':r['ref'],'marker':mk,'pve':'','stat':'','source':key,'grade':'C','curation':'auto (supplementary table)',
      'evidence':f"Supplementary file {r['file']}: {r['row']}",
      'note':("Read automatically from the paper's supplementary table (Nov 2026): "+("the markers in the locus row were placed in silico at their published physical positions" if insil else "physical position given in the locus row")+f" ({r['ref']}). Not yet checked by hand against the paper.")})
    p=d['papers'].get(key)
    if p:
        p['status']='Added to atlas';p['in_atlas']=True
        p['where']='Supplementary table' if not p['where'] else p['where']+'; supplementary table'
out_by={}
for x in out: out_by.setdefault(x['source'],[]).append(x['name'])
for key,ns in out_by.items():
    p=d['papers'].get(key)
    if p: p['notes']=(p['notes']+'; ' if p['notes'] else '')+f"{len(ns)} loci placed from the supplementary tables (downloaded 20 Nov 2026)"
json.dump([{'rec':x} for x in out],open('../pan/need66.json','w'));json.dump(out,open('supprecs.json','w'))
json.dump(d,open(W,'w'),ensure_ascii=False,separators=(',',':'))
import collections;print(len(out),collections.Counter(x['disease'] for x in out),collections.Counter(x['source'] for x in out))
