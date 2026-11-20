import json,re,sys
W='/home/claude/wheatnhr-atlas/web/seed_ad73.json'
d=json.load(open(W));meta=json.load(open('rdfmeta.json'));P={x['h']:x for x in json.load(open('s1.json'))}
fr=json.load(open('../fig/fig_res.json'))
n=max(int(r['id'][2:]) for r in d['reports'])
GD={'yr':'Stripe rust','lr':'Leaf rust','sr':'Stem rust','pm':'Powdery mildew','ml':'Powdery mildew','cb':'Common bunt','bt':'Common bunt'}
def dis(name):
    m=re.match(r'Q?([A-Za-z]{2})',name);return GD.get(m.group(1).lower())
out=[]
for r in fr:
    pp=P[r['h']];m=meta[pp['path']];key='PX_'+r['h']
    if key not in d['sources']:
        au=m['au'];yr=m['yr'] or pp['yr']
        names=', '.join(f"{s} {''.join(w[0] for w in re.split(r'[ .-]+',g) if w)}" for s,g in au[:6])+(' et al.' if len(au)>6 else '')
        short=(au[0][0])+(' et al.' if len(au)>2 else (' and '+au[1][0] if len(au)==2 else ''))+f' ({yr})'
        d['sources'][key]={'citation':f"{names} ({yr}) {m['ti']}. {m['j']}{' '+m['vol'] if m['vol'] else ''}{':'+m['pg'] if m['pg'] else ''}".strip(),'url':'https://doi.org/'+m['doi'] if m['doi'] else '','short':short}
    n+=1;mk=', '.join(x.rstrip('*') for x in r['markers'][:6]);fig=r['fig'].replace('\n',' ')
    out.append({'id':'PX%04d'%n,'name':r['name'],'cls':'QTL' if r['name'].startswith('Q') else 'DG','disease':dis(r['name']),'chr':r['chr'],
      'start':r['start'],'end':r['end'],'peak':None,'ref':r['ref'],'marker':mk,'pve':'','stat':'','source':key,'grade':'C','curation':'auto (paper figure)',
      'evidence':f"{fig}: markers read from the figure ({mk})",
      'note':f"Read automatically from {fig} of the paper (Nov 2026): the markers drawn next to this locus were read from the image and placed in silico at their published physical positions ({r['ref']}). Not yet checked by hand against the paper."})
    p=d['papers'].get(key)
    if p:
        p['status']='Added to atlas';p['in_atlas']=True
        p['notes']=(p['notes']+'; ' if p['notes'] else '')+f"{r['name']} placed from {fig} (marker positions in silico)"
        p['where']=(p['where'] or 'Figures')
json.dump([{'rec':x} for x in out],open('../pan/need66.json','w'));json.dump(out,open('figrecs.json','w'))
json.dump(d,open(W,'w'),ensure_ascii=False,separators=(',',':'))
print(len(out),[x['disease'] for x in out])
