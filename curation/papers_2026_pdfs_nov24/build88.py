import json,re,collections
S='/tmp/claude-0/-home-claude/00583e43-a309-532b-9f1b-7957bcf4ea95/scratchpad'
exec(open('records88.py').read())
sel=json.load(open(S+'/y26/sel.json')); meta=json.load(open(S+'/y26/meta.json'))
YR={v['doi']:v['year'] for v in sel}
LAT=lambda L:[x for x in L if all(ord(c)<0x2E80 for c in x)]
def ini(n):
    p=n.replace('.',' ').split(); return f"{p[-1]} {''.join(x[0] for x in p[:-1] if x[0].isalpha()).upper()}"
def cit(doi):
    m=meta[doi]; au=LAT(m['authors']); a=', '.join(ini(x) for x in au[:6])+(' et al.' if len(au)>6 else '')
    b=m.get('biblio') or {}; vp=''
    if b.get('volume'): vp=f" {b['volume']}"+(f":{b['first_page']}"+(f"-{b['last_page']}" if b.get('last_page') and b['last_page']!=b['first_page'] else '') if b.get('first_page') else '')
    return f"{a} ({YR[doi]}) {re.sub(r'\s+',' ',re.sub('<[^>]+>','',m['title'])).strip()}. {m['journal']}{vp}"
# existing short labels in the atlas and in seeds 84-87
taken=set(x.get('short') for x in json.load(open(S+'/ms/all.json')))
for f in ('seed_ad84.json','seed_ad85.json','seed_ad86.json','seed_ad87.json'):
    taken|={v['short'] for v in json.load(open('/home/claude/wheatnhr-atlas/web/'+f))['sources'].values()}
def short(doi):
    au=LAT(meta[doi]['authors']); base=f"{au[0].split()[-1]} et al. ({YR[doi]})"
    for s in [base]+[base[:-1]+L+')' for L in 'bcdefghij']:
        if s not in taken: taken.add(s); return s
# in-silico placements (spot blotch, Exserohilum)
ins=json.load(open(S+'/p26z/insilico.json'))
DSB='10.1016/j.pmpp.2026.103454'; DES='10.3390/plants15152351'
for x in ins:
    if not x['pos']: continue
    ch,p,ref,srcp=x['pos'][:4]
    tr='/'.join(x['traits']); env='; '.join(x['env'])
    if x['dz']=='SB':
        R.append((DSB,f"MTA SB {x['marker']}",'MTA','Spot blotch',ch,None,None,p,ref,x['marker'],'','FarmCPU','C',
          f"Table {'4' if 'AUDPC' in x['traits'] else '5'}: {x['marker']}, {ch}, {x['cm']} cM; trait {tr}; isolate/year/site {env}",
          f"Spot blotch GWAS ({tr}) across Bipolaris sorokiniana isolates at BHU and BISA, 2019-2020. The paper gives only cM; placed in silico at the published physical position of the marker ({srcp})."))
    else:
        R.append((DES,f"MTA ESB {x['marker']}",'MTA','Exserohilum spike blight',ch,None,None,p,ref,x['marker'],'','','C',
          f"Table 3: {x['marker']}, {ch}, {x['cm']} cM; trait {tr}; {env}",
          f"Exserohilum spike blight (E. rostratum) GWAS in the WAMI panel ({tr}). The paper gives cM and lists candidate genes at Mb equal to the cM value, so the candidate genes were not used; placed in silico at the published physical position of the marker ({srcp})."))
src={};reports=[];n=0;cnt=collections.Counter()
for doi,name,cls,dis,ch,s,e,p,ref,mk,pve,stat,g,ev,note in R:
    key='Y26_'+re.sub(r'[^a-z0-9]','',doi.split('/',1)[1])[:24]
    if key not in src: src[key]={'citation':cit(doi),'url':'https://doi.org/'+doi,'short':short(doi)}
    n+=1;cnt[doi]+=1
    reports.append({'id':f'YP6-{n:03d}','name':name,'cls':cls,'disease':dis,'chr':ch,'start':s,'end':e,'peak':p,'ref':ref,'marker':mk,'pve':pve,'stat':stat,'source':key,'grade':g,
      'curation':'auto (paper full text, in-silico)' if g=='C' else 'auto (paper full text)','evidence':ev,'note':note+' Read from the full-text PDF supplied by the owner on 24 Nov 2026; not yet checked by hand.'})
s84=json.load(open('/home/claude/wheatnhr-atlas/web/seed_ad84.json'))
papers={}
for k,pp in s84['papers'].items():
    d=pp['link'].replace('https://doi.org/','').lower()
    if d in cnt:
        q=dict(pp);q.update(status='Added to atlas',in_atlas=True,where='full text (PDF supplied by owner)',notes=f'{cnt[d]} records added from the full text (seed_ad88)');papers[k]=q
json.dump({'reports':reports,'sources':src,'papers':papers},open('seed88_pre.json','w'),ensure_ascii=False)
print(len(reports),len(src),len(papers));[print(v['short'],cnt[v['url'][16:]]) for v in src.values()]
