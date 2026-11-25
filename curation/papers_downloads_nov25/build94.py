"""Build web/seed_ad94.json from records94.py (old full-text papers in Downloads/Paper_2026/Papers).
Markers are placed on RefSeq v1.0 with the atlas marker table (data/markerpos_all.json.gz) and the consensus marker table of
curation/fhb_metaqtl_cropj2021_esm; two or more placed markers within 30 Mb give an interval, otherwise the first placed marker."""
import json,re,sys,collections
S='/tmp/claude-0/-home-claude/00583e43-a309-532b-9f1b-7957bcf4ea95/scratchpad/p26'
sys.path.insert(0,S); from look2 import pos
exec(open('records94.py').read())
M=json.load(open(S+'/papers_meta.json')); OA=json.load(open(S+'/oa_meta.json'))
A=json.load(open(S+'/atlas_all.json'))
taken={v.get('short') for v in A['sources'].values()}
LAT=lambda L:[x for x in L if all(ord(c)<0x2E80 for c in x)]
def ini(n):
    p=n.replace('.',' ').replace('-',' ').split(); return f"{p[-1]} {''.join(x[0] for x in p[:-1] if x[0].isalpha()).upper()}"
def cit(m):
    o=OA.get(m['doi'] or '',{})
    au=LAT(o.get('authors') or []) or [m['au']]
    a=', '.join(ini(x) for x in au[:6])+(' et al.' if len(au)>6 else '') if o.get('authors') else m['au']
    b=o.get('biblio') or {}; vp=''
    if b.get('volume'): vp=f" {b['volume']}"+(f":{b['first_page']}"+(f"-{b['last_page']}" if b.get('last_page') and b['last_page']!=b['first_page'] else '') if b.get('first_page') else '')
    t=re.sub(r'\s+',' ',re.sub('<[^>]+>','',o.get('title') or m['title'])).strip()
    return f"{a} ({m['year']}) {t}. {o.get('journal') or m['journal']}{vp}", (au[0].split()[-1] if o.get('authors') else m['au'].split(',')[0].strip().title())
def short(last,yr):
    base=f"{last} et al. ({yr})"
    for s in [base]+[base[:-1]+L+')' for L in 'bcdefghij']:
        if s not in taken: taken.add(s); return s
# placements rejected after the arm check (position on the other arm than the paper states)
REJ={'QPm.vt-2BL':'Xbarc200 maps to 2BS (47.6 Mb) but the QTL is on 2BL','Lr (PI 289824, 5BS)':'Xgwm443 maps to 5BL (437.7 Mb) but the gene is on 5BS',
     'Qkb.ksu-6BS.1':'Xwmc105/Xgwm88 map to 6BL (430 Mb) but the QTL is on 6BS','LrZH84':'Xgwm582 maps to distal 1BS (2.6 Mb) while a later paper places LrZH84 on 1BL'}
src={};reports=[];n=0;cnt=collections.Counter();rows=[]
for r in E:
    m=M[r['i']]
    if m['i']==150: m=dict(m,doi='10.1007/s00122-026-05182-0')
    key='P94_'+re.sub(r'[^a-z0-9]','',(m['doi'] or m['title']).split('/',1)[-1].lower())[:24]
    if key not in src:
        c,last=cit(m); src[key]={'citation':c,'url':('https://doi.org/'+m['doi']) if m['doi'] else '','short':short(last,m['year'])}
    P=[(mk,pos(mk,r['chr'])) for mk in r['mk']]; P=[(mk,p) for mk,p in P if p and r['name'] not in REJ]
    s=e=pk=None; g='D'; ref=r['ref'] or 'RefSeq v1.0'; how=''
    if P:
        g='C'; ref='RefSeq v1.0'
        xs=[p[1] for _,p in P]
        if len(P)>1 and max(xs)-min(xs)<=30 and max(xs)>min(xs):
            s,e=round(min(xs),4),round(max(xs),4); how=f"interval {', '.join(f'{mk} {p[1]:.2f} Mb' for mk,p in P)}"
        else:
            pk=round(P[0][1][1],4); how=f"{P[0][0]} at {pk:.2f} Mb"
        how=f" Placed in silico on RefSeq v1.0 from the markers named by the paper ({how}; positions from {P[0][1][3]})."
    elif r['name'] in REJ: how=f" In-silico placement rejected: {REJ[r['name']]}; arm only."
    elif r['mk']: how=f" Markers named ({', '.join(r['mk'])}) have no known physical position; arm only."
    n+=1;cnt[key]+=1
    rec={'id':f'PD6-{n:03d}','name':r['name'],'cls':r['cls'],'disease':r['dis'],'chr':r['chr'],'start':s,'end':e,'peak':pk,'ref':ref,
         'marker':'; '.join(r['mk']),'pve':r['pve'],'stat':'','source':key,'grade':g,
         'curation':'auto (paper full text, in-silico)' if g=='C' else 'auto (paper full text)','evidence':r['ev'],
         'note':(r['note']+how).strip()+' Read from the full-text PDF supplied by the owner on 25 Nov 2026 (Downloads/Paper_2026/Papers); not yet checked by hand.'}
    reports.append(rec); rows.append([rec['id'],m['au'],m['year'],m['doi'] or '',rec['name'],rec['disease'],rec['chr'],g,s or '',e or '',pk or '',rec['marker']])
json.dump({'reports':reports,'sources':src,'papers':{}},open('seed94_pre.json','w'),ensure_ascii=False)
with open('added_records.tsv','w') as f:
    f.write('id\tfirst_author\tyear\tdoi\tname\tdisease\tchr\tgrade\tstart\tend\tpeak\tmarkers\n')
    for x in rows: f.write('\t'.join(map(str,x))+'\n')
print(len(reports),len(src),collections.Counter(r['grade'] for r in reports))
