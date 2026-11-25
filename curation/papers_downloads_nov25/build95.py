"""Build web/seed_ad95.json from records95.py (Downloads/Paper_2026/Paper_2026_1 and Downloads/Paper_2026/Alien).
Positions written by the paper -> grade B (interval) or A (peak); otherwise in silico from named markers (grade C, look2.py);
otherwise arm only (grade D). Records without a chromosome, and papers whose DOI is already an atlas source, are dropped."""
import json,re,sys,collections
S='/tmp/claude-0/-home-claude/00583e43-a309-532b-9f1b-7957bcf4ea95/scratchpad'
sys.path.insert(0,S+'/p26'); from look2 import pos
exec(open('records95.py').read())
F=[o for o in json.load(open(S+'/p27/files.json')) if o['st']=='new']
OA=json.load(open(S+'/p27/oa27.json')); CR={int(k):v for k,v in json.load(open(S+'/p27/cr_nodoi.json')).items()}
A=json.load(open(S+'/p26/atlas_all.json'))
dois=set()
for v in A['sources'].values():
    for d in re.findall(r'10\.\d{4,9}/[^\s"<>]+',json.dumps(v).lower()): dois.add(d.rstrip('.,);'))
for f in ('seed_ad94.json',):
    for v in json.load(open('/home/claude/wheatnhr-atlas/web/'+f))['sources'].values(): dois.add(v['url'].lower()[16:])
taken={v.get('short') for v in A['sources'].values()}|{v['short'] for v in json.load(open('/home/claude/wheatnhr-atlas/web/seed_ad94.json'))['sources'].values()}
BAD={97,101,253,371}   # crossref matches that are wrong
MAN={71:('Leonova I, Börner A, Budashkina E, Kalinina N, Unger O, Röder M, Salina E (2004) Identification of microsatellite markers for a leaf rust resistance gene introgressed into common wheat from Triticum timopheevii. Plant Breeding 123:93-95','10.1046/j.0179-9541.2003.00906.x','Leonova',2004),
     97:('Niranjana M et al. (2017) Cytogenetic analysis and mapping of leaf rust resistance in Aegilops speltoides Tausch derived bread wheat line "Selection 2427" carrying putative gametocidal gene(s). Genome','','Niranjana',2017),
     107:('Rong JK, Millet E, Manisterski J, Feldman M (2000) A new powdery mildew resistance gene: introgression from wild emmer into common wheat and RFLP-based mapping. Euphytica 115:121-126','10.1023/a:1003950431049','Rong',2000),
     253:('Hou L, Jia J, Zhang X, Li X et al. (2016) Molecular mapping of the stripe rust resistance gene Yr69 on wheat chromosome 2AS. Plant Disease','','Hou',2016),
     339:("Singh D, Park RF, Bariana HS, McIntosh RA (2001) Cytogenetic studies in wheat XIX. Chromosome location and linkage studies of a gene for leaf rust resistance in the Australian cultivar 'Harrier'. Plant Breeding 120:7-12",'10.1046/j.1439-0523.2001.00561.x','Singh',2001),
     371:("Wu Y, Wang Y, Yao F, Long L, Li J, Li H et al. (2021) Molecular mapping of a novel QTL conferring adult plant resistance to stripe rust in Chinese wheat landrace 'Guangtoumai'. Plant Disease",'','Wu',2021)}
LAT=lambda L:[x for x in L if all(ord(c)<0x2E80 for c in x)]
def ini(n):
    p=[x.title() if x.isupper() and len(x)>1 else x for x in n.replace('.',' ').replace('-',' ').split()]; return f"{p[-1]} {''.join(x[0] for x in p[:-1] if x[0].isalpha()).upper()}"
def clean(t): return re.sub(r'\s+',' ',re.sub('<[^>]+>','',t or '')).replace('&lt;','').replace('&gt;','').strip()
def meta(k):
    if k in MAN: c,d,last,yr=MAN[k]; return c,d,last,yr
    o=F[k]; d=o.get('doi')
    if d and d in OA and 'err' not in OA[d] and OA[d].get('title'):
        m=OA[d]; au=LAT(m['authors']) or ['Anon']; b=m.get('biblio') or {}; yr=m['year']
    elif k in CR:
        it=CR[k]; d=it['DOI'].lower(); au=[f"{a.get('given','')} {a.get('family','')}".strip() for a in it.get('author',[])] or ['Anon']
        m={'title':it['title'][0],'journal':(it.get('container-title') or [''])[0]}; yr=it['issued']['date-parts'][0][0]
        b={'volume':it.get('volume'),'first_page':(it.get('page') or '').split('-')[0],'last_page':(it.get('page') or '').split('-')[-1]}
    else: return None
    a=', '.join(ini(x) for x in au[:6])+(' et al.' if len(au)>6 else '')
    vp=''
    if b.get('volume'): vp=f" {b['volume']}"+(f":{b['first_page']}"+(f"-{b['last_page']}" if b.get('last_page') and b['last_page']!=b['first_page'] else '') if b.get('first_page') else '')
    au=[(' '.join(x.title() if x.isupper() and len(x)>1 else x for x in y.split())) for y in au]
    return f"{a} ({yr}) {clean(m['title'])}. {clean(m.get('journal') or '')}{vp}".rstrip('. ')+('' if vp else ''), d, au[0].split()[-1], yr
def short(last,yr):
    base=f"{last} et al. ({yr})"
    for s in [base]+[base[:-1]+L+')' for L in 'bcdefghijklm']:
        if s not in taken: taken.add(s); return s
src={};reports=[];rows=[];drop=collections.Counter();n=0;used=collections.Counter()
for r in E:
    k=r['k']
    if not r['chr']: drop['no chromosome']+=1; continue
    if k==101: drop['student essay']+=1; continue
    mt=meta(k)
    if not mt: drop['no metadata']+=1; continue
    cit,doi,last,yr=mt
    if doi and doi.lower() in dois: drop['paper already in atlas']+=1; continue
    key='P95_'+re.sub(r'[^a-z0-9]','',(doi or F[k]['f']).split('/',1)[-1].lower())[:24]
    if key not in src: src[key]={'citation':cit,'url':('https://doi.org/'+doi) if doi else '','short':short(last,yr)}
    s=e=p=None; g='D'; ref=r['ref'] or 'RefSeq v1.0'; how=''
    if r['pref']:
        ref=r['pref']; s,e,p=r['s'],r['en'],r['p']; g='B' if s is not None else 'A'
        if s is not None and p is not None: g='B'
    else:
        P=[(mk,pos(mk,r['chr'])) for mk in r['mk']]; P=[(mk,q) for mk,q in P if q]
        if P:
            g='C'; ref='RefSeq v1.0'; xs=[q[1] for _,q in P]
            if len(P)>1 and 0<max(xs)-min(xs)<=30:
                s,e=round(min(xs),4),round(max(xs),4); how=f"interval {', '.join(f'{mk} {q[1]:.2f} Mb' for mk,q in P)}"
            else:
                p=round(P[0][1][1],4); how=f"{P[0][0]} at {p:.2f} Mb"
            how=f" Placed in silico on RefSeq v1.0 from the markers named by the paper ({how}; positions from {P[0][1][3]})."
        elif r['mk']: how=f" Markers named ({', '.join(r['mk'])}) have no known physical position; arm only."
    n+=1; used[key]+=1
    folder='Alien' if F[k]['f'].startswith('A__') else 'Paper_2026_1'
    rec={'id':f'PE6-{n:03d}','name':r['name'],'cls':r['cls'],'disease':r['dis'],'chr':r['chr'],'start':s,'end':e,'peak':p,'ref':ref,
         'marker':'; '.join(r['mk']),'pve':r['pve'],'stat':'','source':key,'grade':g,
         'curation':{'B':'auto (paper full text)','A':'auto (paper full text)','C':'auto (paper full text, in-silico)','D':'auto (paper full text)'}[g],
         'evidence':r['ev'],'note':(r['note']+how).strip()+f' Read from the full-text PDF supplied by the owner on 25 Nov 2026 (Downloads/Paper_2026/{folder}); not yet checked by hand.'}
    reports.append(rec); rows.append([rec['id'],F[k]['f'][3:],doi or '',rec['name'],rec['disease'],rec['chr'],g,s if s is not None else '',e if e is not None else '',p if p is not None else '',ref])
json.dump({'reports':reports,'sources':src,'papers':{}},open('seed95_pre.json','w'),ensure_ascii=False)
with open('added_records95.tsv','w') as f:
    f.write('id\tfile\tdoi\tname\tdisease\tchr\tgrade\tstart\tend\tpeak\tref\n')
    for x in rows: f.write('\t'.join(map(str,x))+'\n')
# tracker of every file in the two folders
recs=collections.Counter(r['k'] for r in E if r['chr'])
allF=json.load(open(S+'/p27/files.json'))
with open('tracker95.tsv','w') as f:
    f.write('folder\tfile\tdoi\tstatus\trecords_written\n')
    newi={o['f']:i for i,o in enumerate(F)}
    for o in allF:
        st=o['st']; i=newi.get(o['f'])
        if st=='new': st='loci added' if recs.get(i) else 'no locus with a wheat chromosome extracted'
        f.write('\t'.join([('Alien' if o['f'].startswith('A__') else 'Paper_2026_1'),o['f'][3:],o.get('doi') or '',st,str(recs.get(i,0))])+'\n')
print(len(reports),len(src),collections.Counter(r['grade'] for r in reports),dict(drop))
