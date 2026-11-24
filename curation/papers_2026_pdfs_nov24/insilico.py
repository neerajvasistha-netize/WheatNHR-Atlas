import re,sys,json,collections
sys.path.insert(0,'/tmp/claude-0/-home-claude/00583e43-a309-532b-9f1b-7957bcf4ea95/scratchpad/pp')
from mp import look
def pos(m,ch):
    L=[x for x in look(m,ch) if x[2] in ('RefSeq v1.0','RefSeq v2.1')]
    if not L and m.endswith('_dup2'): return pos(m[:-5],ch)
    v1=[x for x in L if x[2]=='RefSeq v1.0']; return (v1 or L or [None])[0]
out={}
# Thirunarayanan spot blotch: Tables 4 (AUDPC) and 5 (IP)
t=open('t/1-s2.0-S0885576526003450-main.txt').read().split('\n')
trait=None
for l in t:
    if l.startswith('Table 4'): trait='AUDPC'
    elif l.startswith('Table 5'): trait='IP'
    elif re.match(r'^Table \d',l): trait=None
    m=re.match(r'\s+(\d[ABD]\.\d+)\s+(\S+)\s+(\d[ABD])\s+(\d+)\s+(\S+)\s+(\d{4})\s+(\S+)',l)
    if trait and m:
        loc,mk,ch,cm,iso,yr,site=m.groups()
        k=('SB',mk,ch); d=out.setdefault(k,{'traits':set(),'env':set(),'cm':cm})
        d['traits'].add(trait); d['env'].add(f'{iso} {yr} {site}')
# Korra Exserohilum: Table 3
t=open('t/plants-15-02351.txt').read().split('\n')
for l in t:
    m=re.match(r'\s+(\S+)\s+(\d[ABD])\s+(\d+)\s+(AUDPC|IP)\s+(BHU|UAS)\s+(\d{4})\s+(TraesCS\w+)',l)
    if m:
        mk,ch,cm,tr,site,yr,g=m.groups(); k=('ES',mk,ch); d=out.setdefault(k,{'traits':set(),'env':set(),'cm':cm,'gene':g})
        d['traits'].add(tr); d['env'].add(f'{site} {yr}')
res=[]
for (dz,mk,ch),d in out.items():
    p=pos(mk,ch); res.append(dict(dz=dz,marker=mk,chr=ch,cm=d['cm'],traits=sorted(d['traits']),env=sorted(d['env']),gene=d.get('gene'),pos=p))
json.dump(res,open('insilico.json','w'),indent=0)
c=collections.Counter((r['dz'],r['pos'] is not None) for r in res);print(c)
